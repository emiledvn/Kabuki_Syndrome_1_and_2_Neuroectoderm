#!/usr/bin/env Rscript
# A09_Motif_Enrichment.R -- TF motif enrichment (monaLisa::calcBinnedMotifEnrR),
# separately for GAINED and LOST Loess DARs per contrast, each tested against
# non-significant peaks as background. Complements A07b's TOBIAS BINDetect
# footprinting layer with a genuinely different, sequence-level question: "is
# this motif's sequence overrepresented in peaks that opened/closed" (this
# script) vs "does measured footprint depth differ at this motif's known
# binding sites between genotypes" (TOBIAS). Same motif database (motifs.meme)
# for both, so results are directly comparable.
#
# monaLisa was chosen over lighter alternatives (universalmotif alone,
# motifmatchr) specifically because calcBinnedMotifEnrR() automatically
# reweights background sequences for GC and k-mer composition differences
# (same approach as Homer findMotifsGenome.pl, see ?calcBinnedMotifEnrR
# Details) -- ATAC peak GC content correlates with motif content, so an
# unmatched background risks spurious "enrichment" that is really a GC
# artifact. monaLisa is not installable via bioconda's default channel
# priority (a transitive dependency, r-tfmpvalue, only has bioconda builds
# pinned to R 3.3/3.4) -- installed 2026-07-20 via
# `mamba install ... --channel-priority flexible`, which picks up
# r-tfmpvalue=1.0.0 from conda-forge instead. universalmotif was added
# alongside it purely as a MEME-format bridge: monaLisa has no MEME reader of
# its own (only TFBSTools/Homer input), so motifs.meme is parsed with
# universalmotif::read_meme() and converted to a TFBSTools PWMatrixList.
# See envs/ks_1_2_r.yml (re-exported after both installs).
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/A09_Motif_Enrichment.R
# Requires: A04b_Csaw_Loess_Norm.R already run (reads its saved full results
# table directly -- A04f's BED files are not needed here, this script derives
# GAINED/LOST/Background itself from the same padj/log2FC/lfc_threshold columns
# A04f uses, so the two stay consistent by construction).
#
# Self-checkpointing: skips entirely if results/tables/A09_summary.tsv exists.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble)
  library(GenomicRanges); library(Rsamtools); library(Biostrings)
  library(TFBSTools); library(universalmotif)
  library(monaLisa)
  library(SummarizedExperiment)
  library(BiocParallel)
  library(ggplot2)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[A09] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/atac"
SESSION_DIR <- "results/Session_info"
SUMMARY     <- file.path(TABLES_DIR, "A09_summary.tsv")
for (d in c(TABLES_DIR, PLOTS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(SUMMARY)) {
  cat("[A09] Already complete (", SUMMARY, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}

FDR       <- cfg$thresholds$fdr
LFC_THR   <- cfg$thresholds$lfc_threshold
CONTRASTS <- cfg$contrasts
CORES     <- cfg$diffbind$cores  # R/BiocParallel-heavy in-memory work, same class of job as A04/A04b -- see pipeline_config.yaml's diffbind.cores comment

# Same "decompress alongside the .gz" convention A06_TOBIAS.sh already uses
# for the reference FASTA (MOODS/samtools-family tools need an uncompressed
# or bgzip-indexed file, not plain gzip); already exists with a .fai here.
FASTA       <- sub("\\.gz$", "", cfg$reference$fasta)
MOTIFS_MEME <- cfg$reference$motifs_meme
if (!file.exists(FASTA))                    stop("[A09] ERROR: ", FASTA, " not found (expected decompressed alongside the .gz -- see A06_TOBIAS.sh for the same convention).")
if (!file.exists(paste0(FASTA, ".fai")))     stop("[A09] ERROR: ", FASTA, ".fai not found -- run samtools faidx on the genome first.")
if (!file.exists(MOTIFS_MEME))               stop("[A09] ERROR: ", MOTIFS_MEME, " not found.")

cat("================================================================================\n")
cat("A09_Motif_Enrichment -- monaLisa, GAINED/LOST Loess DARs vs non-significant background\n")
cat("================================================================================\n\n")

cat("[A09] Loading motifs from ", MOTIFS_MEME, "...\n", sep = "")
motifs_um <- universalmotif::read_meme(MOTIFS_MEME)
# universalmotif's convert_motifs(..., class="TFBSTools-PWMatrix") maps its own
# "name" slot (here: JASPAR accession, e.g. "MA0139.1", guaranteed unique per
# motif) to TFBSTools name(), and "altname" (here: TF gene symbol, e.g. "CTCF")
# to TFBSTools ID() -- the reverse of what monaLisa needs internally.
# calcBinnedMotifEnrR()'s hit-matrix construction requires ID() to be unique
# (crashes with "factor level [n] is duplicated" otherwise), and 51/879 of
# these motifs share a TF name across multiple JASPAR entries (CTCF x3,
# TFAP2A/B/C x3 each, etc.), so ID() as converted is not unique. TFBSTools
# doesn't export an ID<- replacement method to fix this after conversion (only
# the ID() getter), so swap universalmotif's own name/altname slots BEFORE
# conversion instead: after the swap, the guaranteed-unique accession lands in
# TFBSTools ID() and the human-readable TF name lands in name() -- safe for
# name() to have duplicates (display-only, not an internal key).
motifs_um <- lapply(motifs_um, function(mot) {
  acc <- mot["name"]; tfname <- mot["altname"]
  mot["name"] <- tfname
  mot["altname"] <- acc
  mot
})
pwmL <- do.call(TFBSTools::PWMatrixList, universalmotif::convert_motifs(motifs_um, class = "TFBSTools-PWMatrix"))
stopifnot(!any(duplicated(TFBSTools::ID(pwmL))))
cat(sprintf("[A09] %d motifs loaded.\n\n", length(pwmL)))

########
## Restrict to motifs whose TF is actually expressed in this system -- same
## min_mean_expression threshold A07b already uses for exactly this purpose
## ("TF/gene must average > this many normalized counts to count as
## expressed"), applied here for consistency rather than inventing a new one.
## An enriched motif for an unexpressed TF is most likely picking up a
## similarly-shaped motif from an expressed paralog/family member, not real
## signal from that specific factor -- and testing it anyway only costs
## multiple-testing power for nothing.
##
## Composite/dimer motifs (name() holds the TF, e.g. "RXRA::VDR", 76/879 of
## these) require ALL partners expressed, since a functional heterodimer needs
## both proteins present. Matching is case-insensitive because JASPAR mixes
## species capitalization conventions in TF names (e.g. "Arnt", "Pou5f1" are
## mouse-style orthologs of the same human genes we do have expression data
## for -- ARNT, POU5F1); a motif whose name has no match at all in the human
## gene set (e.g. "EWSR1-FLI1", an oncogenic fusion protein, not a normal gene
## product) is correctly dropped as unexpressed/not applicable here.
########

MIN_MEAN_EXPR <- cfg$thresholds$min_mean_expression
norm_counts <- readRDS("results/RDS/R02_normalized_counts.rds")
gene_map    <- readRDS("results/RDS/R02_gene_id_name_map.rds")
mean_expr_by_symbol <- tibble(gene_id = rownames(norm_counts), mean_expr = rowMeans(norm_counts)) %>%
  left_join(gene_map, by = "gene_id") %>%
  filter(!is.na(gene_name)) %>%
  group_by(gene_name) %>% summarise(mean_expr = max(mean_expr), .groups = "drop")  # multi-mapped IDs: max, not sum -- avoid inflating a duplicate-mapped gene's apparent expression
expressed_symbols_upper <- toupper(mean_expr_by_symbol$gene_name[mean_expr_by_symbol$mean_expr > MIN_MEAN_EXPR])

tf_is_expressed <- function(tf_name) {
  parts <- toupper(strsplit(tf_name, "::", fixed = TRUE)[[1]])
  all(parts %in% expressed_symbols_upper)
}
motif_expressed <- vapply(TFBSTools::name(pwmL), tf_is_expressed, logical(1))
cat(sprintf("[A09] TF-expression filter (mean normalized count > %d, matching A07b's threshold): %d/%d motifs kept.\n",
            MIN_MEAN_EXPR, sum(motif_expressed), length(pwmL)))
pwmL <- pwmL[motif_expressed]
cat(sprintf("[A09] %d motifs after expression filter.\n\n", length(pwmL)))

fa <- Rsamtools::FaFile(FASTA)
fa_seqlengths <- GenomeInfoDb::seqlengths(Rsamtools::seqinfo(fa))
BP <- BiocParallel::MulticoreParam(workers = CORES)

peak_id_to_gr <- function(peak_id) {
  m <- regmatches(peak_id, regexec("^(.+):(\\d+)-(\\d+)$", peak_id))
  gr <- GRanges(seqnames = sapply(m, `[`, 2),
                ranges = IRanges(start = as.numeric(sapply(m, `[`, 3)), end = as.numeric(sapply(m, `[`, 4))))
  mcols(gr)$peak_id <- peak_id
  gr
}

run_one_contrast <- function(ct) {
  cat(sprintf("[A09] === %s ===\n", ct$name))
  full_path <- file.path(TABLES_DIR, sprintf("A04b_%s_CsawLoess_full.csv", ct$name))
  if (!file.exists(full_path)) { cat("[A09] Skipping -- ", full_path, " not found.\n"); return(NULL) }
  res <- read_csv(full_path, show_col_types = FALSE)

  peaks_gr <- peak_id_to_gr(res$peak_id)
  # Defensive: a handful of consensus peaks can extend past the actual end of
  # a short/unplaced scaffold (e.g. observed: a peak at KI270336.1:664-1064 on
  # a 1026bp contig) -- Rsamtools::getSeq() errors on these ("record was
  # truncated"). The pipeline already has a min_edge_distance config threshold
  # meant to drop contig-edge peaks (A02/A06), but it evidently wasn't applied
  # upstream of this consensus set; not extractable and not biologically
  # meaningful at these coordinates either way, so drop them here rather than
  # touching upstream peak-calling.
  contig_len <- fa_seqlengths[as.character(seqnames(peaks_gr))]
  valid <- !is.na(contig_len) & end(peaks_gr) <= contig_len
  if (any(!valid)) {
    cat(sprintf("[A09] Dropping %d/%d peaks that extend past their contig end (not extractable): %s\n",
                sum(!valid), length(peaks_gr), paste(head(res$peak_id[!valid], 5), collapse = ", ")))
  }
  res      <- res[valid, ]
  peaks_gr <- peaks_gr[valid]

  direction <- case_when(
    !is.na(res$padj) & res$padj < FDR & res$log2FoldChange >  LFC_THR ~ "GAINED",
    !is.na(res$padj) & res$padj < FDR & res$log2FoldChange < -LFC_THR ~ "LOST",
    TRUE ~ "Background"
  )
  cat(sprintf("[A09] %d GAINED, %d LOST, %d Background peaks\n",
              sum(direction == "GAINED"), sum(direction == "LOST"), sum(direction == "Background")))

  cat("[A09] Extracting peak sequences from genome FASTA...\n")
  seqs <- Biostrings::getSeq(fa, peaks_gr)
  names(seqs) <- mcols(peaks_gr)$peak_id

  results_list <- list()
  for (dir_of_interest in c("GAINED", "LOST")) {
    # Test only against true non-significant background peaks -- EXCLUDE the
    # opposite direction entirely, don't relabel it into "Background". The
    # prior version did `ifelse(direction == dir_of_interest, dir_of_interest,
    # "Background")`, which folded the opposite direction's peaks into the
    # "background" pool too (e.g. LOST peaks counted as background when
    # testing GAINED). Any motif genuinely enriched/depleted in the opposite
    # direction then contaminated the background rate used for this
    # direction's test, biasing log2enr in a motif-dependent way (not a fixed
    # direction) rather than the "peaks with no significant accessibility
    # change" background the header comment describes.
    keep <- direction %in% c(dir_of_interest, "Background")
    bins <- factor(direction[keep], levels = c("Background", dir_of_interest))
    cat(sprintf("[A09] %s vs Background: running calcBinnedMotifEnrR (%d cores, this can take a while)...\n",
                dir_of_interest, CORES))
    se <- calcBinnedMotifEnrR(seqs = seqs[keep], bins = bins, pwmL = pwmL, background = "otherBins",
                               test = "fisher", BPPARAM = BP, verbose = TRUE)

    df <- tibble(
      motif_id      = rowData(se)$motif.id,
      motif_name    = rowData(se)$motif.name,
      motif_pctGC   = rowData(se)$motif.percentGC,
      negLog10P     = assay(se, "negLog10P")[, dir_of_interest],
      negLog10Padj  = assay(se, "negLog10Padj")[, dir_of_interest],
      log2enr       = assay(se, "log2enr")[, dir_of_interest],
      n_fg_with_hit = assay(se, "sumForegroundWgtWithHits")[, dir_of_interest],
      n_bg_with_hit = assay(se, "sumBackgroundWgtWithHits")[, dir_of_interest]
    ) %>% arrange(desc(negLog10Padj))

    write_csv(df, file.path(TABLES_DIR, sprintf("A09_%s_%s_vs_Background.csv", ct$name, dir_of_interest)))
    n_sig <- sum(df$negLog10Padj > -log10(FDR), na.rm = TRUE)
    cat(sprintf("[A09] %s %s: %d motifs significant (FDR < %.2f) of %d tested\n",
                ct$name, dir_of_interest, n_sig, FDR, nrow(df)))

    results_list[[dir_of_interest]] <- df %>% mutate(Contrast = ct$name, Direction = dir_of_interest)

    if (n_sig > 0) {
      # Selected AND ordered by |log2enr| (effect size), not negLog10Padj --
      # with tens of thousands of peaks per bin, statistical power is high
      # enough that ranking by significance alone surfaces lots of small,
      # merely-detectable effects rather than the strongest ones (the same
      # redundancy/power issue flagged for the GO BP GSEA results).
      top <- df %>% filter(negLog10Padj > -log10(FDR)) %>% arrange(desc(abs(log2enr))) %>% head(20)
      p <- ggplot(top, aes(x = log2enr, y = reorder(motif_name, abs(log2enr)), colour = negLog10Padj)) +
        geom_point(size = 3) +
        scale_colour_gradientn(colours = c("#2166ac", "#4393c3", "#d6604d", "#b2182b"), name = expression(-log[10](FDR))) +
        labs(title = sprintf("Motif enrichment: %s, %s vs background", ct$name, dir_of_interest),
             subtitle = "Top 20 by |log2 enrichment| among FDR-significant motifs",
             x = "log2 enrichment", y = NULL) +
        theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
      ggsave(file.path(PLOTS_DIR, sprintf("A09_%s_%s_dotplot.pdf", ct$name, dir_of_interest)), p, width = 8, height = 7)
    }
  }
  bind_rows(results_list)
}

all_results <- bind_rows(lapply(CONTRASTS, run_one_contrast))

summary_df <- all_results %>% group_by(Contrast, Direction) %>%
  summarise(N_tested = n(), N_significant = sum(negLog10Padj > -log10(FDR), na.rm = TRUE), .groups = "drop")
write_csv(summary_df, SUMMARY)
cat("\n[A09] Summary:\n"); print(summary_df)

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "A09_Motif_Enrichment_session_info.txt"))
cat("\n[DONE] A09_Motif_Enrichment complete\n")
