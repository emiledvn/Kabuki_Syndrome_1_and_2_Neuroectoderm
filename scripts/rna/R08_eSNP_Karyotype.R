#!/usr/bin/env Rscript
# R08_eSNP_Karyotype.R -- RNA-seq-based aneuploidy/karyotype QC screen for the iPSC
# lines, using the two methods described to the user for exactly this purpose:
#
#   1. eSNP-Karyotyping (Weissbein, Benvenisty & Ben-David, Stem Cell Reports 2016):
#      bin heterozygous-SNV B-allele frequency (BAF = alt/(ref+alt) read depth) by
#      chromosome arm. A systematic shift away from ~0.5, consistent across
#      replicates of one genotype and absent in WT, is the RNA-seq signature of a
#      large-scale copy-number change (trisomy/monosomy/segmental gain-loss).
#   2. Companion dosage signal (Ben-David et al., Cell Stem Cell 2014): aneuploid
#      regions also show elevated mean expression (gene dosage effect) -- computed
#      here for free from R02's already-normalized VST counts, as a second,
#      independent line of evidence alongside the BAF panel.
#
# This is a QC screen, not a diagnostic test: both methods are threshold-based in
# their own original design, not p-value-based -- flags below are soft flags with
# their supporting numbers reported alongside, same as the rest of this pipeline's
# QC-only checks (see AR10b's header for the same philosophy).
#
# Requires R07_RUN_nfcore_RNAVAR.sh (per-sample filtered VCFs) and R02_DESeq2.R
# (VST counts, gene_id/gene_name map) to have completed first.
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/R08_eSNP_Karyotype.R
#
# Self-checkpointing: skips entirely if results/eSNP_Karyotype_SUMMARY.txt exists.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(stringr)
  library(purrr)
  library(ggplot2)
  library(pheatmap)
  library(DESeq2)   # only for assay()/colData() generics on the saved vsd object
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[R08] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

TABLES_DIR  <- "results/tables"
PLOTS_DIR   <- "results/rna"
RDS_DIR     <- "results/RDS"
SESSION_DIR <- "results/Session_info"
RNAVAR_DIR  <- "data/RNAvar_nfcore_output"
SUMMARY     <- "results/eSNP_Karyotype_SUMMARY.txt"

for (d in c(TABLES_DIR, PLOTS_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(SUMMARY)) {
  cat("[R08] Already complete (", SUMMARY, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}
if (!dir.exists(RNAVAR_DIR)) stop("[R08] ERROR: ", RNAVAR_DIR, " not found -- run R07_RUN_nfcore_RNAVAR.sh first.")
if (!file.exists(file.path(RDS_DIR, "R02_vst_counts.rds"))) stop("[R08] ERROR: missing R02 output -- run R02_DESeq2.R first.")

ESNP_MIN_DEPTH       <- cfg$thresholds$esnp_min_depth
ESNP_MIN_SNPS_PER_ARM <- cfg$thresholds$esnp_min_snps_per_arm
ESNP_BAF_FLAG         <- cfg$thresholds$esnp_baf_deviation_flag
ESNP_MIN_GENES_PER_ARM <- cfg$thresholds$esnp_min_genes_per_arm
FASTA_FAI             <- paste0(sub("\\.gz$", "", cfg$reference$fasta), ".fai")
CHROM_ARMS_BED        <- cfg$reference$chrom_arms_bed
CYTOBAND_URL          <- cfg$reference$cytoband_url
CONTRASTS             <- cfg$contrasts
GENOTYPE_LEVELS        <- c("WT", CONTRASTS[[1]]$treatment, CONTRASTS[[2]]$treatment)
GENOTYPE_COLOR         <- setNames(c("grey40", "#1B9E77", "#7570B3"), GENOTYPE_LEVELS)

cat("================================================================================\n")
cat("R08_eSNP_Karyotype -- RNA-seq aneuploidy/karyotype QC screen\n")
cat("================================================================================\n\n")

########
## 1. Chromosome-arm reference (build once, cache) -- centromere = midpoint of the
## two UCSC "acen"-stain bands per chromosome; arm boundaries = [0, centromere) /
## [centromere, chrom_length]. Chrom lengths from our own pinned FASTA's .fai
## (already on disk, no extra download). Same download-derive-discard-raw idiom as
## A03_QC_ATAC.sh's RefSeq TSS bed. Restricted to autosomes + X: chrY excluded --
## RNA-seq coverage there is sparse/absent and sex-dependent, not informative for a
## heterozygous-SNV BAF scan; chrM excluded (no arms, circular, not nuclear).
########

if (!file.exists(CHROM_ARMS_BED)) {
  cat("[R08] Deriving chromosome-arm reference (one-time) from UCSC cytoBand...\n")
  raw <- tempfile(fileext = ".txt.gz")
  download.file(CYTOBAND_URL, raw, quiet = TRUE, mode = "wb")
  cyto <- read_tsv(raw, col_names = c("CHROM", "start", "end", "band", "stain"), show_col_types = FALSE) %>%
    filter(str_detect(CHROM, "^chr([1-9]|1[0-9]|2[0-2]|X)$"))
  centromere <- cyto %>% filter(stain == "acen") %>%
    group_by(CHROM) %>%
    summarise(centromere = as.integer(round((min(start) + max(end)) / 2)), .groups = "drop")
  fai <- read_tsv(FASTA_FAI, col_names = c("CHROM", "length", "offset", "linebases", "linewidth"), show_col_types = FALSE) %>%
    dplyr::select(CHROM, length) %>%
    filter(CHROM %in% centromere$CHROM)
  cent_len <- inner_join(centromere, fai, by = "CHROM")
  arm_bed <- bind_rows(
    cent_len %>% transmute(CHROM, arm = paste0(CHROM, "_p"), start = 0L, end = centromere),
    cent_len %>% transmute(CHROM, arm = paste0(CHROM, "_q"), start = centromere, end = length)
  ) %>% arrange(CHROM, arm)
  write_tsv(arm_bed, CHROM_ARMS_BED)
  unlink(raw)
  cat(sprintf("[R08] Wrote %s: %d arms across %d chromosomes\n\n", CHROM_ARMS_BED, nrow(arm_bed), length(unique(arm_bed$CHROM))))
}
arm_bed    <- read_tsv(CHROM_ARMS_BED, show_col_types = FALSE)
centromere <- arm_bed %>% filter(str_ends(arm, "_p")) %>% dplyr::select(CHROM, centromere = end)
arm_order  <- arm_bed$arm

assign_arm <- function(chrom_vec, pos_vec) {
  cent <- centromere$centromere[match(chrom_vec, centromere$CHROM)]
  ifelse(is.na(cent), NA_character_, paste0(chrom_vec, if_else(pos_vec <= cent, "_p", "_q")))
}

########
## 2. Locate + parse each RNA sample's rnavar filtered VCF. Plain-text TSV parsing
## (readr::read_tsv, comment = "##") rather than adding VariantAnnotation/vcfR/
## bcftools -- keeps this dependency-free at the filtered-VCF scale we need, and
## avoids the documented slow-solver risk of a fresh conda install into ks_1_2_r.
########

meta <- read_csv("data/sample_metadata.csv", show_col_types = FALSE) %>% filter(assay == "RNA")
cat(sprintf("[R08] %d RNA samples: %s\n", nrow(meta), paste(meta$sample_id, collapse = ", ")))

find_sample_vcf <- function(sample_id) {
  candidates <- list.files(RNAVAR_DIR, pattern = paste0("^", sample_id, ".*\\.vcf\\.gz$"),
                            recursive = TRUE, full.names = TRUE)
  candidates <- candidates[!str_detect(candidates, "\\.g\\.vcf\\.gz$")]
  if (length(candidates) == 1) return(candidates)
  if (length(candidates) > 1) {
    filtered <- candidates[str_detect(candidates, "filter")]
    if (length(filtered) == 1) return(filtered)
  }
  stop(sprintf("[R08] ERROR: expected exactly 1 VCF for %s under %s, found %d:\n  %s",
               sample_id, RNAVAR_DIR, length(candidates), paste(candidates, collapse = "\n  ")))
}

parse_sample_vcf <- function(vcf_path, sample_id) {
  df <- suppressWarnings(read_tsv(vcf_path, comment = "##", show_col_types = FALSE))
  names(df)[names(df) == "#CHROM"] <- "CHROM"
  stopifnot("sample column missing from VCF -- rnavar sample naming may not match sample_metadata.csv sample_id" = sample_id %in% names(df))

  het <- df %>%
    filter(FILTER == "PASS", nchar(REF) == 1, nchar(ALT) == 1, !str_detect(ALT, ",")) %>%
    transmute(CHROM, POS, FORMAT, GTFIELD = .data[[sample_id]])
  if (nrow(het) == 0) return(tibble(sample_id = character(), CHROM = character(), POS = integer(),
                                     ref_depth = integer(), alt_depth = integer(), depth = integer(), baf = double()))

  split_one_format <- function(fmt, gtfield) {
    keys <- str_split(fmt, ":")[[1]]
    parts <- str_split(gtfield, ":")
    gt_i <- match("GT", keys); ad_i <- match("AD", keys)
    tibble(GT = vapply(parts, `[`, character(1), gt_i),
           AD = vapply(parts, `[`, character(1), ad_i))
  }

  het %>%
    group_by(FORMAT) %>%
    group_modify(~ bind_cols(.x, split_one_format(.y$FORMAT[1], .x$GTFIELD))) %>%
    ungroup() %>%
    filter(GT %in% c("0/1", "0|1")) %>%
    separate(AD, into = c("ref_depth", "alt_depth"), sep = ",", convert = TRUE, extra = "drop") %>%
    mutate(depth = ref_depth + alt_depth, baf = alt_depth / depth, sample_id = sample_id) %>%
    filter(depth >= ESNP_MIN_DEPTH) %>%
    dplyr::select(sample_id, CHROM, POS, ref_depth, alt_depth, depth, baf)
}

cat("[R08] Parsing per-sample VCFs...\n")
snvs <- map_dfr(meta$sample_id, function(sid) {
  vcf <- find_sample_vcf(sid)
  d <- parse_sample_vcf(vcf, sid)
  cat(sprintf("  %s: %d PASS het biallelic SNVs (depth >= %d) from %s\n", sid, nrow(d), ESNP_MIN_DEPTH, basename(vcf)))
  d
}) %>%
  mutate(arm = assign_arm(CHROM, POS)) %>%
  filter(!is.na(arm)) %>%
  left_join(meta %>% dplyr::select(sample_id, genotype), by = "sample_id") %>%
  mutate(genotype = factor(genotype, levels = GENOTYPE_LEVELS))

########
## 3. Aggregate per sample x arm, flag candidate arms
########

baf_by_arm <- snvs %>%
  group_by(sample_id, genotype, CHROM, arm) %>%
  summarise(n_snps = n(), median_baf = median(baf), mean_abs_dev = mean(abs(baf - 0.5)), .groups = "drop") %>%
  mutate(arm = factor(arm, levels = arm_order))
write_csv(baf_by_arm, file.path(TABLES_DIR, "R08_BAF_by_arm.csv"))

informative <- baf_by_arm %>% filter(n_snps >= ESNP_MIN_SNPS_PER_ARM)
per_genotype_flag <- informative %>%
  group_by(genotype, arm) %>%
  summarise(n_reps = n(), all_reps_flagged = n() == 3 & all(mean_abs_dev > ESNP_BAF_FLAG),
            mean_of_means = mean(mean_abs_dev), .groups = "drop") %>%
  filter(all_reps_flagged)
wt_flagged_arms <- per_genotype_flag %>% filter(genotype == "WT") %>% pull(arm)
flagged_arms <- per_genotype_flag %>% filter(genotype != "WT", !arm %in% wt_flagged_arms)
write_csv(flagged_arms, file.path(TABLES_DIR, "R08_flagged_arms.csv"))

cat(sprintf("\n[R08] Candidate arms flagged (all 3 reps, mean|BAF-0.5| > %.2f, not also in WT): %d\n", ESNP_BAF_FLAG, nrow(flagged_arms)))
if (nrow(flagged_arms) > 0) print(flagged_arms)

########
## 4. Companion dosage panel (Ben-David et al. 2014 signal): mean-center each
## GENE across samples (subtract that gene's own mean VST across all 9 samples,
## i.e. row-wise on the genes x samples matrix), then average per sample x arm.
## This isolates how much a given sample deviates from the cohort at those
## genes -- the actual dosage-effect signal. Centering the other way (each
## SAMPLE's value against its own genome-wide mean, subtracting colMeans) was
## tried first and produces an arm-level pattern driven almost entirely by
## which arms are gene-dense/highly-expressed in general (chr19, chr12q, etc.
## always red; sparse acrocentric arms always blue) -- IDENTICAL in every
## sample regardless of genotype, since it never removes each gene's own
## typical expression level from the picture. That version is not informative
## for aneuploidy at all and was corrected before this was ever reported.
########

cat("\n[R08] Building gene -> chromosome-arm map from GTF (same idiom as R02's gene_id map)...\n")
GTF <- cfg$reference$gtf
gtf_lines <- system(sprintf("zcat %s | awk -F'\t' '$3==\"gene\"{print $1\"\\t\"$4\"\\t\"$9}'", GTF), intern = TRUE)
gene_pos <- tibble(raw = gtf_lines) %>%
  separate(raw, into = c("CHROM", "start", "attrs"), sep = "\t", extra = "merge") %>%
  mutate(start = as.integer(start),
         gene_id = str_remove(str_match(attrs, 'gene_id "([^"]+)"')[, 2], "\\.\\d+$"),
         arm = assign_arm(CHROM, start)) %>%
  filter(!is.na(arm)) %>%
  dplyr::select(gene_id, arm)

vsd <- readRDS(file.path(RDS_DIR, "R02_vst_counts.rds"))
expr_mat <- assay(vsd)
sample_genotype <- setNames(as.character(colData(vsd)$genotype), colnames(vsd))

centered <- expr_mat - rowMeans(expr_mat)  # gene-wise: R recycles this nrow(expr_mat)-length vector down each column
dosage_long <- as.data.frame(centered) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "centered_vst") %>%
  inner_join(gene_pos, by = "gene_id") %>%
  mutate(genotype = factor(sample_genotype[sample_id], levels = GENOTYPE_LEVELS),
         arm = factor(arm, levels = arm_order))

dosage_by_arm <- dosage_long %>%
  group_by(sample_id, genotype, arm) %>%
  summarise(n_genes = n(), mean_centered_vst = mean(centered_vst), .groups = "drop")
write_csv(dosage_by_arm, file.path(TABLES_DIR, "R08_expression_dosage_by_arm.csv"))
dosage_informative <- dosage_by_arm %>% filter(n_genes >= ESNP_MIN_GENES_PER_ARM)

########
## 5. Plots
########

cat("\n[R08] Plotting...\n")

sample_order <- meta %>% mutate(genotype = factor(genotype, levels = GENOTYPE_LEVELS)) %>%
  arrange(genotype, sample_id) %>% pull(sample_id)

# Genome-wide BAF track, faceted by chromosome, centromere marked -- the classic
# eSNP-Karyotyping ideogram-style figure. size/alpha tuned down: tens of thousands
# of SNVs per sample x 9 samples is a lot of points for a single PDF.
p_genomewide <- ggplot(snvs, aes(x = POS / 1e6, y = baf, color = genotype)) +
  geom_point(size = 0.4, alpha = 0.15) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(data = centromere, aes(xintercept = centromere / 1e6), color = "grey30", linewidth = 0.3) +
  facet_wrap(~CHROM, scales = "free_x", ncol = 4) +
  scale_color_manual(values = GENOTYPE_COLOR, name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(title = "eSNP-Karyotyping: genome-wide B-allele frequency by chromosome",
       subtitle = sprintf("Heterozygous PASS biallelic SNVs, depth >= %d. Vertical line = centromere.", ESNP_MIN_DEPTH),
       x = "Position (Mb)", y = "BAF (alt / (ref+alt) depth)") +
  theme_bw(base_size = 10) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(PLOTS_DIR, "R08_eSNP_BAF_genomewide.pdf"), p_genomewide, width = 16, height = 12)

# Per-sample summary -- one weighted-mean |BAF-0.5| per sample, genome-wide (all
# arms pooled, weighted by each arm's own SNV count). The genome-wide track and
# the 46-arm heatmap are both dense/hard to eyeball for "is any one sample off";
# this collapses the same underlying numbers to 9 points, directly answering that
# question. Dashed line = the flag threshold, for scale -- note this genome-wide
# average is expected to sit near the SAME baseline noise level for every sample
# regardless of genotype (a locus-specific aneuploidy signal would show up as one
# or two genotype's ARMS elevated in R08_BAF_heatmap.pdf, not necessarily in this
# genome-wide average) -- this plot is a quick outlier/QC check, not the primary
# aneuploidy test itself.
sample_summary <- baf_by_arm %>%
  group_by(sample_id, genotype) %>%
  summarise(total_snps = sum(n_snps), mean_dev = weighted.mean(mean_abs_dev, n_snps), .groups = "drop") %>%
  mutate(sample_id = factor(sample_id, levels = rev(sample_order)))
p_summary <- ggplot(sample_summary, aes(x = mean_dev, y = sample_id, color = genotype)) +
  geom_vline(xintercept = ESNP_BAF_FLAG, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  scale_color_manual(values = GENOTYPE_COLOR, name = NULL) +
  labs(title = "eSNP-Karyotyping: per-sample BAF deviation (outlier check)",
       subtitle = sprintf("Weighted mean |BAF-0.5|, all arms pooled (weight = arm's SNV count).\nDashed line = flag threshold (%.2f).", ESNP_BAF_FLAG),
       x = "Weighted mean |BAF - 0.5|", y = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(PLOTS_DIR, "R08_BAF_sample_summary.pdf"), p_summary, width = 8, height = 4.5)
write_csv(sample_summary, file.path(TABLES_DIR, "R08_BAF_sample_summary.csv"))

# Sample x arm heatmaps -- BAF deviation and expression dosage, same grid, side by
# side comparable. Genotype/arm order fixed (not clustered) so genomic order and
# genotype grouping stay legible, same convention as R04's concordant-DEG heatmap.
make_arm_matrix <- function(df, value_col) {
  wide <- df %>% dplyr::select(sample_id, arm, value = all_of(value_col)) %>%
    pivot_wider(names_from = arm, values_from = value)
  mat <- as.matrix(wide %>% dplyr::select(-sample_id))
  rownames(mat) <- wide$sample_id
  mat[sample_order, intersect(arm_order, colnames(mat)), drop = FALSE]
}

row_ann <- data.frame(Genotype = factor(sample_genotype[sample_order], levels = GENOTYPE_LEVELS), row.names = sample_order)
ann_colors <- list(Genotype = GENOTYPE_COLOR)

# Unfiltered (not `informative`) -- shows every cell's raw value regardless of
# n_snps, including below-threshold ones (that threshold still gates the actual
# flagging logic above, unaffected here; this view is for visual inspection of
# what a below-threshold cell's real number looks like, e.g. to compare it
# against its own genotype's other replicates rather than just seeing blank).
# NA_COLOR deliberately dark/unambiguous, not a pale grey -- both palettes below
# run light-to-dark from white (BAF) or through white at the diverging midpoint
# (dosage), so a pale grey NA cell sits almost exactly where a real "near-zero"
# data value would land and reads as "no deviation" instead of "no data". Dark
# grey can't be confused with either palette's own colors at any point in their
# range. Distinct from GENOTYPE_COLOR's WT swatch (grey40) too, and it's in a
# separate annotation strip so the two are never spatially adjacent anyway.
NA_COLOR <- "grey30"

baf_mat <- make_arm_matrix(baf_by_arm, "mean_abs_dev")
pdf(file.path(PLOTS_DIR, "R08_BAF_heatmap.pdf"), width = 12, height = 5)
pheatmap(baf_mat, cluster_rows = FALSE, cluster_cols = FALSE, na_col = NA_COLOR,
         annotation_row = row_ann, annotation_colors = ann_colors,
         color = colorRampPalette(c("white", "#E41A1C"))(50),
         main = sprintf("eSNP-Karyotyping: mean |BAF-0.5| by arm (dark grey = no SNVs observed; flagging requires >=%d informative)", ESNP_MIN_SNPS_PER_ARM),
         fontsize = 8, angle_col = 90)
dev.off()

dosage_mat <- make_arm_matrix(dosage_informative, "mean_centered_vst")
pdf(file.path(PLOTS_DIR, "R08_expression_dosage_heatmap.pdf"), width = 12, height = 5)
pheatmap(dosage_mat, cluster_rows = FALSE, cluster_cols = FALSE, na_col = NA_COLOR,
         annotation_row = row_ann, annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = sprintf("Expression dosage: mean gene-centered VST by arm (dark grey = <%d genes; Ben-David et al. 2014 signal)", ESNP_MIN_GENES_PER_ARM),
         fontsize = 8, angle_col = 90)
dev.off()

########
## 6. Text summary
########

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "R08_eSNP_Karyotype_session_info.txt"))

sink(SUMMARY)
cat("================================================================================\n")
cat("eSNP-KARYOTYPING -- RNA-seq ANEUPLOIDY/KARYOTYPE QC SCREEN\n")
cat("Project: KS1_2_Neuro-ectoderm\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("================================================================================\n\n")

cat("Method: eSNP-Karyotyping (Weissbein, Benvenisty & Ben-David, Stem Cell Reports\n")
cat("2016) -- heterozygous-SNV B-allele frequency binned by chromosome arm, called\n")
cat("from RNA-seq via GATK4 RNA-seq short-variant best practices (nf-core/rnavar).\n")
cat("Companion signal: expression-dosage by arm (Ben-David et al., Cell Stem Cell\n")
cat("2014) -- aneuploid regions show elevated mean expression.\n\n")

cat("Sample breakdown by genotype:\n")
print(meta %>% dplyr::count(genotype))
cat("\n")

cat("Thresholds:\n")
cat(sprintf("  Min depth (ref+alt) at a het SNV: %d\n", ESNP_MIN_DEPTH))
cat(sprintf("  Min informative SNVs per arm: %d\n", ESNP_MIN_SNPS_PER_ARM))
cat(sprintf("  Flag: mean |BAF-0.5| > %.2f in all 3 replicates of one genotype, not also in WT\n\n", ESNP_BAF_FLAG))

cat(sprintf("Candidate arms flagged: %d\n", nrow(flagged_arms)))
if (nrow(flagged_arms) > 0) {
  print(as.data.frame(flagged_arms))
} else {
  cat("  None -- no chromosome arm showed a consistent BAF shift in any mutant\n")
  cat("  genotype beyond WT at the thresholds above.\n")
}
cat("\n")

cat("Caveats:\n")
cat("  - Screening-level QC, not a diagnostic test: n=3 replicates/genotype, soft\n")
cat("    threshold flags (matching eSNP-Karyotyping's own threshold-based design),\n")
cat("    not a hypothesis test with a p-value.\n")
cat("  - chrX ploidy differs by line sex; chrY and chrM excluded (uninformative for\n")
cat("    a heterozygous-SNV BAF scan). Interpret chrX flags with that in mind.\n")
cat("  - Before trusting any flagged arm in a mutant line, confirm WT replicates\n")
cat("    cluster near BAF=0.5 genome-wide in R08_eSNP_BAF_genomewide.pdf (positive\n")
cat("    control that the method itself is unbiased in this dataset).\n")
cat("================================================================================\n")
cat("END OF SUMMARY\n")
sink()

cat("\n[DONE] R08_eSNP_Karyotype complete -- summary: ", SUMMARY, "\n", sep = "")
