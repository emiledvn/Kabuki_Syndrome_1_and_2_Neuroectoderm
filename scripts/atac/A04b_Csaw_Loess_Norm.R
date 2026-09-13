#!/usr/bin/env Rscript
# A04b_Csaw_Loess_Norm.R -- non-linear (loess) alternative to A04's Csaw TMM
# normalization, added after the background-bin-only MA plot for KDM6A_ko_vs_WT
# showed a clear upward tail at high average abundance even outside called peaks --
# i.e. an intensity-dependent (efficiency) bias that a single TMM scale factor per
# sample cannot correct, since TMM only removes a *constant* fold-change across all
# abundances. See docs/A04b_normalization_methodology.md for the full writeup
# (motivation, method, and before/after comparison -- written to double as the
# thesis methods-section draft).
#
# Method: csaw::normOffsets() fits a per-sample loess trend of log-count vs.
# log-average-count on the SAME unfiltered background bins A04's Csaw arm already
# uses for its TMM factors, then interpolates that trend onto the consensus peak
# set via spline interpolation (se.out=<peak SummarizedExperiment>) -- this is
# csaw's own documented mechanism for fitting a trend on one region set (here:
# genome-wide background bins) and applying it to a different one (here: peaks),
# the exact same idiom the csaw docs describe for spike-in-referenced normalization,
# just with background bins standing in for spike-in material. The resulting
# per-peak-per-sample offsets are converted to DESeq2 normalizationFactors (see
# below) and a fresh dispersion/Wald fit is run -- normalizationFactors can't be
# swapped into an already-fitted DESeqDataSet, since dispersion estimation itself
# depends on them.
#
# Peak counts are taken from A04's Csaw-normalized DESeq2$DEdata directly (NOT from
# dba.peakset(bRetrieve=TRUE)), following the same precedent as A05_DESeq2.R -- an
# earlier version of A05 rebuilt counts from dba.peakset() and silently diverged
# from DiffBind's own fit (~1.33x off on most cells, root cause not fully traced).
# Loading DEdata's own counts sidesteps that mismatch entirely.
#
# This is an ADDITIONAL comparison arm, not a replacement -- A04/A05 are untouched.
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/A04b_Csaw_Loess_Norm.R
#
# Self-checkpointing: skips entirely if results/tables/A04b_summary.tsv exists.
# Background-bin counting (shared with A04's TMM arm) is itself cached to
# data/diffbind_output/A04b_background_bins.rds.

suppressPackageStartupMessages({
  library(DiffBind)
  library(csaw)
  library(SummarizedExperiment)
  library(DESeq2)
  library(BiocParallel)
  library(yaml)
  library(dplyr); library(tibble); library(readr)
  library(ggplot2)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[A04b] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

CSAW_DIR    <- "data/diffbind_output/Csaw_Norm"
CSAW_RDS    <- file.path(CSAW_DIR, "diffbind_analyzed.rds")
OUT_BASE    <- "data/diffbind_output"
LOESS_DIR   <- file.path(OUT_BASE, "CsawLoess_Norm")
BIN_CACHE   <- file.path(OUT_BASE, "A04b_background_bins.rds")
PLOTS_DIR   <- "results/atac"
TABLES_DIR  <- "results/tables"
SESSION_DIR <- "results/Session_info"
SUMMARY     <- file.path(TABLES_DIR, "A04b_summary.tsv")

for (d in c(LOESS_DIR, PLOTS_DIR, TABLES_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(SUMMARY)) {
  cat("[A04b] Already complete (", SUMMARY, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}
if (!file.exists(CSAW_RDS)) stop("[A04b] ERROR: ", CSAW_RDS, " not found -- run A04_Diffbind_Compare_NORMS.R first.")

CORES       <- cfg$diffbind$cores
FDR         <- cfg$thresholds$fdr
CONTRASTS   <- cfg$contrasts
register(MulticoreParam(workers = CORES))

cat("================================================================================\n")
cat("A04b_Csaw_Loess_Norm -- non-linear (loess) alternative to Csaw TMM\n")
cat("================================================================================\n\n")

cat("[A04b] Loading Csaw-normalized DiffBind object:", CSAW_RDS, "\n")
atac_csaw  <- readRDS(CSAW_RDS)
dds_tmm    <- atac_csaw$DESeq2$DEdata
peaks_gr   <- dba.peakset(atac_csaw, bRetrieve = TRUE)
if (length(peaks_gr) != nrow(dds_tmm)) stop("[A04b] ERROR: peak count mismatch between dba.peakset() and DESeq2$DEdata.")
rownames(dds_tmm) <- sprintf("%s:%d-%d", as.character(seqnames(peaks_gr)), start(peaks_gr), end(peaks_gr))

raw_peak_counts <- counts(dds_tmm)  # DiffBind's own fitted counts -- authoritative, see header note
sample_ids      <- colnames(dds_tmm)
bam_files       <- atac_csaw$samples$bamReads[match(sample_ids, atac_csaw$samples$SampleID)]

########
## Background bins -- identical construction to A04's Csaw arm (same readParam,
## same bin width, same "drop only structurally-empty bins" filter), cached since
## windowCounts() over 9 BAMs genome-wide is the slow part of this script.
########

if (file.exists(BIN_CACHE)) {
  cat("[A04b] Loading cached background bins:", BIN_CACHE, "\n")
  bin_filtered <- readRDS(BIN_CACHE)
} else {
  cat("[A04b] Counting reads in 10kb background bins (cached afterwards)...\n")
  read_params <- readParam(minq = 30, dedup = FALSE, pe = "both", max.frag = 2000)
  bin_counts  <- windowCounts(bam_files, bin = TRUE, width = 10000, param = read_params,
                               BPPARAM = MulticoreParam(workers = CORES))
  colnames(bin_counts) <- sample_ids
  bin_keep     <- rowSums(assay(bin_counts)) > 0
  bin_filtered <- bin_counts[bin_keep, ]
  saveRDS(bin_filtered, BIN_CACHE)
}

########
## Fit the loess trend on background bins, interpolate onto the peak set.
## normOffsets(object, se.out=<other SE>) computes the trend from `object` and
## applies it via spline interpolation to `se.out` -- csaw's documented mechanism
## for exactly this "fit on one region set, apply to another" case (its own docs
## describe this for spike-in-referenced normalization; background bins play the
## same anchoring role here). Requires se.out$totals == object$totals exactly, so
## we deliberately anchor the peak SE's library sizes to the same background-bin
## totals already used for the fit, keeping both on the same normalization basis.
##
## normOffsets positions each region on its fitted curve using RAW counts, not
## width-normalized density. Background bins are all a fixed 10kb; MACS3 narrowPeak
## calls are far narrower and highly variable in width (~150bp-2kb), so a peak's
## raw count sits at a systematically lower position than a background bin of equal
## true accessibility purely because it spans less genomic sequence -- this
## initially produced a nonsensical correction (verified: peak abundance median
## -0.29 vs. background-bin median 1.26 on the same log2-CPM scale, despite peaks
## being the enriched regions). Fixed by scaling each peak's raw count to a
## per-10kb-equivalent pseudo-count (same density, as if the peak were 10kb wide)
## SOLELY for positioning on the curve -- the resulting offsets are dimensionless
## log-corrections and get applied to the true raw peak counts below, same as
## before; only the curve-lookup step needed the width correction.
########

peak_widths_bp    <- width(peaks_gr)
peak_counts_scaled <- round(raw_peak_counts * (10000 / peak_widths_bp))

peak_se <- SummarizedExperiment(assays = list(counts = peak_counts_scaled))
peak_se$totals <- bin_filtered$totals

cat("[A04b] Fitting loess trend on background bins, interpolating onto", nrow(peak_se), "peaks (width-normalized to 10kb-equivalent density)...\n")
peak_se <- normOffsets(bin_filtered, se.out = peak_se)
offset_matrix <- assay(peak_se, "offset")  # natural-log scale, comparable to log-library-size

## DESeq2 normalizationFactors are linear-scale, per-gene-per-sample multipliers,
## conventionally row-geometric-mean-centered to 1 so they carry only the relative
## sample-to-sample correction (not an overall magnitude shift) -- same convention
## DESeq2 itself uses internally for sizeFactors.
nf <- exp(offset_matrix)
nf <- nf / exp(rowMeans(log(nf)))
stopifnot(all(is.finite(nf)), all(nf > 0))

########
## Fresh DESeqDataSet -- normalizationFactors can't be assigned post-hoc onto an
## already dispersion-fitted object, so this is a new fit, not a patch of dds_tmm.
## Same counts, same colData/design as the TMM fit (colData(dds_tmm) already carries
## DiffBind's own Condition factor + levels), same fitType DiffBind itself used
## (see A05_DESeq2.R header note: DiffBind's pv.DESeq2design runs
## estimateDispersions(fitType="local") + nbinomWaldTest) so the two arms are
## comparable on identical fitting machinery, differing only in normalization.
########

dds_loess <- DESeqDataSetFromMatrix(countData = raw_peak_counts,
                                     colData   = colData(dds_tmm),
                                     design    = ~Condition)
normalizationFactors(dds_loess) <- nf
cat("[A04b] Running DESeq2 (loess-normalized)...\n")
dds_loess <- estimateDispersions(dds_loess, fitType = "local")
dds_loess <- nbinomWaldTest(dds_loess)

dir.create(LOESS_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(dds_loess, file.path(LOESS_DIR, "dds_loess.rds"))
write.csv(data.frame(Sample = sample_ids, RowGeomMeanOffset = colMeans(log(nf))),
          file.path(LOESS_DIR, "loess_offset_summary.csv"), row.names = FALSE, quote = FALSE)

########
## Per-contrast results + before/after MA plots, mirroring A05's extract_results()/
## make_ma_plot() pattern so the two arms are visually and numerically comparable.
########

extract <- function(dds, treatment, reference, norm_label) {
  res <- results(dds, contrast = c("Condition", treatment, reference))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("peak_id") %>%
    mutate(norm = norm_label) %>%
    arrange(padj)
  res_df
}

make_ma <- function(res_df, name, norm_label, suffix) {
  res_plot <- res_df %>% mutate(sig = !is.na(padj) & padj < FDR)
  p <- ggplot(res_plot, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = sig), alpha = 0.5, size = 1) +
    scale_x_log10() +
    scale_color_manual(values = c("TRUE" = "#e31a1c", "FALSE" = "grey70"), name = sprintf("padj < %.2f", FDR)) +
    geom_hline(yintercept = 0, color = "grey30") +
    geom_smooth(method = "loess", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
    theme_bw(base_size = 11) +
    labs(title = paste0(name, " -- ", norm_label), subtitle = "MA plot",
         x = "Mean normalized counts", y = "log2 Fold Change")
  ggsave(file.path(PLOTS_DIR, paste0("A04b_MA_", suffix, ".pdf")), p, width = 7, height = 5)
}

summary_rows <- list()
for (ct in CONTRASTS) {
  res_tmm   <- extract(dds_tmm,   ct$treatment, ct$reference, "TMM")
  res_loess <- extract(dds_loess, ct$treatment, ct$reference, "Loess")

  write_csv(res_loess, file.path(TABLES_DIR, sprintf("A04b_%s_CsawLoess_full.csv", ct$name)))

  make_ma(res_tmm,   ct$name, "TMM (A04)",     sprintf("%s_TMM", ct$name))
  make_ma(res_loess, ct$name, "Loess (A04b)",  sprintf("%s_Loess", ct$name))

  for (norm_label in c("TMM", "Loess")) {
    d <- if (norm_label == "TMM") res_tmm else res_loess
    sig <- d %>% filter(!is.na(padj), padj < FDR)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Contrast = ct$name, Norm = norm_label,
      Gained = sum(sig$log2FoldChange > 0),
      Lost   = sum(sig$log2FoldChange < 0)
    )
  }
  cat(sprintf("[A04b] %s -- TMM: %d gained / %d lost | Loess: %d gained / %d lost\n",
              ct$name,
              sum(res_tmm$padj   < FDR & res_tmm$log2FoldChange   > 0, na.rm = TRUE),
              sum(res_tmm$padj   < FDR & res_tmm$log2FoldChange   < 0, na.rm = TRUE),
              sum(res_loess$padj < FDR & res_loess$log2FoldChange > 0, na.rm = TRUE),
              sum(res_loess$padj < FDR & res_loess$log2FoldChange < 0, na.rm = TRUE)))
}

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, SUMMARY, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n[A04b] Summary (TMM vs Loess):\n")
print(summary_df)

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "A04b_Csaw_Loess_Norm_session_info.txt"))
cat("\n[DONE] A04b_Csaw_Loess_Norm complete\n")
