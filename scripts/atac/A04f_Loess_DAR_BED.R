#!/usr/bin/env Rscript
# A04f_Loess_DAR_BED.R -- GAINED/LOST BED files for the Loess-normalized DARs
# (A04b), matching the convention A04_Diffbind_Compare_NORMS.R already uses for
# its TMM/Default/Background arms (FDR-significant, split by sign of log2FC,
# written via rtracklayer::export.bed). A04b itself only writes the full
# (unfiltered) results table -- these per-direction BED files didn't exist yet
# for Loess, needed for direction-specific downstream work (e.g. motif
# enrichment split by GAINED vs LOST).
#
# Read-only against A04b's output; does not refit anything, does not touch
# TOBIAS or its inputs.
#
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/A04f_Loess_DAR_BED.R
# Requires: A04b_Csaw_Loess_Norm.R already run.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
  library(GenomicRanges); library(rtracklayer); library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[A04f] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)
FDR <- cfg$thresholds$fdr
LFC_THR <- cfg$thresholds$lfc_threshold
CONTRASTS <- cfg$contrasts

TABLES_DIR <- "results/tables"
OUT_DIR    <- "data/diffbind_output/CsawLoess_Norm"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("A04f_Loess_DAR_BED -- GAINED/LOST BED files for Loess-normalized DARs\n")
cat("================================================================================\n\n")

peak_id_to_gr <- function(peak_id) {
  m <- regmatches(peak_id, regexec("^(.+):(\\d+)-(\\d+)$", peak_id))
  GRanges(seqnames = sapply(m, `[`, 2),
          ranges = IRanges(start = as.numeric(sapply(m, `[`, 3)), end = as.numeric(sapply(m, `[`, 4))))
}

for (ct in CONTRASTS) {
  full_path <- file.path(TABLES_DIR, sprintf("A04b_%s_CsawLoess_full.csv", ct$name))
  if (!file.exists(full_path)) { cat("[A04f] Skipping ", ct$name, " -- ", full_path, " not found.\n"); next }

  res <- read_csv(full_path, show_col_types = FALSE) %>%
    filter(!is.na(padj), padj < FDR, abs(log2FoldChange) > LFC_THR)
  gained <- res %>% filter(log2FoldChange > 0)
  lost   <- res %>% filter(log2FoldChange < 0)

  export.bed(peak_id_to_gr(gained$peak_id), file.path(OUT_DIR, sprintf("%s_GAINED.bed", ct$name)))
  export.bed(peak_id_to_gr(lost$peak_id),   file.path(OUT_DIR, sprintf("%s_LOST.bed", ct$name)))

  cat(sprintf("[A04f] %s (Loess, padj < %.2f, |log2FC| > %.2f): %d GAINED, %d LOST -> %s\n",
              ct$name, FDR, LFC_THR, nrow(gained), nrow(lost), OUT_DIR))
}

cat("\n[DONE] A04f_Loess_DAR_BED complete\n")
