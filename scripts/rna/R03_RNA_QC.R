#!/usr/bin/env Rscript
# R03_RNA_QC.R -- RNA QC summary, consolidating nf-core/rnaseq's own already-computed
# QC (STAR alignment stats, Picard duplication) rather than recomputing by hand --
# same philosophy as A03_QC_ATAC.sh for the ATAC arm. No source script existed to
# migrate this from (unlike most of this pipeline), so it's written fresh against
# nf-core/rnaseq 3.26.0's standard MultiQC output structure.
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/R03_RNA_QC.R
#
# Column names in nf-core's multiqc_general_stats.txt are matched flexibly
# (case-insensitive substring on tool + metric name) rather than hardcoded exactly,
# since exact column naming has changed across MultiQC/nf-core versions in the past
# and this project has already been burned once (see A03/README "Methods notes") by
# trusting an unverified upstream metric without checking its exact provenance --
# if the expected columns aren't found, this fails with a clear message showing
# what *was* found, rather than silently producing an empty/wrong table.
#
# Self-checkpointing: skips entirely if results/tables/R03_qc_summary.tsv exists.
# Requires R01_RUN_nfcore_RNA.sh to have completed first.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[R03] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

MULTIQC_DIR <- "data/RNA_nfcore_output/multiqc/star_salmon/multiqc_report_data"
GENERAL_STATS <- file.path(MULTIQC_DIR, "multiqc_general_stats.txt")
TABLES_DIR <- "results/tables"
SESSION_DIR <- "results/Session_info"
OUT <- file.path(TABLES_DIR, "R03_qc_summary.tsv")

for (d in c(TABLES_DIR, SESSION_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (file.exists(OUT)) {
  cat("[R03] Already complete (", OUT, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}
if (!dir.exists(MULTIQC_DIR)) stop("[R03] ERROR: ", MULTIQC_DIR, " not found -- run R01_RUN_nfcore_RNA.sh first.")
if (!file.exists(GENERAL_STATS)) stop("[R03] ERROR: ", GENERAL_STATS, " not found. nf-core/rnaseq's MultiQC output structure may differ from what this script expects -- check ", MULTIQC_DIR, " manually and update the path here.")

MIN_UNIQUE_MAPPED <- cfg$thresholds$rna_qc_min_unique_mapped_pct

cat("================================================================================\n")
cat("R03_RNA_QC\n")
cat("================================================================================\n\n")

raw <- read_tsv(GENERAL_STATS, show_col_types = FALSE)
if (!"Sample" %in% names(raw)) stop("[R03] ERROR: no 'Sample' column in ", GENERAL_STATS, ". Columns found: ", paste(names(raw), collapse = ", "))

find_col <- function(df, ...) {
  patterns <- list(...)
  for (p in patterns) {
    hit <- names(df)[str_detect(tolower(names(df)), p)]
    if (length(hit) >= 1) return(hit[1])
  }
  NA_character_
}

col_unique_mapped <- find_col(raw, "star.*uniquely_mapped_percent", "uniquely_mapped_percent")
col_dup           <- find_col(raw, "picard.*percent_duplication", "percent_duplication")
col_gc            <- find_col(raw, "fastqc.*percent_gc", "percent_gc")

if (is.na(col_unique_mapped)) {
  stop("[R03] ERROR: could not find a STAR uniquely-mapped-percent column in ", GENERAL_STATS,
       ". Columns found: ", paste(names(raw), collapse = ", "),
       "\n[R03] Update find_col() patterns above once the real column name is known.")
}

qc <- raw %>%
  transmute(
    Sample = Sample,
    Unique_mapped_pct = .data[[col_unique_mapped]],
    Duplication_pct   = if (!is.na(col_dup)) .data[[col_dup]] else NA_real_,
    GC_pct            = if (!is.na(col_gc)) .data[[col_gc]] else NA_real_
  ) %>%
  # multiqc_general_stats.txt also carries per-read FastQC-only rows (e.g. "<sample>
  # Read 1/2") that have no STAR alignment stats -- drop those so each real sample
  # is counted once.
  filter(!is.na(Unique_mapped_pct)) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  mutate(Unique_mapped_OK = ifelse(Unique_mapped_pct > MIN_UNIQUE_MAPPED, "PASS", "WARN"))

if (nrow(qc) == 0) stop("[R03] ERROR: parsed 0 samples from ", GENERAL_STATS, " -- check the file manually.")

write_tsv(qc, OUT)
cat(sprintf("[R03] %d samples summarized\n", nrow(qc)))
cat(sprintf("[R03] Mean unique-mapped: %.1f%% (threshold: >%.0f%%, reference value not independently verified -- see config comment)\n",
            mean(qc$Unique_mapped_pct, na.rm = TRUE), MIN_UNIQUE_MAPPED))
if (any(qc$Unique_mapped_OK == "WARN")) {
  cat("[R03] WARN samples:\n")
  print(qc %>% filter(Unique_mapped_OK == "WARN") %>% select(Sample, Unique_mapped_pct))
}

writeLines(capture.output(sessionInfo()), file.path(SESSION_DIR, "R03_RNA_QC_session_info.txt"))
cat("\n[DONE] R03_RNA_QC complete -- ", OUT, "\n", sep = "")
