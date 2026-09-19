#!/usr/bin/env Rscript
# A04_Diffbind_Compare_NORMS.R -- DiffBind differential accessibility, comparing
# three normalization strategies (Default, Background, Csaw TMM) so the choice of
# normalization is a documented decision, not an assumption.
# Adapted from ED26_001_ATAC_c/Narrow_allReads/scripts/03_diffbind_allReads.R and
# 04_diffbind_csaw_allReads.R (merged into one script here, one per code convention).
# Env: ks_1_2_r. Run from the repo root: Rscript scripts/atac/A04_Diffbind_Compare_NORMS.R
#
# Self-checkpointing: skips entirely if results/tables/A04_summary.tsv exists.
# The expensive peak-counting step is itself cached to
# data/diffbind_output/diffbind_counted_base.rds, so an interruption after that
# point doesn't force a recount.

suppressPackageStartupMessages({
  library(DiffBind)
  library(csaw)
  library(edgeR)
  library(rtracklayer)
  library(BiocParallel)
  library(yaml)
})

CONFIG <- "config/pipeline_config.yaml"
if (!file.exists(CONFIG)) stop("[A04] ERROR: ", CONFIG, " not found -- run this script from the repo root.")
cfg <- yaml::read_yaml(CONFIG)

METADATA   <- "data/sample_metadata.csv"
BAM_DIR    <- "data/ATAC_nfcore_output/bwa/merged_library"
PEAK_DIR   <- "data/atac_peaks"
OUT_BASE   <- "data/diffbind_output"
TABLES_DIR <- "results/tables"
RDS_BASE   <- file.path(OUT_BASE, "diffbind_counted_base.rds")
SUMMARY    <- file.path(TABLES_DIR, "A04_summary.tsv")

if (file.exists(SUMMARY)) {
  cat("[A04] Already complete (", SUMMARY, " exists). Skipping. Delete that file to force a rerun.\n", sep = "")
  quit(save = "no", status = 0)
}
if (!file.exists(METADATA)) stop("[A04] ERROR: ", METADATA, " not found.")
if (!dir.exists(BAM_DIR))   stop("[A04] ERROR: ", BAM_DIR, " not found -- run A01_RUN_nfcore_ATAC.sh first.")
if (!dir.exists(PEAK_DIR))  stop("[A04] ERROR: ", PEAK_DIR, " not found -- run A02_MACS3_narrow.sh first.")

dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

CORES     <- cfg$diffbind$cores
FDR       <- cfg$thresholds$fdr
CONTRASTS <- cfg$contrasts   # list of {name, treatment, reference}, shared with RNA (R02/R04)

register(MulticoreParam(workers = CORES))

########
## Build the DiffBind samplesheet directly from sample_metadata.csv (single
## source of truth) instead of maintaining a separate diffbind_samplesheet.csv.
########

meta      <- read.csv(METADATA, stringsAsFactors = FALSE)
atac_meta <- meta[meta$assay == "ATAC", ]

samplesheet <- data.frame(
  SampleID   = atac_meta$sample_id,
  Condition  = atac_meta$genotype,
  Replicate  = atac_meta$replicate,
  bamReads   = file.path(BAM_DIR, paste0(atac_meta$sample_id, "_REP1.mLb.clN.sorted.bam")),
  Peaks      = file.path(PEAK_DIR, paste0(atac_meta$sample_id, "_peaks.narrowPeak")),
  PeakCaller = "bed",
  stringsAsFactors = FALSE
)

missing_bams  <- samplesheet$bamReads[!file.exists(samplesheet$bamReads)]
missing_peaks <- samplesheet$Peaks[!file.exists(samplesheet$Peaks)]
if (length(missing_bams)  > 0) stop("[A04] ERROR: missing BAM(s):\n",  paste(missing_bams,  collapse = "\n"))
if (length(missing_peaks) > 0) stop("[A04] ERROR: missing peak file(s):\n", paste(missing_peaks, collapse = "\n"))

########
## Load or count
########

if (file.exists(RDS_BASE)) {
  cat("[A04] Loading pre-counted object:", RDS_BASE, "\n")
  atac <- readRDS(RDS_BASE)
} else {
  cat("[A04] Counting reads in peaks (slow step; cached afterwards)...\n")
  atac <- dba(sampleSheet = samplesheet)
  atac$config$cores <- CORES
  atac <- dba.count(atac, bUseSummarizeOverlaps = TRUE)
  saveRDS(atac, RDS_BASE)
  cat("[A04] Saved:", RDS_BASE, "\n")
}
atac$config$cores <- CORES

########
## Normalization + DESeq2, one normalization method at a time
########

run_and_save <- function(dba_obj, norm_type, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("==============================================\n")
  cat(" DiffBind --", norm_type, "\n")
  cat("==============================================\n")

  dba_obj <- if (norm_type == "Background") dba.normalize(dba_obj, background = TRUE) else dba.normalize(dba_obj)

  for (ct in CONTRASTS) dba_obj <- dba.contrast(dba_obj, contrast = c("Condition", ct$treatment, ct$reference))

  cat("[A04] Running DESeq2...\n")
  dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2)

  for (i in seq_along(CONTRASTS)) {
    ct  <- CONTRASTS[[i]]
    res <- dba.report(dba_obj, contrast = i, th = FDR, bUsePval = FALSE, bCounts = TRUE)
    write.csv(as.data.frame(res), file.path(TABLES_DIR, sprintf("A04_DiffBind_%s_%s.csv", ct$name, norm_type)),
              quote = FALSE, row.names = FALSE)
    export.bed(res[res$Fold > 0], file.path(out_dir, sprintf("%s_GAINED.bed", ct$name)))
    export.bed(res[res$Fold < 0], file.path(out_dir, sprintf("%s_LOST.bed", ct$name)))
  }

  saveRDS(dba_obj, file.path(out_dir, "diffbind_analyzed.rds"))
  cat("[A04] ", norm_type, " complete\n\n")
  invisible(dba_obj)
}

run_and_save(atac, "Default",    file.path(OUT_BASE, "Default_Norm"))
run_and_save(atac, "Background", file.path(OUT_BASE, "Background_Norm"))

########
## Csaw TMM normalization (background 10kb bins, independent of peak calls)
########

cat("[A04] Csaw: counting reads in 10kb background bins...\n")

bam_files  <- atac$samples$bamReads
sample_ids <- atac$samples$SampleID

read_params <- readParam(minq = 30, dedup = FALSE, pe = "both", max.frag = 2000)
bin_counts  <- windowCounts(bam_files, bin = TRUE, width = 10000, param = read_params,
                             BPPARAM = MulticoreParam(workers = CORES))
colnames(bin_counts) <- sample_ids

## TMM factors must be computed on UNFILTERED bins (csaw's composition-bias
## normalization relies on "most bins are background"; the csaw book explicitly
## warns that computing factors from high-abundance/enriched windows and applying
## them to background regions gives incorrect normalization -- filterWindowsGlobal's
## enrichment filter belongs on the later differential-window TEST, not here). Only
## drop bins with zero coverage in every sample -- structurally uninformative for
## TMM, not an enrichment-based filter.
bin_keep     <- rowSums(assay(bin_counts)) > 0
bin_filtered <- bin_counts[bin_keep, ]
y            <- asDGEList(bin_filtered)
y            <- calcNormFactors(y, method = "TMM")
norm_factors <- y$samples$norm.factors
lib_sizes    <- y$samples$lib.size

csaw_out <- file.path(OUT_BASE, "Csaw_Norm")
dir.create(csaw_out, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(Sample = sample_ids, LibSize = lib_sizes, TMM_Factor = norm_factors),
          file.path(csaw_out, "csaw_tmm_norm_factors.csv"), row.names = FALSE, quote = FALSE)

effective_lib_sizes <- lib_sizes * norm_factors
atac_csaw <- dba.normalize(atac, normalize = DBA_NORM_LIB, library = effective_lib_sizes, background = FALSE)

for (ct in CONTRASTS) atac_csaw <- dba.contrast(atac_csaw, contrast = c("Condition", ct$treatment, ct$reference))

cat("[A04] Running DESeq2 (Csaw-normalized)...\n")
atac_csaw <- dba.analyze(atac_csaw, method = DBA_DESEQ2)

for (i in seq_along(CONTRASTS)) {
  ct  <- CONTRASTS[[i]]
  res <- dba.report(atac_csaw, contrast = i, th = FDR, bUsePval = FALSE, bCounts = TRUE)
  write.csv(as.data.frame(res), file.path(TABLES_DIR, sprintf("A04_DiffBind_%s_Csaw.csv", ct$name)),
            quote = FALSE, row.names = FALSE)
  export.bed(res[res$Fold > 0], file.path(csaw_out, sprintf("%s_GAINED.bed", ct$name)))
  export.bed(res[res$Fold < 0], file.path(csaw_out, sprintf("%s_LOST.bed", ct$name)))
}
saveRDS(atac_csaw, file.path(csaw_out, "diffbind_analyzed.rds"))

########
## Summary across all three normalizations
########

summary_rows <- list()
for (ct in CONTRASTS) {
  for (norm in c("Default", "Background", "Csaw")) {
    f <- file.path(TABLES_DIR, sprintf("A04_DiffBind_%s_%s.csv", ct$name, norm))
    if (file.exists(f)) {
      d <- read.csv(f)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Contrast = ct$name, Norm = norm,
        Gained = sum(d$FDR < FDR & d$Fold > 0),
        Lost   = sum(d$FDR < FDR & d$Fold < 0)
      )
    }
  }
}
summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, SUMMARY, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n[A04] Summary (Default vs Background vs Csaw):\n")
print(summary_df)

cat("\n[DONE] A04_Diffbind_Compare_NORMS complete\n")
