#!/usr/bin/env Rscript
# 03_define_background_regions.R -- defines the "background" region set for
# Panel B (footprint-vs-background binding map): consensus ATAC peaks that
# have NO TOBIAS-called footprint, downsampled to match the size of the
# bound_highconf_canon set.
#
# NOTE ON PROVENANCE: the ORIGINAL background set used for the published
# figure (bg_fixed.bed) was defined in an interactive exploration session
# and never captured as a script -- its exact selection logic is lost. This
# script defines a fresh, sound, and fully documented replacement using the
# same conceptual definition ("consensus peak, not footprinted"), but it is
# NOT guaranteed to reproduce the original figure pixel-for-pixel. Treat any
# rerun of this panel as scientifically equivalent, not a byte-identical
# reproduction, and say so explicitly wherever this figure is described.
#
# Method: consensus ATAC peaks (data/tobias_output/consensus_peaks.bed, the
# same set TOBIAS BINDetect itself scans) minus any peak overlapping a
# TOBIAS high-confidence WT-bound site (bound_highconf_canon.bed from
# 02_define_bound_highconf_canon.sh), then downsampled with a fixed seed to
# match bound_highconf_canon's region count.
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

REPO_ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), "..", "..", ".."))
setwd(REPO_ROOT)

SEED <- 42
CONSENSUS_BED <- "data/tobias_output/consensus_peaks.bed"
BOUND_BED     <- "results/atac/tornado_intermediate/bound_highconf_canon.bed"
OUT_DIR       <- "results/atac/tornado_intermediate"
OUT           <- file.path(OUT_DIR, "bg_regions.bed")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(CONSENSUS_BED)) {
  stop(sprintf("[03] %s not found -- run scripts/atac/A06_TOBIAS.sh first (it builds the consensus peak set).", CONSENSUS_BED))
}
if (!file.exists(BOUND_BED)) {
  stop(sprintf("[03] %s not found -- run 02_define_bound_highconf_canon.sh first.", BOUND_BED))
}

consensus <- import(CONSENSUS_BED, format = "bed")
bound     <- import(BOUND_BED, format = "bed")

unbound <- consensus[!overlapsAny(consensus, bound)]
cat(sprintf("[03] %d consensus peaks, %d overlap a footprint, %d unbound candidates\n",
            length(consensus), length(consensus) - length(unbound), length(unbound)))

n_target <- length(unique(bound))
n_target <- min(n_target, length(unbound))

set.seed(SEED)
bg <- unbound[sample(length(unbound), n_target)]
bg <- sort(bg)

export(bg, OUT, format = "bed")
cat(sprintf("[03] wrote %d background regions (seed=%d) -> %s\n", length(bg), SEED, OUT))
