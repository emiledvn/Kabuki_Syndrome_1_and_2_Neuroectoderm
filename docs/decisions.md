# Methodology choices

This repository contains only the single analysis path actually used for the
reported results. Where multiple methods were genuinely evaluated, this page
records which one was used and why, and where the full comparison lives.

## Consensus peak set and TMM comparison script: restored, not excluded

`scripts/atac/A04_Diffbind_Compare_NORMS.R` — which builds the DiffBind
consensus peak set from the per-sample MACS3 peaks and runs the
Default/Background/Csaw-TMM comparison arms — was missed by the initial
curation pass and absent from this repository even though `A04b`'s loess
arm (below) depends directly on its `Csaw_Norm/diffbind_analyzed.rds`
output. Unlike the genuinely-superseded arms documented elsewhere on this
page, this script is a live, still-necessary dependency, not a rejected
alternative; it has been added back to `scripts/atac/` and wired into
`EXECUTE_PIPELINE.sh` ahead of `A04b`.

## ATAC-seq normalization: loess, not TMM

Standard TMM normalization produced opposite-direction, strongly asymmetric
differential-accessibility calls between the two genotype contrasts. A
non-linear (loess) alternative was evaluated instead and adopted after a
background-bin-vs-peak diagnostic and a WT-vs-WT null check. Full comparison
and numbers: `docs/A04b_normalization_methodology.md`. Only the loess arm is
present in this repository; TMM is described in that document as the
comparison baseline.

## TF binding: TOBIAS, not chromVAR, as the primary call

chromVAR was evaluated as an independent, per-sample cross-check on TF
binding activity. Its test showed wildly asymmetric statistical power between
the two contrasts (2 of 612 motifs significant in KDM6A_ko_vs_WT vs. 407 of
612 in KMT2D_Het_vs_WT — a >100x difference), which is not a plausible
biological result and points to a power/assumptions problem with the method
on this dataset rather than a real effect. TOBIAS's own footprint-based
binding-significance call was used instead throughout. See
`scripts/atac/A07b_TOBIAS_DOWNSTREAM.R` for the full reasoning. chromVAR is
not included in this repository.

## TF regulatory network: CollecTRI, not ANANSE

ANANSE infers regulatory networks from chromatin accessibility alone. Applied
to this dataset, the resulting networks were judged too weak to support the
specific TF-target claims made in the manuscript (ATAC signal alone,
without an orthogonal expression-based prior, was not sufficiently
discriminating). The CollecTRI-based network — which combines TOBIAS
footprint calls with curated TF-target regulon information — was used
instead. ANANSE is not included in this repository.

## Figure 4 network's edge table: extracted, not excluded wholesale

The working repository's TF→peak→gene edge computation (TOBIAS per-site
footprint calls intersected with DAR/DEG-anchored gene links) lived in the
same script as an igraph-based network visualization that was superseded by
the CollecTRI-based network below. Only the visualization was dropped; the
edge-table computation itself is still the CollecTRI network's actual input
and has been kept as its own script, `scripts/integration/AR06b_Direct_Site_Edges.R`.

## Figure 4 network layout: pinned to the exact published version

The CollecTRI network's force-simulation layout script
(`scripts/integration/AR11l_ForceSim_network.R`) went through several
cosmetic revisions after the figure was finalized (larger labels, stronger
label decluttering). This repository pins that script to the exact parameter
version that produced the published PDF (internally versioned "v14" in the
working repository's history), not the latest cosmetic revision. The
network's underlying node layout is unaffected either way — deterministic,
fixed seed — only label rendering changed between versions.

## Figure 4 tornado plot: implementation notes

- The background-region set for the footprint-vs-background panel is defined
  by `scripts/publication_figures/tornado/03_define_background_regions.R`
  (consensus ATAC peaks with no TOBIAS footprint, fixed-seed sampling).
- The upstream "high-confidence WT-bound sites" file
  (`all_WT_bound_highconf.bed`) that both tornado panels depend on is an
  input to `scripts/publication_figures/tornado/02_define_bound_highconf_canon.sh`
  — see that script's header for what it expects.

See `FIGURES.md` for the full figure-to-script mapping.
