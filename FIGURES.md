# Figure manifest

Maps every panel in the article to the exact script that generates it. Folder
names below are content-based (`genotype_validation`, not `Figure1`), so a
future renumbering of the manuscript's figures only means editing this table
— no files move, nothing needs re-running.

Where a panel has no generating script (a photograph, a Sanger trace, a
hand-drawn schematic, an externally supplied table), it's listed as a static
asset, not a pipeline output.

## Genotype / line validation (current manuscript: Figure 1)

| Panel | Script | Output folder |
|---|---|---|
| KMT2D/KDM6A allele quantification | `scripts/genotyping/AR10c_Sane_Mutated_Panel.R` | `results/figures/genotype_validation/` |
| Mutation read-evidence pileups (KMT2D, KDM6A) | `scripts/genotyping/AR10_Mutation_Locus.py` | `results/figures/genotype_validation/` |
| Sanger trace, immunofluorescence images, microscopy panels | *static assets* — imaging/sequencing-facility output, not pipeline-generated | `results/figures/genotype_validation/assets/` |

## RNA-seq transcriptional phenotype (current manuscript: Figure 2)

| Panel | Script | Output folder |
|---|---|---|
| Volcano plots (KDM6A_ko vs WT, KMT2D_Het vs WT), PCA, BMP/WNT signalling barplot, GO BP dotplot, patterning/lineage gene barplot | `scripts/publication_figures/AR07_Publication_Figures.R` | `results/figures/rna_transcriptional_phenotype/` |
| Cross-genotype convergence scatter, up/down Venn diagrams | `scripts/rna/R04_RNA_Overlap.R` | `results/rna/` |
| *(one panel currently has no traced generating script — see note below)* | — | — |

## ATAC–RNA concordance (current manuscript: Figure 3)

| Panel | Script | Output folder |
|---|---|---|
| ATAC volcano plots (KDM6A_ko vs WT, KMT2D_Het vs WT) | `scripts/atac/A05_DESeq2.R` | `results/atac/` |
| DAR genomic annotation | `scripts/integration/AR01_ATAC_RNA_Integration.R` | `results/integration/` |
| Four-way Venn (up/down, both contrasts, both assays) | `scripts/integration/AR04_FourWay_Venn.R` | `results/integration/` |
| ATAC–RNA concordant heatmap | `scripts/integration/AR05_ATAC_RNA_Concordant_Heatmap.R` | `results/integration/` |
| GO enrichment, 103-gene concordant set | `scripts/integration/AR07b_save_103gene_GO_table.R` + `scripts/integration/AR02_GO_Plots.R` | `results/integration/` |

## TF footprinting and regulatory network (current manuscript: Figure 4)

| Panel | Script | Output folder |
|---|---|---|
| TOBIAS binding heatmaps + volcano plots (both contrasts) | `scripts/atac/A07b_TOBIAS_DOWNSTREAM.R` | `results/atac/` |
| CollecTRI TF regulatory network (both contrasts, force-simulation layout) | `scripts/integration/AR06b_Direct_Site_Edges.R` (TF→peak→gene edge table, extracted from an otherwise-superseded script — see `docs/decisions.md`) → `AR11_CollecTRI_fetch_regulons.R` → `AR11d_CollecTRI_dense_with_candidates.R` → `AR11l_ForceSim_network.R` (uses `AR11_network_plot_helpers.R`) | `results/integration/` |
| Tornado plot — both-contrast-significant TFs, footprint vs. background | `scripts/publication_figures/tornado/` (see that folder's own `README.md` for the run order; depends on a public ChIP reanalysis and public ENCODE tracks, fetched by `00_fetch_encode_tracks.sh` / documented in `01_fetch_akiyama_chip.md`) | `results/figures/tf_footprinting_network/` |
| TOBIAS method schematic (`TOBIAS_NETWORK_ORGA.svg`) | *static asset* — hand-drawn diagram, not pipeline output | `results/figures/tf_footprinting_network/assets/` |

## Supplementary

| Panel | Script | Output folder |
|---|---|---|
| ATAC-seq library QC panel | `scripts/atac/A03_QC_ATAC.sh` | `results/atac/` |
| RNA-seq library QC panel | `scripts/rna/R03_RNA_QC.R` | `results/rna/` |
| eSNP/BAF karyotype panels, expression dosage heatmap | `scripts/rna/R08_eSNP_Karyotype.R` | `results/rna/` |
| ATAC–RNA scatter, anterior patterning loci | `scripts/publication_figures/AR07_Publication_Figures.R` | `results/figures/supplementary/` |
| gRNA off-target table, other manuscript-assembly screenshots | *static assets* — wet-lab/manuscript-assembly material, not pipeline output | `results/figures/supplementary/assets/` |

## Notes

- **One Figure 2 panel currently has no traced generating script.** It isn't
  reproduced by any script in this repository as of this writing; if you need
  it, ask the corresponding author rather than assuming it can be regenerated
  from the pipeline as-is.
- **The tornado plot's background-region panel is a documented
  reconstruction, not a byte-identical rerun** of the original: the exact
  region-selection logic used for the published figure was never scripted and
  the intermediate files are gone. `scripts/publication_figures/tornado/03_define_background_regions.R`
  implements a sound, clearly commented standalone method instead — see that
  script's header for exactly what it does and why.
- **The tornado plot's underlying "high-confidence WT-bound sites" input
  (`all_WT_bound_highconf.bed`) also has no recoverable generating script.**
  A genuine attempt was made to reverse-engineer the exact per-TF selection
  rule against the still-existing merged file (its score column was traced to
  `KDM6A_ko_vs_WT`'s BINDetect `WT_score`, and several selection hypotheses —
  a fixed score threshold, top-N-by-score, one-site-per-peak — were tested
  and ruled out), but the actual rule could not be independently derived.
  Rather than ship a guess that might silently change the figure,
  `scripts/publication_figures/tornado/02_define_bound_highconf_canon.sh`
  documents this honestly as an open gap.
- **Every ATAC differential-accessibility result in this repo uses the loess
  normalization arm** (`scripts/atac/A04b_Csaw_Loess_Norm.R`), the one
  reported in the paper. An earlier TMM arm was evaluated and rejected; see
  `docs/A04b_normalization_methodology.md` for the full comparison — it is not
  used by any figure here.
- **chromVAR was evaluated as a cross-check for TF binding calls and
  rejected** (wildly asymmetric power between contrasts) in favor of TOBIAS's
  own binding-significance flag; see `scripts/atac/A07b_TOBIAS_DOWNSTREAM.R`
  for the reasoning. Not part of any figure here.
- **ANANSE was evaluated as a network-inference method and dropped** in favor
  of the CollecTRI-based network (ATAC-only signal was judged too weak on its
  own). Not present in this repository at all.
