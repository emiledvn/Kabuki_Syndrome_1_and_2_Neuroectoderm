# Kabuki Syndrome (KMT2D/KDM6A) Neuroectoderm 

RNA-seq and ATAC-seq analysis of Day 8 neuroectoderm differentiation, comparing
wild-type (WT) against two independent Kabuki-syndrome-model genotypes:
**KMT2D_Het** (heterozygous frameshift, KS1) and **KDM6A_ko** (knockout, KS2),
each represented by 3 independent clones. Two nf-core pipelines (RNA-seq,
ATAC-seq) feed a shared set of downstream differential expression/accessibility,
transcription-factor footprinting (TOBIAS), and RNA–ATAC integration steps.

We also developed a Network approach Based on TOBIAS and RNA-seq

## Reproducibility

Every tool/package version is pinned (`environments/*.yml` conda
environments, exact container/Nextflow versions, a specific GENCODE release
with checksums), stochastic steps are seeded, and figure-generating scripts
are pinned to the exact parameter version that produced the published panel
(see `FIGURES.md`). This was a genuine best effort at full end-to-end
reproducibility, not a guarantee — if you hit a discrepancy, please open an
issue with your environment and command.

## Repository layout

```
scripts/
├── atac/                 ATAC-seq processing, QC, differential accessibility, TOBIAS footprinting
├── rna/                  RNA-seq processing, QC, differential expression, karyotype/dosage checks
├── integration/          ATAC-RNA concordance, GO enrichment, CollecTRI TF regulatory network
├── genotyping/           Allele quantification and mutation-evidence panels (Figure 1)
├── publication_figures/  Final manuscript panel assembly (see FIGURES.md)
└── EXECUTE_PIPELINE.sh   Runs every stage in order

results/
├── atac/, rna/, integration/   Full pipeline output (not just what made the paper)
└── figures/                    The curated panels that appear in the manuscript

config/            Pipeline parameters, genotype contrasts, GEO sample table
environments/      Pinned conda environment specifications
docs/              Methods writeups and normalization/methodology decision records
input/             nf-core samplesheets (relative paths — portable across machines)
```

## Quickstart

1. Create the pinned environments (one per pipeline stage):
   ```bash
   for env in environments/*.yml; do conda env create -f "$env"; done
   ```
2. Place raw FASTQs under `data/fastq_rna/` and `data/fastq_atac/` (see
   **Data availability** below) — filenames must match `input/*_samplesheet.csv`.
3. Run the pipeline:
   ```bash
   bash scripts/EXECUTE_PIPELINE.sh
   ```
   Each stage is self-checkpointing; re-running after an interruption resumes
   rather than restarting. The final step (`scripts/publication_figures/`)
   assembles the exact manuscript panels; the tornado plot in
   `scripts/publication_figures/tornado/` additionally needs a public ChIP
   reanalysis and ENCODE tracks — see that folder's own `README.md`.

## Data availability

Raw and processed data are deposited at GEO under accession **[pending —
GSE accession to be added once assigned]**. Sample-to-accession mapping is in
`config/sample_metadata.csv`. Two external, public datasets are also used
(TF-network and tornado-plot panels only, not the core RNA/ATAC analysis):

- ChIP-seq reanalysis of Akiyama et al. (bioRxiv 2025.07.03.663017 /
  *Development* 2026, PMID 41906541), GEO **GSE301295**.
- ENCODE H1-hESC histone ChIP tracks (H3K4me1, H3K27ac, H3K27me3) — accessions
  in `scripts/publication_figures/tornado/00_fetch_encode_tracks.sh`.

## Figures

`FIGURES.md` maps every manuscript panel to the exact script that produces it,
independent of the manuscript's current figure numbering (which may still
change during review).

## License

MIT — see `LICENSE`.

## Citation

**[citation to be added on publication]**
