# Methods

Sourced from the pipeline's own config/scripts/environment files, not recalled
from memory — values below are cited against their source. Sections needing
detail this repository doesn't have on record are marked **[NEEDS INPUT]**
rather than guessed.

## Cell lines and experimental design

Three genotypes were compared: wild-type (WT) and two independent
Kabuki-syndrome-model lines, **KMT2D_Het** (heterozygous frameshift,
c.1301_1302insA / WT, p.L434fs / WT — KS1) and **KDM6A_ko**
(nonsense/knockout, c.306_318delinsTAAgaggatcc, p.N103\* — KS2), each
represented by 3 independent clones (source: `config/sample_metadata.csv` —
clone IDs WT_1/2/3, KMT2D_11_A/12_i/14_B, KDM6A_e2/e9/e23). RNA integrity
numbers (RIN) for the RNA-seq samples ranged 7.3–10.0.
**[NEEDS INPUT: differentiation protocol / stage reference — Day 8
neuroectoderm per project context, but the specific published
protocol/citation isn't recorded in any pipeline file]**.

RNA-seq: 2×151bp paired-end. ATAC-seq: 2×76bp paired-end (both confirmed via
FastQC/MultiQC). **[NEEDS INPUT: sequencing platform/instrument, library prep
kit for both assays]**.

## Reference genome and annotation

GENCODE release 47, GRCh38 primary assembly (source:
`config/pipeline_config.yaml`), ENSG gene IDs, chr-prefixed hg38 coordinates.
Genome and GTF fetched directly from EBI's GENCODE mirror with MD5-verified
downloads.

## RNA-seq processing and differential expression

Reads were processed with **nf-core/rnaseq v3.26.0** (STAR alignment + Salmon
quantification, standard pipeline defaults). Gene-level counts were
pre-filtered to genes with ≥10 total counts across all samples before
differential testing. Differential expression was fit per contrast (KMT2D_Het
vs WT; KDM6A_ko vs WT) with **DESeq2 v1.50.2**, design `~ genotype`, Wald
test, BH-adjusted p-values. **No log2FC shrinkage** (apeglm/ashr) was applied
to the reported fold-changes: shrinkage was evaluated against this dataset and
rejected for DEG calling because it over-suppressed low-baseMean
neuroectoderm-patterning/pluripotency TFs at n=3/genotype (e.g. EMX1, LHX2,
SIX3, FOXG1, POU5F1, HOXB1 collapsed from raw log2FC of -2.1 to -4.4 down to
-0.05 to -0.24 post-shrinkage) — these are exactly the genes where a modest
absolute expression change is expected to be functionally consequential, so
shrinkage toward zero would suppress real signal rather than correct noise.
Raw log2FoldChange + padj was used throughout for DEG calling; any rank-based
analysis (GSEA) ranks by the Wald statistic, never by LFC, for the same
reason.

Significance threshold throughout (RNA and ATAC alike, see below): **padj <
0.05 and |log2FC| > 1** (2-fold) — a single threshold, not a two-tier
relaxed/strict system.

GO Biological Process enrichment: **clusterProfiler v4.18.4** `enrichGO()`,
BH-adjusted, explicit background universe = all genes passing the count
pre-filter (not clusterProfiler's whole-genome default). GSEA: **fgsea
v1.36.2**, ranked by DESeq2's Wald statistic, against GO Biological Process,
MSigDB Hallmark, and MSigDB Reactome (via **msigdbr v26.1.0**).

## ATAC-seq processing, QC, and differential accessibility

Reads were processed with **nf-core/atacseq v2.1.2**; peaks called with
**MACS3 v3.0.4** (narrow peak mode). QC: FRiP, peak count, and library
complexity from ataqv's own output where available; TSS enrichment was
recomputed directly against a curated RefSeq TSS set (not GENCODE's), since
ENCODE's published TSS-enrichment cutoff table (<5 Concerning, 5–7
Acceptable, >7 Ideal) is calibrated against RefSeq's curated TSS set.

**Normalization**: csaw/DiffBind background-bin TMM normalization was
compared against a non-linear (Loess, `csaw::normOffsets`) alternative (full
methodology and diagnostics in `docs/A04b_normalization_methodology.md`).
Loess was selected as the reported arm after: (1) background-bin-restricted
MA plots showing a real abundance-dependent trend under TMM not present under
Loess; (2) a WT-vs-WT null check establishing the noise floor for this
pipeline's small-n replicate-split behavior; (3) an extrapolation/sparse-
support check confirming the correction isn't fabricated in a data-sparse
abundance region. TMM is retained in `docs/A04b_normalization_methodology.md`
only as the comparison baseline that motivated this decision — it is not
used by any reported result.

Differential accessibility: **DESeq2 v1.50.2** fit directly on DiffBind's own
fitted counts (never reconstructed from `dba.peakset()`, which was found to
diverge from DiffBind's own fitted counts by ~1.33× in an earlier version of
this pipeline). Same threshold as RNA: padj < 0.05, |log2FC| > 1.

## Cross-genotype convergence (RNA and ATAC)

For both DEGs and DARs, overlap between the two genotypes' significant sets
was quantified two ways: **Fisher's exact test odds ratio** (kept only as an
effect-size summary, since Fisher's own p-value assumes gene/region-
independent Bernoulli trials — violated by co-regulated genes/co-accessible
regions — making it markedly anti-conservative at this scale), and a
**baseMean-binned permutation null** (10,000 reps, 20 bins, fixed seed for
reproducibility) as the actual significance test. Genome-wide Spearman
correlation (log2FC vs log2FC across all tested genes/regions, not just the
significant ones) is reported as a complementary, threshold-free convergence
measure.

## TF/motif-level analyses

**TOBIAS v0.17.3** (`BINDetect`) measured per-site footprint depth at
predicted binding sites, pooled pseudo-bulk BAM per genotype. Differential
binding threshold: p < 0.05 and |binding score change| > 0.1. Binding-
direction sign convention: positive = higher in WT, negative = higher in the
mutant.

**chromVAR** (per-sample motif accessibility deviation) was evaluated as an
independent cross-check and rejected: its t-test showed wildly asymmetric
power between contrasts (2/612 significant motifs in KDM6A_ko_vs_WT vs
407/612 in KMT2D_Het_vs_WT, a >100× difference) — see
`scripts/atac/A07b_TOBIAS_DOWNSTREAM.R` for the full reasoning. Not used in
any reported result; excluded from this repository's automated pipeline.

**ANANSE** (network inference from ATAC accessibility alone) was evaluated
for the TF regulatory network and dropped in favor of the CollecTRI-based
network below — ATAC-only signal was judged too weak on its own. Not present
in this repository.

## ATAC-RNA integration

Peak-to-gene annotation: **rGREAT v2.12.2**, basal+extension regulatory-
domain model (basal domain -5kb/+1kb from TSS, extension up to 1Mb,
truncated at neighboring basal domains). Applied only to called GAINED/LOST
DAR peaks, not the full peak set.

Per-gene ATAC summary: the single strongest assigned peak (max |log2FC|),
with both fold-change and FDR taken from that same peak — not a
sum-of-log-fold-changes/min-FDR aggregation across all assigned peaks (an
earlier version summed log-fold-changes, which implicitly models peak
effects as multiplying together and inflates apparent effect size for genes
assigned more peaks by GREAT's model — itself correlated with local gene
density, not independent regulatory element count).

ATAC-RNA concordance was classified per gene as Concordant (both
ATAC- and RNA-significant, same direction), Discordant (both significant,
opposite direction), ATAC-only, RNA-only, or Not significant.

## TF regulatory network

A peak-resolved TF→target-gene network was built from **TOBIAS's own
per-site footprint calls** (not monaLisa's set-level enrichment, which only
says a motif is enriched somewhere across an entire GAINED/LOST peak set, not
at which specific peak). Each individually-scored TF binding site was
intersected (via genomic-range overlap) against the full peak-resolved
DAR-to-gene link table, giving specific TF→peak→gene edges, collapsed to one
edge per (TF, gene) pair and tiered by independent corroborating evidence.
RNA evidence is a corroborating tier on top of this ATAC/motif-derived
backbone, not the signal that defines an edge.

## Pipeline management and reproducibility

The full analysis is implemented as a single version-controlled pipeline
rather than a set of interactive scripts, orchestrated end-to-end by one
driver (`scripts/EXECUTE_PIPELINE.sh`) that runs all steps in dependency
order across five stages: ATAC-seq processing, RNA-seq processing,
RNA–ATAC integration, genotype/allele validation, and manuscript figure
assembly. Every step is self-checkpointing (skipped if its declared output
already exists), so an interrupted run resumes rather than restarts, and any
individual step or downstream sub-chain can be re-run in isolation (e.g.
after a threshold change) via `--from`/`--only` flags without re-executing
the full pipeline. All statistical thresholds, reference-genome paths,
genotype contrasts, and compute-resource limits are centralized in a single
config file (`config/pipeline_config.yaml`) that every script reads — no
threshold or path is hand-edited inside an individual script.

Software is version-pinned throughout. nf-core pipelines (rnaseq v3.26.0,
atacseq v2.1.2, rnavar v1.3.0) ran under Nextflow v25.10.2 with
Apptainer v1.4.5-containerized processes. Custom R and Python analyses ran
in five dedicated conda environments, each specified in a checked-in
`environments/*.yml` file re-exported after every package addition and
verified against the tool versions actually present at runtime
(`environments/00_Environment_setup_versions.sh`, logged to
`results/Session_info/`). The reference genome/annotation (GENCODE release
47, GRCh38) and other externally-sourced references (RefSeq TSS
coordinates, GATK known-sites bundles, Roadmap ChromHMM segmentations) were
fetched from their canonical source URLs, recorded in
`config/pipeline_config.yaml` alongside MD5 checksums where feasible.
Stochastic procedures (the baseMean-binned cross-genotype permutation test;
tornado-plot background-region sampling) used fixed random seeds.

Every analytical choice made among multiple candidate methods (loess vs.
TMM accessibility normalization; TOBIAS vs. chromVAR TF-binding calls;
CollecTRI vs. ANANSE network inference) is documented alongside the
rejected alternative and the diagnostic that motivated the decision
(`docs/decisions.md`, `docs/A04b_normalization_methodology.md`), and every
manuscript panel is mapped to the exact script and parameter version that
produced it (`FIGURES.md`), independent of the manuscript's own figure
numbering.

### Reproducing the analysis

The companion repository (see Code availability) contains everything needed
to reproduce all analyses from raw FASTQs: pinned environment
specifications, the config file, and the orchestration script. In brief:
create the five conda environments from `environments/*.yml`, place raw
FASTQs under `data/fastq_rna/` and `data/fastq_atac/` matching
`input/*_samplesheet.csv`, and run `bash scripts/EXECUTE_PIPELINE.sh` (see
the repository README for the full walkthrough and for the one manual
external step, a public ChIP re-analysis feeding the Figure 4 tornado
plot).
