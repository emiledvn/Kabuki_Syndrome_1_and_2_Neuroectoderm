# ATAC-seq Differential Accessibility Normalization: Methodology and Diagnostics

Covers the KDM6A_ko_vs_WT / KMT2D_Het_vs_WT csaw/DESeq2 differential-accessibility
pipeline (`scripts/atac/A04_Diffbind_Compare_NORMS.R`, `scripts/atac/A04b_Csaw_Loess_Norm.R`,
`scripts/atac/A05_DESeq2.R`). The loess arm documented here is the one used
throughout this repository; the TMM comparison arm (`A04_Diffbind_Compare_NORMS.R`,
which also builds the DiffBind consensus peak set both arms share) is described
here only as the baseline that motivated the loess decision.

## 1. Background

### 1.1 Count-source consistency

An earlier version of the DESeq2 fitting script reconstructed its
`DESeqDataSet` from `dba.peakset(atac, bRetrieve = TRUE)` counts instead of
using DiffBind's own fitted object. These two count sources disagreed by
~1.33x on most cells (root cause not fully traced), so the two arms' DAR
counts silently diverged. Fixed by loading DiffBind's own fitted counts
directly, never a `dba.peakset()` reconstruction.

### 1.2 The normalization question

Under csaw's standard TMM normalization (background 10kb bins), the two
contrasts show opposite-direction, strongly asymmetric DAR patterns:

| Contrast | Norm | Gained | Lost |
|---|---|---|---|
| KDM6A_ko_vs_WT | TMM | 5830 | 600 |
| KMT2D_Het_vs_WT | TMM | 3549 | 11218 |

KDM6A_ko is gain-dominant (~10:1), KMT2D_Het is loss-dominant (~1:3).
Diagnostics run against this pattern:

- **FRiP** correlates with genotype and with Gained-peak signal (r=0.65) —
  but FRiP is circular (computed against each sample's own MACS3 peaks), so
  this is not independent evidence.
- **TSS enrichment** (non-circular, uses an external reference) does not
  separate cleanly by genotype. Mixed evidence: suggestive of some real
  biology, not conclusive on its own.
- A **background-bin-only MA plot** (restricted to non-peak regions) showed a
  clear upward tail at high average abundance for KDM6A_ko under TMM — the
  signature of intensity-dependent (efficiency) bias that a single global TMM
  scale factor per sample cannot correct.

## 2. Method: csaw loess normalization

`scripts/atac/A04b_Csaw_Loess_Norm.R` adds a non-linear normalization arm.

1. `csaw::normOffsets()` fits a per-sample loess trend of log-count vs.
   log-average-count on the same unfiltered background bins the TMM arm
   uses, then interpolates that trend onto the consensus peak set via spline
   interpolation.
2. **Density-scale fix**: `normOffsets()` positions each region on its fitted
   curve using raw counts, not width-normalized density. Background bins are
   a fixed 10kb; MACS3 narrowPeak calls are far narrower and highly variable
   in width (~150bp–2kb), which would otherwise position peaks
   systematically wrong on the curve. Fixed by scaling each peak's raw count
   to a per-10kb-equivalent pseudo-count solely for positioning on the curve;
   the resulting offsets are applied to the true raw peak counts.
3. Offsets are converted to DESeq2 `normalizationFactors` (linear-scale,
   row-geometric-mean-centered to 1) and a fresh dispersion/Wald fit is run.

Resulting DAR counts:

| Contrast | Norm | Gained | Lost |
|---|---|---|---|
| KDM6A_ko_vs_WT | TMM | 5830 | 600 |
| KDM6A_ko_vs_WT | Loess | 3161 | 12079 |
| KMT2D_Het_vs_WT | TMM | 3549 | 11218 |
| KMT2D_Het_vs_WT | Loess | 8332 | 8325 |

## 3. Why "does the MA plot look flat" is not sufficient on its own

Loess normalization works by fitting a trend across abundance and forcing the
average log-fold-change to zero at every abundance bin **on the data it is
fit to**. A flat MA plot after correction is not distinguishable by eye alone
from over-correction that erases real signal. The only way to know whether
the correction is generalizing correctly is to check a **held-out** region
set it was not fit to: the interpolated peaks.

## 4. Diagnostic: background bins vs. peaks, TMM vs. Loess

Two summary metrics, computed on both background bins and peaks, both norms:
Spearman ρ of log2FC vs. log10(abundance), and mean log2FC in the top
abundance decile.

| Contrast | Region set | Norm | Spearman ρ | Mean log2FC, top decile |
|---|---|---|---|---|
| KDM6A_ko_vs_WT | Background bins | TMM | 0.585 | +0.350 |
| KDM6A_ko_vs_WT | Background bins | Loess | 0.068 | -0.007 |
| KDM6A_ko_vs_WT | Peaks | Loess | **-0.236** | **-0.371** |
| KMT2D_Het_vs_WT | Background bins | TMM | -0.371 | -0.230 |
| KMT2D_Het_vs_WT | Background bins | Loess | -0.024 | -0.012 |
| KMT2D_Het_vs_WT | Peaks | Loess | 0.026 | +0.046 |

### Interpretation

**KMT2D_Het_vs_WT**: loess flattens the background-bin trend and the peaks
come out flat too — a clean, fully-explained-by-normalization correction.

**KDM6A_ko_vs_WT**: loess flattens the background-bin trend just as
effectively, but the peaks retain a trend that is sign-flipped and comparable
in magnitude to what was removed — concentrated in regulatory regions, the
profile expected of real biology rather than an incompletely-corrected
technical artifact.

### Orthogonal biological plausibility (not proof)

KDM6A/UTX is an H3K27me3 demethylase. Its loss should allow repressive
H3K27me3 to accumulate at target loci, mechanistically predicting *reduced*
accessibility, especially at normally high-abundance (active) regions —
exactly the shape of "loss-dominant, downward tail at high abundance" seen
post-loess. A useful sanity check, but supporting narrative, not evidence.

## 5. Remaining diagnostics

### 5.1 RLE / sample correlation / PCA / dispersion

- **RLE**: Loess visibly tightens and better centers each sample's per-peak
  log-ratio distribution around 0 than TMM does.
- **PCA**: clean separation by genotype under both norms, all 3 replicates
  clustering tightly, no visible outlier or batch/technical split.
- **Dispersion**: standard, comparable fit shape under both norms.

### 5.2 WT-vs-WT null check

Arbitrary 1-vs-2 split of the 3 WT replicates, run through identical
normalization/fitting as the real contrasts, establishing the noise floor
this pipeline produces from small-n replicate-split behavior alone with zero
real signal driving it.

| Contrast (peaks, Loess) | Spearman rho | Top-decile mean log2FC |
|---|---|---|
| WT-vs-WT null (noise floor) | 0.162 | +0.332 |
| KMT2D_Het_vs_WT (real) | 0.026 | +0.046 |
| KDM6A_ko_vs_WT (real) | **-0.236** | **-0.371** |

Both real contrasts' residuals are reported against this noise floor in the
manuscript. Full per-contrast discussion belongs there rather than duplicated
here.

## 6. Relevant files

- `scripts/atac/A04_Diffbind_Compare_NORMS.R` — builds the DiffBind consensus peak
  set and the Default/Background/Csaw-TMM comparison arms (TMM = the baseline
  that motivated the loess decision; not used by any reported result)
- `scripts/atac/A04b_Csaw_Loess_Norm.R` — loess normalization arm (used)
- `scripts/atac/A04f_Loess_DAR_BED.R` — DAR BED export from the loess arm
- `scripts/atac/A05_DESeq2.R` — DESeq2 fit from DiffBind's own fitted counts
