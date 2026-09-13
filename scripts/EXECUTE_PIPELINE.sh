#!/usr/bin/env bash
# EXECUTE_PIPELINE.sh -- runs the full pipeline end to end, in order. Every step
# is self-checkpointing (skips if its declared output already exists), so
# rerunning this after an interruption just resumes where it left off.
#
# Usage:
#   bash scripts/EXECUTE_PIPELINE.sh                # full run / resume
#   bash scripts/EXECUTE_PIPELINE.sh --from A05      # rerun A05 onward (e.g. after
#                                                     #   changing a threshold)
#   bash scripts/EXECUTE_PIPELINE.sh --only A04b     # rerun just one step
#
# Steps run strictly in the order below by default. A01-A09 and R01-R04 are
# independent arms that CAN run concurrently (we did so manually for the actual
# overnight launch) -- to do that here, invoke this script twice with --only
# targeting different steps in separate screens, or call the individual step
# scripts directly (always still supported on its own).
#
# Five stages: ATAC (A0*) -> RNA (R0*) -> ATAC/RNA integration (AR0*) ->
# genotyping/allele validation (AR10*) -> publication figure assembly (AR07
# top-level + the tornado/ chain).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SCREEN_SESSION="ks_1_2_pipeline"
if [[ -z "${STY:-}" ]]; then
    echo "[EXECUTE_PIPELINE] Not inside screen -- relaunching in screen session '${SCREEN_SESSION}' (reattach with: screen -r ${SCREEN_SESSION})."
    mkdir -p logs
    exec screen -L -Logfile "logs/EXECUTE_PIPELINE_screen_$(date +%Y%m%d_%H%M%S).log" -S "${SCREEN_SESSION}" bash "$0" "$@"
fi

# R-based steps need the ks_1_2_r conda env explicitly -- a bare `Rscript` here
# would silently resolve to the system R (/usr/local/bin/Rscript, no pinned
# packages installed) instead, since this script itself is never conda-activated.
# `conda run --no-capture-output` activates just for that one call and streams
# stdout/stderr live (without --no-capture-output, conda run buffers/redirects
# output, which breaks tee/screen logging downstream).
CONDA_BASE="${KS12_CONDA_BASE:-/storage/volume01/emile/miniconda3}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
RSCRIPT_KS12="conda run --no-capture-output -n ks_1_2_r Rscript"
PYTHON_VIZ="conda run --no-capture-output -n viz_env python"

declare -A STEP_CMD=(
    [00]="bash environments/00_Environment_setup_versions.sh"
    [A01]="bash scripts/atac/A01_RUN_nfcore_ATAC.sh"
    [A02]="bash scripts/atac/A02_MACS3_narrow.sh"
    [A03]="bash scripts/atac/A03_QC_ATAC.sh"
    [A04b]="${RSCRIPT_KS12} scripts/atac/A04b_Csaw_Loess_Norm.R"
    [A04f]="${RSCRIPT_KS12} scripts/atac/A04f_Loess_DAR_BED.R"
    [A05]="${RSCRIPT_KS12} scripts/atac/A05_DESeq2.R"
    [A06]="bash scripts/atac/A06_TOBIAS.sh"
    [A07b]="${RSCRIPT_KS12} scripts/atac/A07b_TOBIAS_DOWNSTREAM.R"
    [A08]="bash scripts/atac/A08_ATAC_Viz.sh"
    [A09]="${RSCRIPT_KS12} scripts/atac/A09_Motif_Enrichment.R"
    [R01]="bash scripts/rna/R01_RUN_nfcore_RNA.sh"
    [R07]="bash scripts/rna/R07_RUN_nfcore_RNAVAR.sh"
    [R02]="${RSCRIPT_KS12} scripts/rna/R02_DESeq2.R"
    [R03]="${RSCRIPT_KS12} scripts/rna/R03_RNA_QC.R"
    [R04]="${RSCRIPT_KS12} scripts/rna/R04_RNA_Overlap.R"
    [R08]="${RSCRIPT_KS12} scripts/rna/R08_eSNP_Karyotype.R"
    [AR01]="${RSCRIPT_KS12} scripts/integration/AR01_ATAC_RNA_Integration.R"
    [AR11]="${RSCRIPT_KS12} scripts/integration/AR11_CollecTRI_fetch_regulons.R"
    # AR06b extracts just the direct_site_edges computation from the source
    # repo's AR06_TF_Regulatory_Network.R (that script's own network
    # visualization was dropped, superseded by the CollecTRI/AR11 network
    # actually used in the paper -- but this edge-table computation is a
    # separate, still-used dependency of AR11d). Needs A06 (TOBIAS) + R02
    # (DESeq2), both already run earlier above.
    [AR06b]="${RSCRIPT_KS12} scripts/integration/AR06b_Direct_Site_Edges.R"
    [AR11d]="${RSCRIPT_KS12} scripts/integration/AR11d_CollecTRI_dense_with_candidates.R"
    [AR04]="${RSCRIPT_KS12} scripts/integration/AR04_FourWay_Venn.R"
    [AR05]="${RSCRIPT_KS12} scripts/integration/AR05_ATAC_RNA_Concordant_Heatmap.R"
    [AR07b]="${RSCRIPT_KS12} scripts/integration/AR07b_save_103gene_GO_table.R"
    [AR02]="${RSCRIPT_KS12} scripts/integration/AR02_GO_Plots.R"
    [AR03]="${RSCRIPT_KS12} scripts/atac/AR03_ATAC_Overlap.R"
    [R05]="${RSCRIPT_KS12} scripts/rna/R05_fgsea.R"
    [R06]="${RSCRIPT_KS12} scripts/rna/R06_Generate_Report.R"
    [AR10_KMT2D]="${PYTHON_VIZ} scripts/genotyping/AR10_Mutation_Locus.py --locus KMT2D"
    [AR10_KDM6A]="${PYTHON_VIZ} scripts/genotyping/AR10_Mutation_Locus.py --locus KDM6A"
    [AR10c]="${RSCRIPT_KS12} scripts/genotyping/AR10c_Sane_Mutated_Panel.R"
    [AR07_pub]="${RSCRIPT_KS12} scripts/publication_figures/AR07_Publication_Figures.R"
)
# AR02 reads AR04's GO table (AR04_GO_venn4_intersections.csv) -- AR04 (and
# AR01, which AR04 itself needs) must run before AR02, not after.
# R07 (eSNP-Karyotyping's variant-calling arm) is independent of the DESeq2 arm --
# only needs raw fastqs + reference -- so it's slotted right after R01. R08 (the
# BAF/dosage analysis) needs both R07's VCFs and R02's VST counts, so it's slotted
# after R05, once both are guaranteed done.
# AR07b needs AR05's concordant-heatmap table, so it runs right after AR05.
# The genotyping/allele-validation steps (AR10*) and the final publication-figure
# assembly (AR07_pub) are new stages added for the publication repo -- the
# source repo ran these by hand rather than through this orchestration script.
# The tornado/ chain (Figure4's TF-binding tornado plot) is intentionally NOT in
# STEP_ORDER: 01_fetch_akiyama_chip.md is a manual/external ChIP-realignment
# step, so that whole chain is run by hand, in order, after this script
# completes -- see scripts/publication_figures/tornado/README.md.
STEP_ORDER=(00 A01 A02 A03 A04b A04f A05 A06 A07b A08 A09 R01 R07 R02 R03 R04 R05 R08 AR01 AR11 AR06b AR11d AR04 AR05 AR07b AR02 AR03 R06 AR10_KMT2D AR10_KDM6A AR10c AR07_pub)

FROM=""
ONLY=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from) FROM="$2"; shift 2 ;;
        --only) ONLY="$2"; shift 2 ;;
        *) echo "[EXECUTE_PIPELINE] ERROR: unknown argument: $1"; exit 1 ;;
    esac
done

is_valid_step() { [[ -n "${STEP_CMD[$1]:-}" ]]; }

if [[ -n "$ONLY" ]]; then
    is_valid_step "$ONLY" || { echo "[EXECUTE_PIPELINE] ERROR: unknown step '$ONLY'. Valid: ${STEP_ORDER[*]}"; exit 1; }
    RUN_STEPS=("$ONLY")
elif [[ -n "$FROM" ]]; then
    is_valid_step "$FROM" || { echo "[EXECUTE_PIPELINE] ERROR: unknown step '$FROM'. Valid: ${STEP_ORDER[*]}"; exit 1; }
    RUN_STEPS=()
    FOUND=0
    for s in "${STEP_ORDER[@]}"; do
        [[ "$s" == "$FROM" ]] && FOUND=1
        [[ "$FOUND" -eq 1 ]] && RUN_STEPS+=("$s")
    done
else
    RUN_STEPS=("${STEP_ORDER[@]}")
fi

echo "[EXECUTE_PIPELINE] Steps to run: ${RUN_STEPS[*]}"
echo "[EXECUTE_PIPELINE] Started: $(date)"

for STEP in "${RUN_STEPS[@]}"; do
    echo ""
    echo "=============================================="
    echo " [EXECUTE_PIPELINE] ${STEP} -- $(date)"
    echo "=============================================="
    if ! eval "${STEP_CMD[$STEP]}"; then
        echo "[EXECUTE_PIPELINE] FAILED at step ${STEP} -- $(date)"
        echo "[EXECUTE_PIPELINE] Fix the issue, then resume with: bash scripts/EXECUTE_PIPELINE.sh --from ${STEP}"
        exit 1
    fi
done

echo ""
echo "[EXECUTE_PIPELINE] All steps complete -- $(date)"
echo "[EXECUTE_PIPELINE] NOTE: the Figure4 tornado plot (scripts/publication_figures/tornado/) is"
echo "  not included above -- it requires a manual external ChIP realignment step first. See"
echo "  scripts/publication_figures/tornado/README.md to run it."
