#!/usr/bin/env bash
# A02_MACS3_narrow.sh -- custom MACS3 narrow peak calling + bigWig generation,
# reading the deduplicated/filtered BAMs produced by A01 (nf-core/atacseq).
# Adapted from ED26_001_ATAC_c/Narrow_allReads/scripts/02_MACS3_allReads.sh.
#
# Self-checkpointing: skips entirely if the whole-step marker already exists;
# within a run, each sample's peak-calling and bigWig generation are themselves
# skipped individually if their output already exists, so an interrupted run
# resumes without redoing finished samples.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SCREEN_SESSION="ks1_2_A02"
if [[ -z "${STY:-}" ]]; then
    echo "[A02] Not inside screen -- relaunching in screen session '${SCREEN_SESSION}' (reattach with: screen -r ${SCREEN_SESSION})."
    mkdir -p logs
    exec screen -L -Logfile "logs/A02_screen_$(date +%Y%m%d_%H%M%S).log" -S "${SCREEN_SESSION}" bash "$0" "$@"
fi

CONFIG="config/pipeline_config.yaml"
SAMPLESHEET="input/atacseq_samplesheet.csv"
BAM_DIR="data/ATAC_nfcore_output/bwa/merged_library"
OUT_DIR="data/atac_peaks"
BW_DIR="data/atac_bigwigs"
LOGDIR="logs"
DONE_MARKER="${OUT_DIR}/.pipeline_complete"
mkdir -p "$LOGDIR"

if [[ -f "$DONE_MARKER" ]]; then
    echo "[A02] Already complete (${DONE_MARKER} exists). Skipping. Delete that file to force a rerun."
    exit 0
fi

[[ -f "$CONFIG" ]] || { echo "[A02] ERROR: missing ${CONFIG}."; exit 1; }
[[ -f "$SAMPLESHEET" ]] || { echo "[A02] ERROR: missing ${SAMPLESHEET}."; exit 1; }
[[ -d "$BAM_DIR" ]] || { echo "[A02] ERROR: missing ${BAM_DIR}. Run A01_RUN_nfcore_ATAC.sh (nf-core/atacseq) to completion first."; exit 1; }

CONDA_BASE="${KS12_CONDA_BASE:-/storage/volume01/emile/miniconda3}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
# conda's own activate/deactivate machinery (and some envs' activate.d hooks, e.g.
# atac-narrow's compiler activation scripts referencing $AR/$CC/...) are not
# nounset-safe -- set -u must be off across every conda activate/deactivate call,
# not just this one, or it silently kills whichever subshell hits it (see A02 bug
# writeup in README "Methods notes").
set +u
conda activate nf-run   # provides yq, for reading config below
set -u

GSIZE=$(yq -r '.reference.effective_genome_size' "$CONFIG")
MAX_PARALLEL=$(yq -r '.resources.macs_max_parallel' "$CONFIG")
BW_THREADS=$(yq -r '.resources.bigwig_threads' "$CONFIG")
set +u
conda deactivate
set -u

mkdir -p "$OUT_DIR" "$BW_DIR"

# Sample list is derived from the samplesheet, not hardcoded -- single source of truth.
mapfile -t SAMPLES < <(tail -n +2 "$SAMPLESHEET" | cut -d',' -f1)
echo "[A02] ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

process_sample() {
    local SAMPLE="$1"
    local BAM="${BAM_DIR}/${SAMPLE}_REP1.mLb.clN.sorted.bam"
    local PEAK="${OUT_DIR}/${SAMPLE}_peaks.narrowPeak"
    local BW="${BW_DIR}/${SAMPLE}.bw"

    if [[ ! -f "$BAM" ]]; then
        echo "[A02][ERROR] BAM not found for ${SAMPLE}: ${BAM}"
        return 1
    fi

    echo "[A02] Processing ${SAMPLE} -- $(date)"
    source "${CONDA_BASE}/etc/profile.d/conda.sh"

    set +u
    conda activate atac-narrow
    set -u
    if [[ ! -f "$PEAK" ]]; then
        macs3 callpeak \
            -t "$BAM" \
            -f BAMPE \
            --nomodel \
            --nolambda \
            --keep-dup all \
            -q 0.05 \
            -g "$GSIZE" \
            -n "$SAMPLE" \
            --outdir "$OUT_DIR" \
            2>&1 | tee "${LOGDIR}/A02_macs3_${SAMPLE}.log"
        [[ -f "$PEAK" ]] || { echo "[A02][ERROR] macs3 did not produce ${PEAK} for ${SAMPLE}"; return 1; }
    else
        echo "[A02] Peaks already exist for ${SAMPLE}, skipping."
    fi
    set +u
    conda deactivate
    set -u

    set +u
    conda activate viz_env
    set -u
    if [[ ! -f "$BW" ]]; then
        bamCoverage \
            --bam "$BAM" \
            --outFileName "$BW" \
            --outFileFormat bigwig \
            --normalizeUsing RPGC \
            --effectiveGenomeSize "$GSIZE" \
            --binSize 10 \
            --centerReads \
            --extendReads \
            --ignoreDuplicates \
            --ignoreForNormalization chrM \
            --numberOfProcessors "$BW_THREADS" \
            2>&1 | tee "${LOGDIR}/A02_bamcoverage_${SAMPLE}.log"
        [[ -f "$BW" ]] || { echo "[A02][ERROR] bamCoverage did not produce ${BW} for ${SAMPLE}"; return 1; }
    else
        echo "[A02] bigWig already exists for ${SAMPLE}, skipping."
    fi
    set +u
    conda deactivate
    set -u

    echo "[A02] ${SAMPLE} complete -- $(date)"
}
export -f process_sample
export CONDA_BASE BAM_DIR OUT_DIR BW_DIR LOGDIR GSIZE BW_THREADS

declare -A PIDS
for SAMPLE in "${SAMPLES[@]}"; do
    process_sample "$SAMPLE" &
    PIDS["$SAMPLE"]=$!
    while [[ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]]; do sleep 15; done
done

FAILED=()
for SAMPLE in "${!PIDS[@]}"; do
    wait "${PIDS[$SAMPLE]}" || FAILED+=("$SAMPLE")
done

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo "[A02][ERROR] ${#FAILED[@]} sample(s) failed: ${FAILED[*]}"
    exit 1
fi

touch "$DONE_MARKER"
echo "[A02] Complete -- $(date). Peaks: ${OUT_DIR}  BigWigs: ${BW_DIR}"
