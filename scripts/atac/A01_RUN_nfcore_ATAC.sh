#!/usr/bin/env bash
# A01_RUN_nfcore_ATAC.sh -- runs nf-core/atacseq via Nextflow + Apptainer.
#
# Self-checkpointing: skips if data/ATAC_nfcore_output/.pipeline_complete already
# exists (delete it to force a rerun). Requires input/atacseq_samplesheet.csv and
# the shared genome reference (config/pipeline_config.yaml) to exist first.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SCREEN_SESSION="ks1_2_A01"
if [[ -z "${STY:-}" ]]; then
    echo "[A01] Not inside screen -- relaunching in screen session '${SCREEN_SESSION}' (reattach with: screen -r ${SCREEN_SESSION})."
    mkdir -p logs
    exec screen -L -Logfile "logs/A01_screen_$(date +%Y%m%d_%H%M%S).log" -S "${SCREEN_SESSION}" bash "$0" "$@"
fi

CONFIG="config/pipeline_config.yaml"
SAMPLESHEET="input/atacseq_samplesheet.csv"
OUTDIR="data/ATAC_nfcore_output"
WORKDIR="data/work_atac"
DONE_MARKER="${OUTDIR}/.pipeline_complete"
LOGDIR="logs"
mkdir -p "$LOGDIR"

if [[ -f "$DONE_MARKER" ]]; then
    echo "[A01] Already complete (${DONE_MARKER} exists). Skipping. Delete that file to force a rerun."
    exit 0
fi

[[ -f "$CONFIG" ]] || { echo "[A01] ERROR: missing ${CONFIG}."; exit 1; }
[[ -f "$SAMPLESHEET" ]] || { echo "[A01] ERROR: missing ${SAMPLESHEET}. Generate it with: python3 data/data_import_local/make_samplesheets.py"; exit 1; }

CONDA_BASE="${KS12_CONDA_BASE:-/storage/volume01/emile/miniconda3}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate nf-run

FASTA="${REPO_ROOT}/$(yq -r '.reference.fasta' "$CONFIG")"
GTF="${REPO_ROOT}/$(yq -r '.reference.gtf' "$CONFIG")"
GSIZE=$(yq -r '.reference.effective_genome_size' "$CONFIG")
[[ -f "$FASTA" ]] || { echo "[A01] ERROR: missing genome FASTA at ${FASTA}. See ${CONFIG}:reference.fasta_url."; exit 1; }
[[ -f "$GTF" ]]   || { echo "[A01] ERROR: missing GTF at ${GTF}. See ${CONFIG}:reference.gtf_url."; exit 1; }

MAX_CPUS=$(yq -r '.resources.max_cpus' "$CONFIG")
MAX_MEMORY=$(yq -r '.resources.max_memory' "$CONFIG")
MAX_TIME=$(yq -r '.resources.max_time' "$CONFIG")
export KS12_QUEUE_SIZE=$(yq -r '.resources.queue_size' "$CONFIG")
export KS12_ALIGN_MAX_FORKS=$(yq -r '.resources.align_max_forks' "$CONFIG")
export KS12_ALIGN_CPUS=$(yq -r '.resources.align_cpus' "$CONFIG")
REVISION=$(yq -r '.nfcore.atacseq_revision' "$CONFIG")

mkdir -p "$OUTDIR" "$WORKDIR"

TOWER_ARGS=()
if [[ -n "${TOWER_ACCESS_TOKEN:-}" ]]; then
    TOWER_ARGS=(-with-tower)
    echo "[A01] TOWER_ACCESS_TOKEN set -- reporting to Seqera Tower."
else
    echo "[A01] TOWER_ACCESS_TOKEN not set -- running without Tower monitoring."
fi

LOG="${LOGDIR}/A01_nextflow_$(date +%Y%m%d_%H%M%S).log"
echo "[A01] Starting nf-core/atacseq ${REVISION} -- $(date)"
echo "[A01] fasta=${FASTA} gtf=${GTF} outdir=${OUTDIR}"
echo "[A01] Log: ${LOG}"

set +e
nextflow run nf-core/atacseq \
    -revision "$REVISION" \
    -profile apptainer \
    -c config/nextflow_resources.config \
    -w "$WORKDIR" \
    -resume \
    "${TOWER_ARGS[@]}" \
    --input        "$SAMPLESHEET" \
    --outdir       "$OUTDIR" \
    --fasta        "$FASTA" \
    --gtf          "$GTF" \
    --macs_gsize   "$GSIZE" \
    --max_cpus     "$MAX_CPUS" \
    --max_memory   "$MAX_MEMORY" \
    --max_time     "$MAX_TIME" \
    2>&1 | tee "$LOG"
STATUS=${PIPESTATUS[0]}
set -e

if [[ $STATUS -ne 0 ]]; then
    echo "[A01] FAILED (exit ${STATUS}) -- $(date). Check ${LOG}. Rerun this script to resume (-resume uses Nextflow's own cache)."
    exit "$STATUS"
fi

touch "$DONE_MARKER"
echo "[A01] Complete -- $(date). Marker: ${DONE_MARKER}"
