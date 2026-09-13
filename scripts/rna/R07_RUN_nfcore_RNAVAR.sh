#!/usr/bin/env bash
# R07_RUN_nfcore_RNAVAR.sh -- runs nf-core/rnavar via Nextflow + Apptainer to call
# per-sample RNA-seq variants (GATK4 RNA-seq short-variant best practices: 2-pass
# STAR, MarkDuplicates, SplitNCigarReads, BQSR, HaplotypeCaller, hard-filtering).
# Feeds R08_eSNP_Karyotype.R's BAF-by-chromosome-arm aneuploidy screen -- see that
# script's header and docs/ for the eSNP-Karyotyping method (Weissbein et al. 2016).
#
# Structural clone of R01_RUN_nfcore_RNA.sh -- same samplesheet, same reference,
# same screen/self-checkpoint/Tower conventions -- plus a one-time download of the
# GATK BQSR known-sites bundle (dbSNP + Mills/1000G known indels) this pipeline
# didn't previously need.
#
# Self-checkpointing: skips if data/RNAvar_nfcore_output/.pipeline_complete already
# exists (delete it to force a rerun). Requires input/rnaseq_samplesheet.csv and
# the shared genome reference (config/pipeline_config.yaml) to exist first.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SCREEN_SESSION="ks1_2_R07"
if [[ -z "${STY:-}" ]]; then
    echo "[R07] Not inside screen -- relaunching in screen session '${SCREEN_SESSION}' (reattach with: screen -r ${SCREEN_SESSION})."
    mkdir -p logs
    exec screen -L -Logfile "logs/R07_screen_$(date +%Y%m%d_%H%M%S).log" -S "${SCREEN_SESSION}" bash "$0" "$@"
fi

CONFIG="config/pipeline_config.yaml"
SAMPLESHEET="input/rnaseq_samplesheet.csv"
OUTDIR="data/RNAvar_nfcore_output"
DONE_MARKER="${OUTDIR}/.pipeline_complete"
LOGDIR="logs"
mkdir -p "$LOGDIR"

if [[ -f "$DONE_MARKER" ]]; then
    echo "[R07] Already complete (${DONE_MARKER} exists). Skipping. Delete that file to force a rerun."
    exit 0
fi

[[ -f "$CONFIG" ]] || { echo "[R07] ERROR: missing ${CONFIG}."; exit 1; }
[[ -f "$SAMPLESHEET" ]] || { echo "[R07] ERROR: missing ${SAMPLESHEET}. Generate it with: python3 data/data_import_local/make_samplesheets.py"; exit 1; }

CONDA_BASE="${KS12_CONDA_BASE:-/storage/volume01/emile/miniconda3}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate nf-run

FASTA="${REPO_ROOT}/$(yq -r '.reference.fasta' "$CONFIG")"
GTF="${REPO_ROOT}/$(yq -r '.reference.gtf' "$CONFIG")"
[[ -f "$FASTA" ]] || { echo "[R07] ERROR: missing genome FASTA at ${FASTA}. See ${CONFIG}:reference.fasta_url."; exit 1; }
[[ -f "$GTF" ]]   || { echo "[R07] ERROR: missing GTF at ${GTF}. See ${CONFIG}:reference.gtf_url."; exit 1; }

# nf-core/rnavar's param schema requires an uncompressed --gtf (regex ^\S+\.gtf$ --
# unlike nf-core/rnaseq, which R01 already feeds the .gz directly). Use the
# already-on-disk uncompressed FASTA (data/reference already keeps both forms for
# other reasons) and derive a one-time uncompressed GTF copy alongside the pinned
# .gz (kept, not replaced -- every other script in this pipeline still reads the
# .gz via the shared reference.gtf config key).
FASTA_PLAIN="${FASTA%.gz}"
[[ -f "$FASTA_PLAIN" ]] || { echo "[R07] ERROR: expected uncompressed FASTA at ${FASTA_PLAIN} (nf-core/rnaseq's R01 run should have left it)."; exit 1; }
GTF_PLAIN="${GTF%.gz}"
if [[ ! -f "$GTF_PLAIN" ]]; then
    echo "[R07] Deriving uncompressed GTF (one-time, for rnavar's schema only): ${GTF_PLAIN}"
    gunzip -k -c "$GTF" > "$GTF_PLAIN"
fi

# R07-specific resource overrides (rnavar_*), deliberately separate from the shared
# resources.max_cpus/max_memory/queue_size/align_* R01/A01/etc. use -- see
# pipeline_config.yaml's comment there for why (CUVERTINO_Reanalysis precedent).
MAX_CPUS=$(yq -r '.resources.rnavar_max_cpus' "$CONFIG")
MAX_MEMORY=$(yq -r '.resources.rnavar_max_memory' "$CONFIG")
MAX_TIME=$(yq -r '.resources.max_time' "$CONFIG")
export KS12_QUEUE_SIZE=$(yq -r '.resources.rnavar_queue_size' "$CONFIG")
export KS12_ALIGN_MAX_FORKS=$(yq -r '.resources.rnavar_align_max_forks' "$CONFIG")
export KS12_ALIGN_CPUS=$(yq -r '.resources.rnavar_align_cpus' "$CONFIG")
REVISION=$(yq -r '.nfcore.rnavar_revision' "$CONFIG")
RAM_WORKDIR_MIN_FREE_GB=$(yq -r '.resources.rnavar_ram_workdir_min_free_gb' "$CONFIG")

# Nextflow's work dir is where nearly all of a run's actual I/O happens (staged
# inputs, sorted/indexed BAMs, per-task GATK scratch) -- putting it on /dev/shm
# (tmpfs, RAM-backed) instead of NFS-backed data/work_rnavar removes that traffic
# from the network mount entirely. Same trick + same trade-off CUVERTINO_Reanalysis's
# 03_run_nfcore_rna.sh uses: tmpfs is volatile (a reboot mid-run loses -resume-ability
# -- a full restart, not a correctness risk), acceptable given the free space checked
# below. Final --outdir output still lands on NFS as normal, only the disposable
# intermediate cache moves to RAM. Falls back to NFS-backed data/work_rnavar if
# /dev/shm is missing or has less than rnavar_ram_workdir_min_free_gb free.
if [[ -d /dev/shm ]] && [[ $(df --output=avail /dev/shm | tail -1) -gt $((RAM_WORKDIR_MIN_FREE_GB * 1024 * 1024)) ]]; then
    WORKDIR="/dev/shm/${USER}_ks1_2_work_rnavar"
    echo "[R07] Using RAM-backed work dir: ${WORKDIR} (clean up manually with 'rm -rf ${WORKDIR}' once results are confirmed good -- persists in RAM until removed or reboot)."
else
    WORKDIR="data/work_rnavar"
    echo "[R07] /dev/shm not usable (missing or <${RAM_WORKDIR_MIN_FREE_GB}G free) -- falling back to NFS-backed work dir: ${WORKDIR}"
fi

mkdir -p "$OUTDIR" "$WORKDIR"

########
## One-time GATK BQSR known-sites bundle download (dbSNP + Mills/1000G indels).
## Not part of the shared reference.* download pattern (fasta/gtf) because these
## are ~1.6GB combined and used only by this one arm -- kept out of the routine
## reference set so a fresh checkout of this repo doesn't imply pulling them.
## curl -C - allows resuming a partial download across reruns (these files are
## large enough that a mid-download interruption is realistic).
########

# nf-core/rnavar 1.3.0 requires Nextflow >=25.10.4; the shared nf-run conda env is
# pinned to 25.10.2 (used live by other analysis/ED2x_* projects too -- not ours to
# bump in place). A project-scoped binary (gitignored under data/*, doesn't touch
# shared infra) is installed on demand instead and preferred over conda's if present.
NEXTFLOW_BIN="nextflow"
if [[ -x "${REPO_ROOT}/data/tools/nextflow" ]]; then
    NEXTFLOW_BIN="${REPO_ROOT}/data/tools/nextflow"
    echo "[R07] Using project-local Nextflow: $("$NEXTFLOW_BIN" -version 2>/dev/null | grep version)"
fi

DBSNP="${REPO_ROOT}/$(yq -r '.reference.gatk_dbsnp_vcf' "$CONFIG")"
DBSNP_URL=$(yq -r '.reference.gatk_dbsnp_vcf_url' "$CONFIG")
KNOWN_INDELS="${REPO_ROOT}/$(yq -r '.reference.gatk_known_indels_vcf' "$CONFIG")"
KNOWN_INDELS_URL=$(yq -r '.reference.gatk_known_indels_vcf_url' "$CONFIG")
GATK_BUNDLE_DIR="$(dirname "$DBSNP")"
mkdir -p "$GATK_BUNDLE_DIR"

download_with_index() {
    local url="$1" dest="$2" label="$3"
    if [[ ! -f "$dest" ]]; then
        echo "[R07] Downloading ${label} -- $(date)"
        curl -L --retry 3 --retry-delay 5 -C - -o "$dest" "$url"
        echo "[R07] ${label} md5: $(md5sum "$dest" | awk '{print $1}')  (not pre-pinned -- see ${CONFIG}:reference comments; add it there once confirmed stable)"
    else
        echo "[R07] ${label} already present: ${dest}"
    fi
    if [[ ! -f "${dest}.tbi" ]]; then
        echo "[R07] Downloading ${label} index -- $(date)"
        curl -L --retry 3 --retry-delay 5 -C - -o "${dest}.tbi" "${url}.tbi"
    fi
}

download_with_index "$DBSNP_URL" "$DBSNP" "GATK dbSNP known-sites"
download_with_index "$KNOWN_INDELS_URL" "$KNOWN_INDELS" "GATK Mills/1000G known-indels"

TOWER_ARGS=()
if [[ -n "${TOWER_ACCESS_TOKEN:-}" ]]; then
    TOWER_ARGS=(-with-tower)
    echo "[R07] TOWER_ACCESS_TOKEN set -- reporting to Seqera Tower."
else
    echo "[R07] TOWER_ACCESS_TOKEN not set -- running without Tower monitoring."
fi

# --dbsnp/--known_indels enable GATK's BaseRecalibrator/ApplyBQSR (BQSR) steps --
# deliberately NOT --skip_baserecalibration, per the full-best-practices choice for
# this arm. --aligner/--generate_gvcf/--remove_duplicates left at pipeline defaults:
# we want per-sample filtered VCFs (not gVCFs/joint genotyping) for R08's per-sample
# BAF scan.
LOG="${LOGDIR}/R07_nextflow_$(date +%Y%m%d_%H%M%S).log"
echo "[R07] Starting nf-core/rnavar ${REVISION} -- $(date)"
echo "[R07] fasta=${FASTA_PLAIN} gtf=${GTF_PLAIN} outdir=${OUTDIR}"
echo "[R07] dbsnp=${DBSNP} known_indels=${KNOWN_INDELS}"
echo "[R07] Log: ${LOG}"

set +e
"$NEXTFLOW_BIN" run nf-core/rnavar \
    -revision "$REVISION" \
    -profile apptainer \
    -c config/nextflow_resources.config \
    -w "$WORKDIR" \
    -resume \
    "${TOWER_ARGS[@]}" \
    --input          "$SAMPLESHEET" \
    --outdir         "$OUTDIR" \
    --fasta          "$FASTA_PLAIN" \
    --gtf            "$GTF_PLAIN" \
    --dbsnp          "$DBSNP" \
    --dbsnp_tbi      "${DBSNP}.tbi" \
    --known_indels     "$KNOWN_INDELS" \
    --known_indels_tbi "${KNOWN_INDELS}.tbi" \
    --max_cpus       "$MAX_CPUS" \
    --max_memory     "$MAX_MEMORY" \
    --max_time       "$MAX_TIME" \
    2>&1 | tee "$LOG"
STATUS=${PIPESTATUS[0]}
set -e

if [[ $STATUS -ne 0 ]]; then
    echo "[R07] FAILED (exit ${STATUS}) -- $(date). Check ${LOG}. Rerun this script to resume (-resume uses Nextflow's own cache)."
    exit "$STATUS"
fi

touch "$DONE_MARKER"
echo "[R07] Complete -- $(date). Marker: ${DONE_MARKER}"
