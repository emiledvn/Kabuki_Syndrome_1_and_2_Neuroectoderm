#!/usr/bin/env bash
# A06_TOBIAS.sh -- transcription-factor footprinting via TOBIAS ATACorrect ->
# FootprintScores -> BINDetect, on a CONSENSUS peak set (union of all samples'
# narrowPeaks, edge-filtered) -- not on differential-accessibility-derived peaks.
#
# Adapted from ED26_001_ATAC_c/Narrow_allReads/scripts/05_tobias_consensus.sh,
# which itself superseded two earlier, circular versions in the same directory
# (05_tobias_allReads.sh and 05_b_run_tobias_strict_combined.sh): both of those
# pre-filtered --peaks by the very differential-accessibility call (DiffBind DARs)
# that TOBIAS's own differential-binding result would then seem to "confirm" --
# testing TF binding only in regions already selected for differing chromatin
# signal between the same two conditions. The consensus approach avoids that: TOBIAS
# runs genome-wide across a neutral peak set, and any DAR/expression attribution
# happens as a separate downstream step (A07b) on the resulting differential-binding
# calls, not as a pre-filter baked into the footprinting input itself.
#
# Self-checkpointing: skips entirely if data/tobias_output/.pipeline_complete exists.
# Every sub-step (pooling, ATACorrect, FootprintScores, BINDetect) is itself
# skipped individually if its output already exists.
set -euo pipefail
renice -n 10 $$ > /dev/null 2>&1 || true
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SCREEN_SESSION="ks1_2_A06"
if [[ -z "${STY:-}" ]]; then
    echo "[A06] Not inside screen -- relaunching in screen session '${SCREEN_SESSION}' (reattach with: screen -r ${SCREEN_SESSION})."
    mkdir -p logs
    exec screen -L -Logfile "logs/A06_screen_$(date +%Y%m%d_%H%M%S).log" -S "${SCREEN_SESSION}" bash "$0" "$@"
fi

CONFIG="config/pipeline_config.yaml"
METADATA="data/sample_metadata.csv"
BAM_DIR="data/ATAC_nfcore_output/bwa/merged_library"
PEAK_DIR="data/atac_peaks"
OUT_BASE="data/tobias_output"
POOL_DIR="${OUT_BASE}/pooled_bams"
CONSENSUS="${OUT_BASE}/consensus_peaks.bed"
LOGDIR="logs"
DONE_MARKER="${OUT_BASE}/.pipeline_complete"
mkdir -p "$LOGDIR" "$OUT_BASE" "$POOL_DIR"

if [[ -f "$DONE_MARKER" ]]; then
    echo "[A06] Already complete (${DONE_MARKER} exists). Skipping. Delete that file to force a rerun."
    exit 0
fi

[[ -f "$CONFIG" ]]   || { echo "[A06] ERROR: missing ${CONFIG}."; exit 1; }
[[ -f "$METADATA" ]] || { echo "[A06] ERROR: missing ${METADATA}."; exit 1; }
[[ -d "$BAM_DIR" ]]  || { echo "[A06] ERROR: missing ${BAM_DIR}. Run A01_RUN_nfcore_ATAC.sh first."; exit 1; }
[[ -d "$PEAK_DIR" ]] || { echo "[A06] ERROR: missing ${PEAK_DIR}. Run A02_MACS3_narrow.sh first."; exit 1; }

CONDA_BASE="${KS12_CONDA_BASE:-/storage/volume01/emile/miniconda3}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate nf-run   # provides yq -- tobias_env does not have it

MOTIFS="${REPO_ROOT}/$(yq -r '.reference.motifs_meme' "$CONFIG")"
FASTA_GZ="${REPO_ROOT}/$(yq -r '.reference.fasta' "$CONFIG")"
EDGE_DIST=$(yq -r '.thresholds.min_edge_distance' "$CONFIG")
CORES=$(yq -r '.resources.tobias_cores' "$CONFIG")

# Read all per-contrast values now too, while yq is available -- used later
# once tobias_env is active instead.
N_CONTRASTS=$(yq -r '.contrasts | length' "$CONFIG")
CONTRAST_NAMES=(); CONTRAST_TREATMENTS=(); CONTRAST_REFERENCES=()
for ((i=0; i<N_CONTRASTS; i++)); do
    CONTRAST_NAMES+=("$(yq -r ".contrasts[$i].name" "$CONFIG")")
    CONTRAST_TREATMENTS+=("$(yq -r ".contrasts[$i].treatment" "$CONFIG")")
    CONTRAST_REFERENCES+=("$(yq -r ".contrasts[$i].reference" "$CONFIG")")
done
conda deactivate

conda activate tobias_env
for tool in TOBIAS samtools bedtools; do
    command -v "$tool" &>/dev/null || { echo "[A06] ERROR: $tool not found in tobias_env"; exit 1; }
done
[[ -f "$MOTIFS" ]] || { echo "[A06] ERROR: missing ${MOTIFS}."; exit 1; }

# TOBIAS needs an uncompressed, indexed FASTA. We keep one shared copy in
# data/reference/ derived from our own pinned reference (not nf-core's internal
# per-run copy) so "one genome reference used throughout" stays literally true.
GENOME="${FASTA_GZ%.gz}"
FAI="${GENOME}.fai"
if [[ ! -f "$GENOME" ]]; then
    echo "[A06] Unzipping shared reference FASTA (one-time, cached for future steps)..."
    gunzip -k -c "$FASTA_GZ" > "${GENOME}.tmp"
    mv "${GENOME}.tmp" "$GENOME"
fi
[[ -f "$FAI" ]] || samtools faidx "$GENOME"

########
## 1. Consensus peak set = union of all narrowPeaks, edge-filtered
########
if [[ ! -s "$CONSENSUS" ]]; then
    echo "[A06] Building consensus peak set from ${PEAK_DIR}/*_peaks.narrowPeak"
    NP=$(ls "${PEAK_DIR}"/*_peaks.narrowPeak 2>/dev/null | wc -l)
    [[ "$NP" -eq 0 ]] && { echo "[A06] ERROR: no narrowPeak files in ${PEAK_DIR}"; exit 1; }

    cat "${PEAK_DIR}"/*_peaks.narrowPeak \
        | cut -f1-3 \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - > "${OUT_BASE}/consensus_raw.bed"

    awk -v edge="$EDGE_DIST" 'BEGIN{OFS="\t"} NR==FNR{len[$1]=$2; next}
        ($1 in len) && ($2>=edge) && ($3<=(len[$1]-edge)) {print}' \
        "$FAI" "${OUT_BASE}/consensus_raw.bed" > "$CONSENSUS"
    rm -f "${OUT_BASE}/consensus_raw.bed"
    echo "[A06] Consensus peaks (edge-filtered): $(wc -l < "$CONSENSUS")"
else
    echo "[A06] Consensus exists: $(wc -l < "$CONSENSUS") peaks"
fi

########
## 2. Pool BAMs per genotype (genotype list derived from sample_metadata.csv)
########
mapfile -t GENOTYPES < <(tail -n +2 "$METADATA" | awk -F',' '$3=="ATAC"{print $4}' | sort -u)
echo "[A06] Genotypes: ${GENOTYPES[*]}"

pool_condition() {
    local GENO="$1"
    local POOLED="${POOL_DIR}/${GENO}_pooled.bam"
    if [[ -f "${POOLED}.bai" ]]; then
        echo "[A06] Pooled BAM exists: $POOLED"
        return 0
    fi
    mapfile -t REPS < <(tail -n +2 "$METADATA" | awk -F',' -v g="$GENO" '$3=="ATAC" && $4==g{print $1}')
    local BAMS=()
    for s in "${REPS[@]}"; do BAMS+=("${BAM_DIR}/${s}_REP1.mLb.clN.sorted.bam"); done
    for b in "${BAMS[@]}"; do [[ -f "$b" ]] || { echo "[A06] ERROR: missing $b"; exit 1; }; done
    echo "[A06] Pooling ${GENO}: ${#BAMS[@]} BAM(s)"
    samtools merge -f -@ "$CORES" "$POOLED" "${BAMS[@]}"
    samtools index "$POOLED" -@ "$CORES"
}
for GENO in "${GENOTYPES[@]}"; do pool_condition "$GENO"; done

########
## 3. TOBIAS per contrast (config-driven, shared with A04/R02/R04), on the same
##    consensus peak set
########
run_tobias() {
    local LABEL="$1" REF_BAM="$2" TREAT_BAM="$3" TREAT_NAME="$4" REF_NAME="$5" OUT_DIR="$6"
    mkdir -p "${OUT_DIR}/ATACorrect" "${OUT_DIR}/Footprints" "${OUT_DIR}/BINDetect"

    local REF_CORR="${OUT_DIR}/ATACorrect/${REF_NAME}_${LABEL}_corrected.bw"
    local TR_CORR="${OUT_DIR}/ATACorrect/${LABEL}_corrected.bw"
    local REF_FP="${OUT_DIR}/Footprints/${REF_NAME}_${LABEL}_footprints.bw"
    local TR_FP="${OUT_DIR}/Footprints/${LABEL}_footprints.bw"

    [[ ! -f "$REF_CORR" ]] && TOBIAS ATACorrect --bam "$REF_BAM" --genome "$GENOME" \
        --peaks "$CONSENSUS" --outdir "${OUT_DIR}/ATACorrect" \
        --prefix "${REF_NAME}_${LABEL}" --cores "$CORES"
    [[ ! -f "$TR_CORR" ]] && TOBIAS ATACorrect --bam "$TREAT_BAM" --genome "$GENOME" \
        --peaks "$CONSENSUS" --outdir "${OUT_DIR}/ATACorrect" \
        --prefix "${LABEL}" --cores "$CORES"

    [[ ! -f "$REF_FP" ]] && TOBIAS FootprintScores --signal "$REF_CORR" \
        --regions "$CONSENSUS" --output "$REF_FP" --cores "$CORES"
    [[ ! -f "$TR_FP" ]] && TOBIAS FootprintScores --signal "$TR_CORR" \
        --regions "$CONSENSUS" --output "$TR_FP" --cores "$CORES"

    echo "[A06] BINDetect: ${REF_NAME} vs ${TREAT_NAME} -- consensus peaks, full motif set"
    TOBIAS BINDetect \
        --motifs "$MOTIFS" \
        --signals "$REF_FP" "$TR_FP" \
        --cond_names "$REF_NAME" "$TREAT_NAME" \
        --genome "$GENOME" \
        --peaks "$CONSENSUS" \
        --outdir "${OUT_DIR}/BINDetect" \
        --cores "$CORES"
}

for ((i=0; i<N_CONTRASTS; i++)); do
    NAME="${CONTRAST_NAMES[$i]}"
    TREATMENT="${CONTRAST_TREATMENTS[$i]}"
    REFERENCE="${CONTRAST_REFERENCES[$i]}"
    echo ""
    echo "=============================================="
    echo " ${NAME}: ${REFERENCE} vs ${TREATMENT}"
    echo "=============================================="
    run_tobias \
        "$TREATMENT" \
        "${POOL_DIR}/${REFERENCE}_pooled.bam" \
        "${POOL_DIR}/${TREATMENT}_pooled.bam" \
        "$TREATMENT" \
        "$REFERENCE" \
        "${OUT_BASE}/${NAME}"
done

touch "$DONE_MARKER"
echo ""
echo "[A06] Complete -- $(date). Results: ${OUT_BASE}/"
echo "[A06] NOTE: DAR/expression attribution happens in A07b_TOBIAS_DOWNSTREAM (not yet written)."
conda deactivate
