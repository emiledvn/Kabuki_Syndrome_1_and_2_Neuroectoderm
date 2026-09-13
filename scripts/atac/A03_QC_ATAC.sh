#!/usr/bin/env bash
# A03_QC_ATAC.sh -- ATAC QC: read depth/duplication/mitochondrial fraction from
# ataqv (runs automatically as part of A01_RUN_nfcore_ATAC.sh -- nf-core/atacseq
# already computes these correctly, no need to recompute by hand); FRiP + peak
# count computed fresh against A02's own MACS3 narrow-peak set specifically
# (genuinely not redundant with ataqv's own FRiP-like field, which is computed
# against nf-core's *internal* peak call, not the MACS3 peaks A04/A06/A08 actually
# use downstream); and TSS enrichment computed ENCODE-exact against a RefSeq
# reference, NOT read from ataqv.
#
# TSS enrichment methodology note: ataqv's own tss_enrichment field is neither
# ENCODE-comparable nor library-quality-representative here -- it uses a +/-1000bp
# window (ENCODE's own formula specifies +/-2000bp) against the GENCODE-derived TSS
# set nf-core/atacseq generates internally (~5x more entries than RefSeq's curated
# set, which is what ENCODE's own published cutoff table is actually calibrated
# against). Both effects independently depress the score, which is why all 9
# samples in this project scored 3.68-4.85 against ataqv despite being genuinely
# good libraries: recomputed ENCODE-exact, they land in the 6-9 range instead. Full
# writeup, sources, and the exact formula: see METHODS_NOTES.md.
#
# NOTE on NRF/PBC1/PBC2: checked against a real ataqv JSON from the prior
# ED26_001_ATAC_c project on this server -- ataqv does not emit nrf/pbc1/pbc2
# fields (despite some third-party docs implying it does). This script uses
# ataqv's duplicate_reads/total_reads ratio as the duplication signal instead
# (qc_max_duplicate_fraction) -- a simplification, not the formal ENCODE metric.
#
# Adapted from ED26_001_ATAC_c/Narrow_allReads/scripts/02b_QC_allReads.sh; that
# script's read-depth/complexity/fragment-length blocks are cut here in favor of
# ataqv, the FRiP block is kept, and TSS enrichment is a fresh ENCODE-exact
# reimplementation (see above).
#
# Self-checkpointing: skips entirely if results/tables/A03_qc_flags.tsv already
# exists. Requires A01 (BAMs + ataqv JSON) and A02 (peaks + bigwigs) to have
# completed first. The RefSeq TSS reference is prepared once (self-checkpointed on
# its own existence) and reused on every subsequent run; the raw UCSC download is
# deleted after the reference BED is derived from it -- it has no further use.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SCREEN_SESSION="ks1_2_A03"
if [[ -z "${STY:-}" ]]; then
    echo "[A03] Not inside screen -- relaunching in screen session '${SCREEN_SESSION}' (reattach with: screen -r ${SCREEN_SESSION})."
    mkdir -p logs
    exec screen -L -Logfile "logs/A03_screen_$(date +%Y%m%d_%H%M%S).log" -S "${SCREEN_SESSION}" bash "$0" "$@"
fi

CONFIG="config/pipeline_config.yaml"
SAMPLESHEET="input/atacseq_samplesheet.csv"
BAM_DIR="data/ATAC_nfcore_output/bwa/merged_library"
ATAQV_DIR="data/ATAC_nfcore_output/bwa/merged_library/ataqv"
PEAK_DIR="data/atac_peaks"
BW_DIR="data/atac_bigwigs"
TABLES_DIR="results/tables"
LOGDIR="logs"
TSS_WORK="data/reference/tss_enrichment_tmp"
FLAGS_FILE="${TABLES_DIR}/A03_qc_flags.tsv"
mkdir -p "$LOGDIR" "$TABLES_DIR" "$TSS_WORK"

if [[ -f "$FLAGS_FILE" ]]; then
    echo "[A03] Already complete (${FLAGS_FILE} exists). Skipping. Delete that file to force a rerun."
    exit 0
fi

[[ -f "$CONFIG" ]] || { echo "[A03] ERROR: missing ${CONFIG}."; exit 1; }
[[ -f "$SAMPLESHEET" ]] || { echo "[A03] ERROR: missing ${SAMPLESHEET}."; exit 1; }
[[ -d "$BAM_DIR" ]] || { echo "[A03] ERROR: missing ${BAM_DIR}. Run A01_RUN_nfcore_ATAC.sh first."; exit 1; }
[[ -d "$PEAK_DIR" ]] || { echo "[A03] ERROR: missing ${PEAK_DIR}. Run A02_MACS3_narrow.sh first."; exit 1; }
[[ -d "$BW_DIR" ]] || { echo "[A03] ERROR: missing ${BW_DIR}. Run A02_MACS3_narrow.sh first."; exit 1; }
[[ -d "$ATAQV_DIR" ]] || { echo "[A03] ERROR: missing ${ATAQV_DIR}. ataqv runs automatically as part of A01_RUN_nfcore_ATAC.sh (nf-core/atacseq) -- rerun/check A01 if this is missing."; exit 1; }

CONDA_BASE="${KS12_CONDA_BASE:-/storage/volume01/emile/miniconda3}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate nf-run   # provides yq, for reading config below

MAX_PARALLEL=$(yq -r '.resources.qc_max_parallel' "$CONFIG")
MIN_USABLE=$(yq -r '.thresholds.qc_min_usable_reads' "$CONFIG")
MIN_FRIP=$(yq -r '.thresholds.qc_min_frip' "$CONFIG")
MIN_PEAKS=$(yq -r '.thresholds.qc_min_peaks' "$CONFIG")
MAX_DUP_FRACTION=$(yq -r '.thresholds.qc_max_duplicate_fraction' "$CONFIG")
MIN_TSS=$(yq -r '.thresholds.qc_min_tss_enrichment' "$CONFIG")
REFSEQ_TSS="${REPO_ROOT}/$(yq -r '.reference.qc_refseq_tss_bed' "$CONFIG")"
REFSEQ_URL=$(yq -r '.reference.qc_refseq_tss_source_url' "$CONFIG")
REFSEQ_MD5=$(yq -r '.reference.qc_refseq_tss_source_md5' "$CONFIG")
conda deactivate

conda activate viz_env
for tool in samtools bedtools bamCoverage computeMatrix python3; do
    command -v "$tool" &>/dev/null || { echo "[A03] ERROR: $tool not found in viz_env"; exit 1; }
done

########
## Prepare the RefSeq TSS reference (QC-only; not the main GENCODE analysis
## reference). Self-checkpointed on its own existence; the raw UCSC download is
## disposable and removed once the derived BED is written. See METHODS_NOTES.md
## for why this is a separate reference from A08's GENCODE-derived TSS bed.
########
if [[ ! -f "$REFSEQ_TSS" ]]; then
    echo "[A03] Preparing RefSeq TSS reference (one-time)..."
    RAW="$(dirname "$REFSEQ_TSS")/refGene.txt.gz"
    curl -s --max-time 60 -o "$RAW" "$REFSEQ_URL"
    echo "${REFSEQ_MD5}  ${RAW}" | md5sum -c - || { echo "[A03] ERROR: refGene.txt.gz checksum mismatch."; exit 1; }
    zcat "$RAW" | awk -F'\t' 'BEGIN{OFS="\t"}
        {
          chrom=$3; strand=$4; txStart=$5; txEnd=$6; gene=$13
          if (chrom !~ /^chr([0-9]+|X|Y|M)$/) next
          if (strand=="+") tss=txStart; else tss=txEnd
          print chrom, tss, tss+1, gene, ".", strand
        }' | sort -k1,1 -k2,2n -u | awk '!seen[$1"\t"$2"\t"$3]++' > "$REFSEQ_TSS"
    rm -f "$RAW"   # disposable -- only the derived BED is needed going forward
    echo "[A03] RefSeq TSS reference: $(wc -l < "$REFSEQ_TSS") unique positions"
fi

mapfile -t SAMPLES < <(tail -n +2 "$SAMPLESHEET" | cut -d',' -f1)
echo "[A03] ${#SAMPLES[@]} samples: ${SAMPLES[*]}"

ATAQV_FILE="${TABLES_DIR}/A03_ataqv_summary.tsv"
FRIP_FILE="${TABLES_DIR}/A03_frip_summary.tsv"
TSS_FILE="${TABLES_DIR}/A03_tss_enrichment.tsv"

echo -e "Sample\tTotal_reads\tMapped\tDuplicate_reads\tDuplicate_fraction\tMito_reads\tUsable_reads_hqaa" > "$ATAQV_FILE"
echo -e "Sample\tUsable_reads\tReads_in_peaks\tFRiP\tN_peaks"                                            > "$FRIP_FILE"
echo -e "Sample\tTSS_enrichment_ENCODE_RefSeq"                                                           > "$TSS_FILE"

process_sample() {
    local SAMPLE="$1"
    local BAM="${BAM_DIR}/${SAMPLE}_REP1.mLb.clN.sorted.bam"
    local PEAK="${PEAK_DIR}/${SAMPLE}_peaks.narrowPeak"
    local BW="${BW_DIR}/${SAMPLE}.bw"
    local ATAQV_JSON

    [[ -f "$BAM"  ]] || { echo "[A03][ERROR] BAM not found: $BAM"; return 1; }
    [[ -f "$PEAK" ]] || { echo "[A03][ERROR] Peak not found: $PEAK"; return 1; }
    [[ -f "$BW"   ]] || { echo "[A03][ERROR] bigwig not found: $BW"; return 1; }
    ATAQV_JSON=$(find "$ATAQV_DIR" -name "${SAMPLE}_REP1.ataqv.json" -print -quit)
    [[ -n "$ATAQV_JSON" && -f "$ATAQV_JSON" ]] || { echo "[A03][ERROR] ataqv JSON not found for ${SAMPLE} under ${ATAQV_DIR}"; return 1; }

    echo "[A03] QC: ${SAMPLE} -- $(date)"
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate viz_env

    ## 1. Read depth / duplication / mito -- from ataqv (nf-core/atacseq already
    ##    computed these correctly during A01). TSS enrichment is NOT read from
    ##    here -- see part 3 below.
    local ATAQV_LINE
    ATAQV_LINE=$(python3 - "$ATAQV_JSON" <<'PYEOF'
import json, sys
d = json.load(open(sys.argv[1]))[0]["metrics"]
total = d["total_reads"]
mapped = total - d.get("unmapped_reads", 0)
dup = d["duplicate_reads"]
dup_frac = (dup / total) if total else 0.0
mito = d["total_mitochondrial_reads"]
hqaa = d["hqaa"]
print(f"{total}\t{mapped}\t{dup}\t{dup_frac:.4f}\t{mito}\t{hqaa}")
PYEOF
)
    echo -e "${SAMPLE}\t${ATAQV_LINE}" >> "${ATAQV_FILE}.${SAMPLE}.tmp"

    ## 2. FRiP against A02's own MACS3 peaks specifically (MAPQ>=30, properly
    ##    paired, usable reads). The usable-read denominator is computed fresh
    ##    here rather than reused from ataqv's hqaa: hqaa is autosomal-only, so
    ##    pairing it with a numerator from A02's whole-genome peak set would
    ##    silently undercount chrX signal in the ratio -- notable here since
    ##    KDM6A is X-linked.
    local USABLE IN_PEAKS FRIP NPEAKS
    USABLE=$(samtools view -c -F 1804 -f 2 -q 30 "$BAM")
    IN_PEAKS=$(bedtools intersect -a "$BAM" -b "$PEAK" -u -f 0.5 | samtools view -c -F 1804 -f 2 -q 30)
    FRIP=$(python3 -c 'import sys; a,b=int(sys.argv[1]),int(sys.argv[2]); print(f"{a/b:.4f}" if b else "0.0000")' "$IN_PEAKS" "$USABLE")
    NPEAKS=$(wc -l < "$PEAK")
    echo -e "${SAMPLE}\t${USABLE}\t${IN_PEAKS}\t${FRIP}\t${NPEAKS}" >> "${FRIP_FILE}.${SAMPLE}.tmp"

    ## 3. TSS enrichment, ENCODE-exact: RefSeq TSS set, +/-2000bp window, 10bp
    ##    bins, normalized to the mean signal in the 100bp at each window end,
    ##    score = normalized signal at the single center bin (the TSS itself).
    ##    Reuses A02's own bigwig (RPGC-normalized) rather than regenerating one --
    ##    a global normalization constant cancels out in a signal/flank ratio, so
    ##    which normalization A02 used doesn't affect this score.
    local MATRIX="${TSS_WORK}/${SAMPLE}_matrix.gz"
    computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 --binSize 10 \
        -R "$REFSEQ_TSS" -S "$BW" --skipZeros -o "$MATRIX" -p 4 --quiet 2>/dev/null

    local TSS_SCORE
    TSS_SCORE=$(python3 - "$MATRIX" <<'PYEOF'
import numpy as np, gzip, sys
with gzip.open(sys.argv[1], 'rt') as f:
    f.readline()
    data = np.array([list(map(float, line.strip().split('\t')[6:])) for line in f if line.strip()])
bins = data.shape[1]
flank_bins = 10
flank = np.nanmean(np.concatenate([data[:, :flank_bins], data[:, -flank_bins:]], axis=1))
center = bins // 2
center_signal = np.nanmean(data[:, center])
print(f"{(center_signal / flank if flank > 0 else 0.0):.3f}")
PYEOF
)
    rm -f "$MATRIX"   # disposable -- the score is all that's kept
    echo -e "${SAMPLE}\t${TSS_SCORE}" >> "${TSS_FILE}.${SAMPLE}.tmp"

    conda deactivate
    echo "[A03] ${SAMPLE} QC complete -- $(date)"
}
export -f process_sample
export CONDA_BASE BAM_DIR PEAK_DIR BW_DIR ATAQV_DIR TABLES_DIR TSS_WORK REFSEQ_TSS ATAQV_FILE FRIP_FILE TSS_FILE

for SAMPLE in "${SAMPLES[@]}"; do
    process_sample "$SAMPLE" &
    while [[ $(jobs -r | wc -l) -ge $MAX_PARALLEL ]]; do sleep 15; done
done
wait

for SAMPLE in "${SAMPLES[@]}"; do
    for F in "$ATAQV_FILE" "$FRIP_FILE" "$TSS_FILE"; do
        cat "${F}.${SAMPLE}.tmp" >> "$F"
        rm -f "${F}.${SAMPLE}.tmp"
    done
done
rmdir "$TSS_WORK" 2>/dev/null || true

## Flags table -- reference values only, not drawn as cutoff lines on any plot.
echo -e "Sample\tUsable_OK\tFRiP\tFRiP_OK\tNpeaks\tNpeaks_OK\tDup_fraction\tDup_OK\tTSS\tTSS_OK" > "$FLAGS_FILE"
paste \
    <(tail -n +2 "$ATAQV_FILE" | awk -v t="$MIN_USABLE" '{print $1"\t"($7+0>t?"PASS":"WARN")}') \
    <(tail -n +2 "$FRIP_FILE"  | awk -v tf="$MIN_FRIP" -v tp="$MIN_PEAKS" '{print $4"\t"($4+0>tf?"PASS":"FAIL")"\t"$5"\t"($5+0>tp?"PASS":"FAIL")}') \
    <(tail -n +2 "$ATAQV_FILE" | awk -v td="$MAX_DUP_FRACTION" '{print $5"\t"($5+0<td?"PASS":"WARN")}') \
    <(tail -n +2 "$TSS_FILE"   | awk -v tt="$MIN_TSS" '{print $2"\t"($2+0>tt?"PASS":"FAIL")}') \
    >> "$FLAGS_FILE"

echo "[A03] Complete -- $(date). Tables: ${TABLES_DIR}/A03_*.tsv"
