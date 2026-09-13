#!/usr/bin/env bash
# A08_ATAC_Viz.sh -- tornado plots (deepTools heatmaps) over Loess-normalized
# DARs (A04b/A04f) and their union/shared/TSS-distance splits. TMM (A04) is
# kept only as the comparison baseline that justifies preferring Loess (see
# docs/A04b_normalization_methodology.md), not fed into this or any other
# downstream ATAC step.
# Adapted from ED26_001_ATAC_c/Narrow_allReads/scripts/06_tornado_replot.sh.
#
# Self-checkpointing: skips entirely if results/atac/A08_union_DARs_by_distance.pdf
# (the last plot produced) already exists. Individual intermediate files (merged
# bigwigs, BED files, matrices) are each skipped if already present.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SCREEN_SESSION="ks1_2_A08"
if [[ -z "${STY:-}" ]]; then
    echo "[A08] Not inside screen -- relaunching in screen session '${SCREEN_SESSION}' (reattach with: screen -r ${SCREEN_SESSION})."
    mkdir -p logs
    exec screen -L -Logfile "logs/A08_screen_$(date +%Y%m%d_%H%M%S).log" -S "${SCREEN_SESSION}" bash "$0" "$@"
fi

CONFIG="config/pipeline_config.yaml"
METADATA="data/sample_metadata.csv"
BW_DIR="data/atac_bigwigs"
DAR_DIR="data/diffbind_output/CsawLoess_Norm"
WORK_DIR="data/tornado_output"
PLOTS_DIR="results/atac"
LOGDIR="logs"
LAST_PLOT="${PLOTS_DIR}/A08_union_DARs_by_distance.pdf"
mkdir -p "$LOGDIR" "$WORK_DIR" "$PLOTS_DIR"

if [[ -f "$LAST_PLOT" ]]; then
    echo "[A08] Already complete (${LAST_PLOT} exists). Skipping. Delete that file to force a rerun."
    exit 0
fi

[[ -f "$CONFIG" ]]   || { echo "[A08] ERROR: missing ${CONFIG}."; exit 1; }
[[ -f "$METADATA" ]] || { echo "[A08] ERROR: missing ${METADATA}."; exit 1; }
[[ -d "$BW_DIR" ]]   || { echo "[A08] ERROR: missing ${BW_DIR}. Run A02_MACS3_narrow.sh first."; exit 1; }
[[ -d "$DAR_DIR" ]]  || { echo "[A08] ERROR: missing ${DAR_DIR}. Run A04b_Csaw_Loess_Norm.R and A04f_Loess_DAR_BED.R first."; exit 1; }

CONDA_BASE="${KS12_CONDA_BASE:-/storage/volume01/emile/miniconda3}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate nf-run   # provides yq -- viz_env does not have it

TSS="${REPO_ROOT}/$(yq -r '.reference.tss_bed' "$CONFIG")"
CORES=$(yq -r '.resources.tornado_cores' "$CONFIG")
# Genome file (chrom order) for bedtools sort/closest -- karyotypic order
# (chr1,chr2,...,chr10,...), not plain lexicographic ("chr1","chr10","chr2",...
# sort -k1,1 alone would produce). bedtools closest -d requires -a/-b to share
# one consistent order; the reference .fai already has it and doubles as a
# valid bedtools genome file (chrom, length as its first two columns).
FASTA_GZ="${REPO_ROOT}/$(yq -r '.reference.fasta' "$CONFIG")"
FAI="${FASTA_GZ%.gz}.fai"

N_CONTRASTS=$(yq -r '.contrasts | length' "$CONFIG")
CONTRAST_NAMES=(); CONTRAST_TREATMENTS=(); CONTRAST_REFERENCES=()
for ((i=0; i<N_CONTRASTS; i++)); do
    CONTRAST_NAMES+=("$(yq -r ".contrasts[$i].name" "$CONFIG")")
    CONTRAST_TREATMENTS+=("$(yq -r ".contrasts[$i].treatment" "$CONFIG")")
    CONTRAST_REFERENCES+=("$(yq -r ".contrasts[$i].reference" "$CONFIG")")
done
conda deactivate

conda activate viz_env
[[ -f "$TSS" ]] || { echo "[A08] ERROR: missing ${TSS}."; exit 1; }
[[ -f "$FAI" ]] || { echo "[A08] ERROR: missing ${FAI} (genome file for bedtools sort/closest)."; exit 1; }

FLANK=2000
BIN=50
HEIGHT=30
WIDTH=15
# Sequential, not diverging: RdYlBu_r's low end is a solid medium blue, not
# pale, so faint/background regions never actually looked faint. YlOrRd runs
# near-white (low) -> yellow -> orange -> red (high), same "strong = red"
# reading as before but with real visual contrast between faint and strong.
CMAP=YlOrRd
# TSS-distance cutoff for the proximal/distal split (A08_union_DARs_by_distance.pdf
# only) -- widened from an earlier 2kb, which put only 331/7471 (4.4%) of union
# GAINED DARs into that bin, making it too sparse to read as a heatmap panel.
# Not "promoter": bedtools closest reports distance to the nearest TSS, which
# for many of these DARs is not necessarily THE gene's own promoter -- "proximal"
# says only what's actually being measured.
PROXIMAL_DIST=5000
PROXIMAL_LABEL="$(( PROXIMAL_DIST / 1000 ))kb"
# Fixed z-scale across all sample columns (KDM6A_ko/KMT2D_Het/WT) so signal
# intensity is directly comparable between genotypes -- without this,
# plotHeatmap auto-scales each column to its own max, which can visually
# flatten a real difference in accessibility between genotypes.
# Z_MAX=4 was ~the 95th percentile of union-matrix bin values (already
# clipping ~5% of bins to full saturation). Widened to 12 after comparing
# 4/6/8/10/12 side by side across all 3 genotypes on the same region set
# (KMT2D_Het_vs_WT GAINED): at this range WT (genuinely weaker signal at
# GAINED DARs, by construction) fades furthest while KDM6A_ko/KMT2D_Het still
# retain visible structure -- the widened scale doesn't just make everything
# fainter, it makes the real WT-vs-mutant accessibility gap read more clearly.
# True max in this data is 244 (a rare outlier bin), so scaling to raw max
# was never an option.
Z_MIN=0
Z_MAX=12
XLABEL="distance from peak center (bp)"

########
## STEP 0: per-genotype merged bigwig + per-genotype sample bigwig list
########
echo "[A08] Merging bigwigs per genotype..."
mapfile -t GENOTYPES < <(tail -n +2 "$METADATA" | awk -F',' '$3=="ATAC"{print $4}' | sort -u)
# WT first, everything else in its existing (alphabetical) order after it --
# all downstream ALL_MERGED_BWS/ALL_MERGED_LABELS arrays are built by looping
# over GENOTYPES, so this one reorder makes WT the first column (and hence
# deepTools sample index 1) everywhere merged-genotype heatmaps are plotted.
REORDERED_GENOTYPES=("WT")
for g in "${GENOTYPES[@]}"; do [[ "$g" != "WT" ]] && REORDERED_GENOTYPES+=("$g"); done
GENOTYPES=("${REORDERED_GENOTYPES[@]}")

declare -A GENO_BWS   # space-separated per-sample bigwig paths, one entry per genotype
declare -A GENO_LABELS

for GENO in "${GENOTYPES[@]}"; do
    mapfile -t SAMPLES < <(tail -n +2 "$METADATA" | awk -F',' -v g="$GENO" '$3=="ATAC" && $4==g{print $1}')
    BWS=()
    LABELS=()
    for s in "${SAMPLES[@]}"; do BWS+=("${BW_DIR}/${s}.bw"); LABELS+=("$s"); done
    GENO_BWS["$GENO"]="${BWS[*]}"
    GENO_LABELS["$GENO"]="${LABELS[*]}"

    MERGED="${WORK_DIR}/${GENO}_merged.bw"
    if [[ ! -f "$MERGED" ]]; then
        bigwigAverage -b "${BWS[@]}" -o "$MERGED" -p "$CORES"
    fi
done

########
## STEP 1: union / shared / TSS-distance BEDs, across all configured contrasts
########
echo "[A08] Building union/shared/distance BED files..."

UNION_GAINED_PARTS=()
UNION_LOST_PARTS=()
for ((i=0; i<N_CONTRASTS; i++)); do
    NAME="${CONTRAST_NAMES[$i]}"
    UNION_GAINED_PARTS+=("${DAR_DIR}/${NAME}_GAINED.bed")
    UNION_LOST_PARTS+=("${DAR_DIR}/${NAME}_LOST.bed")
done

cat "${UNION_GAINED_PARTS[@]}" | bedtools sort -g "$FAI" -i - | bedtools merge > "${WORK_DIR}/union_GAINED.bed"
cat "${UNION_LOST_PARTS[@]}"   | bedtools sort -g "$FAI" -i - | bedtools merge > "${WORK_DIR}/union_LOST.bed"

if [[ "$N_CONTRASTS" -ge 2 ]]; then
    NAME0="${CONTRAST_NAMES[0]}"
    NAME1="${CONTRAST_NAMES[1]}"
    bedtools intersect -a "${DAR_DIR}/${NAME0}_GAINED.bed" -b "${DAR_DIR}/${NAME1}_GAINED.bed" -u > "${WORK_DIR}/shared_GAINED.bed"
    bedtools intersect -a "${DAR_DIR}/${NAME0}_LOST.bed"   -b "${DAR_DIR}/${NAME1}_LOST.bed"   -u > "${WORK_DIR}/shared_LOST.bed"
fi

TSS_SORTED="${WORK_DIR}/tss_sorted.bed"
bedtools sort -g "$FAI" -i "$TSS" > "$TSS_SORTED"

bedtools closest -g "$FAI" -a "${WORK_DIR}/union_GAINED.bed" -b "$TSS_SORTED" -d | awk -v OFS="\t" -v t="$PROXIMAL_DIST" '$NF <= t {print $1,$2,$3}' | sort -u > "${WORK_DIR}/union_GAINED_proximal.bed"
bedtools closest -g "$FAI" -a "${WORK_DIR}/union_GAINED.bed" -b "$TSS_SORTED" -d | awk -v OFS="\t" -v t="$PROXIMAL_DIST" '$NF >  t {print $1,$2,$3}' | sort -u > "${WORK_DIR}/union_GAINED_distal.bed"
bedtools closest -g "$FAI" -a "${WORK_DIR}/union_LOST.bed"   -b "$TSS_SORTED" -d | awk -v OFS="\t" -v t="$PROXIMAL_DIST" '$NF <= t {print $1,$2,$3}' | sort -u > "${WORK_DIR}/union_LOST_proximal.bed"
bedtools closest -g "$FAI" -a "${WORK_DIR}/union_LOST.bed"   -b "$TSS_SORTED" -d | awk -v OFS="\t" -v t="$PROXIMAL_DIST" '$NF >  t {print $1,$2,$3}' | sort -u > "${WORK_DIR}/union_LOST_distal.bed"

echo "[A08] Union GAINED: $(wc -l < "${WORK_DIR}/union_GAINED.bed")  Union LOST: $(wc -l < "${WORK_DIR}/union_LOST.bed")"

plot_tornado() {
    local NAME="$1" MATRIX="$2" OUT_PDF="$3"; shift 3
    local SAMPLES_LABEL=("$@")
    plotHeatmap -m "$MATRIX" -o "$OUT_PDF" \
        --colorMap "$CMAP" --sortUsing mean --sortUsingSamples 1 --sortRegions descend \
        --zMin "$Z_MIN" --zMax "$Z_MAX" --xAxisLabel "$XLABEL" \
        --samplesLabel "${SAMPLES_LABEL[@]}" \
        --regionsLabel "GAINED" "LOST" \
        --heatmapHeight "$HEIGHT" --heatmapWidth 3 \
        --whatToShow "heatmap and colorbar"
    echo "[A08] Saved: ${OUT_PDF}"
}

########
## STEP 2: per-contrast tornados (reference vs treatment, per-sample bigwigs)
########
for ((i=0; i<N_CONTRASTS; i++)); do
    NAME="${CONTRAST_NAMES[$i]}"
    TREATMENT="${CONTRAST_TREATMENTS[$i]}"
    REFERENCE="${CONTRAST_REFERENCES[$i]}"
    GAINED="${DAR_DIR}/${NAME}_GAINED.bed"
    LOST="${DAR_DIR}/${NAME}_LOST.bed"

    echo "[A08] Per-sample tornado: ${NAME}"
    MATRIX="${WORK_DIR}/matrix_${NAME}_persample.gz"
    read -ra REF_BWS <<< "${GENO_BWS[$REFERENCE]}"
    read -ra TR_BWS  <<< "${GENO_BWS[$TREATMENT]}"
    read -ra REF_LBL <<< "${GENO_LABELS[$REFERENCE]}"
    read -ra TR_LBL  <<< "${GENO_LABELS[$TREATMENT]}"
    computeMatrix reference-point \
        -S "${REF_BWS[@]}" "${TR_BWS[@]}" \
        -R "$GAINED" "$LOST" \
        --referencePoint center -a "$FLANK" -b "$FLANK" --binSize "$BIN" \
        --missingDataAsZero -o "$MATRIX" -p "$CORES"
    plot_tornado "$NAME" "$MATRIX" "${PLOTS_DIR}/A08_${NAME}_persample.pdf" "${REF_LBL[@]}" "${TR_LBL[@]}"

    echo "[A08] Merged-genotype tornado: ${NAME}"
    MERGED_MATRIX="${WORK_DIR}/matrix_${NAME}_merged.gz"
    ALL_MERGED_BWS=()
    ALL_MERGED_LABELS=()
    for GENO in "${GENOTYPES[@]}"; do ALL_MERGED_BWS+=("${WORK_DIR}/${GENO}_merged.bw"); ALL_MERGED_LABELS+=("$GENO"); done
    computeMatrix reference-point \
        -S "${ALL_MERGED_BWS[@]}" \
        -R "$GAINED" "$LOST" \
        --referencePoint center -a "$FLANK" -b "$FLANK" --binSize "$BIN" \
        --missingDataAsZero -o "$MERGED_MATRIX" -p "$CORES"
    plotHeatmap -m "$MERGED_MATRIX" -o "${PLOTS_DIR}/A08_${NAME}_merged.pdf" \
        --colorMap "$CMAP" --sortUsing mean --sortUsingSamples 1 --sortRegions descend \
        --zMin "$Z_MIN" --zMax "$Z_MAX" --xAxisLabel "$XLABEL" \
        --samplesLabel "${ALL_MERGED_LABELS[@]}" \
        --regionsLabel "GAINED" "LOST" \
        --heatmapHeight "$HEIGHT" --heatmapWidth "$WIDTH" \
        --whatToShow "heatmap and colorbar"
    echo "[A08] Saved: ${PLOTS_DIR}/A08_${NAME}_merged.pdf"

    echo "[A08] Merged-genotype tornado by TSS distance: ${NAME} (proximal <=${PROXIMAL_LABEL} / distal >${PROXIMAL_LABEL})"
    GAINED_SORTED="${WORK_DIR}/${NAME}_GAINED_sorted.bed"
    LOST_SORTED="${WORK_DIR}/${NAME}_LOST_sorted.bed"
    bedtools sort -g "$FAI" -i "$GAINED" > "$GAINED_SORTED"
    bedtools sort -g "$FAI" -i "$LOST"   > "$LOST_SORTED"
    bedtools closest -g "$FAI" -a "$GAINED_SORTED" -b "$TSS_SORTED" -d | awk -v OFS="\t" -v t="$PROXIMAL_DIST" '$NF <= t {print $1,$2,$3}' | sort -u > "${WORK_DIR}/${NAME}_GAINED_proximal.bed"
    bedtools closest -g "$FAI" -a "$GAINED_SORTED" -b "$TSS_SORTED" -d | awk -v OFS="\t" -v t="$PROXIMAL_DIST" '$NF >  t {print $1,$2,$3}' | sort -u > "${WORK_DIR}/${NAME}_GAINED_distal.bed"
    bedtools closest -g "$FAI" -a "$LOST_SORTED"   -b "$TSS_SORTED" -d | awk -v OFS="\t" -v t="$PROXIMAL_DIST" '$NF <= t {print $1,$2,$3}' | sort -u > "${WORK_DIR}/${NAME}_LOST_proximal.bed"
    bedtools closest -g "$FAI" -a "$LOST_SORTED"   -b "$TSS_SORTED" -d | awk -v OFS="\t" -v t="$PROXIMAL_DIST" '$NF >  t {print $1,$2,$3}' | sort -u > "${WORK_DIR}/${NAME}_LOST_distal.bed"
    echo "[A08] ${NAME} GAINED: $(wc -l < "${WORK_DIR}/${NAME}_GAINED_proximal.bed") proximal / $(wc -l < "${WORK_DIR}/${NAME}_GAINED_distal.bed") distal | LOST: $(wc -l < "${WORK_DIR}/${NAME}_LOST_proximal.bed") proximal / $(wc -l < "${WORK_DIR}/${NAME}_LOST_distal.bed") distal"

    DIST_MATRIX="${WORK_DIR}/matrix_${NAME}_by_distance.gz"
    computeMatrix reference-point \
        -S "${ALL_MERGED_BWS[@]}" \
        -R "${WORK_DIR}/${NAME}_GAINED_proximal.bed" "${WORK_DIR}/${NAME}_GAINED_distal.bed" \
           "${WORK_DIR}/${NAME}_LOST_proximal.bed" "${WORK_DIR}/${NAME}_LOST_distal.bed" \
        --referencePoint center -a "$FLANK" -b "$FLANK" --binSize "$BIN" \
        --missingDataAsZero -o "$DIST_MATRIX" -p "$CORES"
    plotHeatmap -m "$DIST_MATRIX" -o "${PLOTS_DIR}/A08_${NAME}_by_distance.pdf" \
        --colorMap "$CMAP" --sortUsing mean --sortUsingSamples 1 --sortRegions descend \
        --zMin "$Z_MIN" --zMax "$Z_MAX" --xAxisLabel "$XLABEL" \
        --samplesLabel "${ALL_MERGED_LABELS[@]}" \
        --regionsLabel "GAINED proximal (<=${PROXIMAL_LABEL})" "GAINED distal (>${PROXIMAL_LABEL})" \
                       "LOST proximal (<=${PROXIMAL_LABEL})" "LOST distal (>${PROXIMAL_LABEL})" \
        --heatmapHeight "$HEIGHT" --heatmapWidth "$WIDTH" \
        --whatToShow "heatmap and colorbar"
    echo "[A08] Saved: ${PLOTS_DIR}/A08_${NAME}_by_distance.pdf"
done

########
## STEP 3: shared / union / union-by-distance tornados (merged genotypes, all contrasts)
########
ALL_MERGED_BWS=()
ALL_MERGED_LABELS=()
for GENO in "${GENOTYPES[@]}"; do ALL_MERGED_BWS+=("${WORK_DIR}/${GENO}_merged.bw"); ALL_MERGED_LABELS+=("$GENO"); done

if [[ -s "${WORK_DIR}/shared_GAINED.bed" || -s "${WORK_DIR}/shared_LOST.bed" ]]; then
    echo "[A08] Shared-DAR tornado"
    computeMatrix reference-point -S "${ALL_MERGED_BWS[@]}" \
        -R "${WORK_DIR}/shared_GAINED.bed" "${WORK_DIR}/shared_LOST.bed" \
        --referencePoint center -a "$FLANK" -b "$FLANK" --binSize "$BIN" \
        --missingDataAsZero -o "${WORK_DIR}/matrix_shared.gz" -p "$CORES"
    plotHeatmap -m "${WORK_DIR}/matrix_shared.gz" -o "${PLOTS_DIR}/A08_shared_DARs_merged.pdf" \
        --colorMap "$CMAP" --sortUsing mean --sortUsingSamples 1 --sortRegions descend \
        --zMin "$Z_MIN" --zMax "$Z_MAX" --xAxisLabel "$XLABEL" \
        --samplesLabel "${ALL_MERGED_LABELS[@]}" --regionsLabel "Shared GAINED" "Shared LOST" \
        --heatmapHeight "$HEIGHT" --heatmapWidth "$WIDTH" --whatToShow "heatmap and colorbar"
    echo "[A08] Saved: ${PLOTS_DIR}/A08_shared_DARs_merged.pdf"
fi

echo "[A08] Union-DAR tornado"
computeMatrix reference-point -S "${ALL_MERGED_BWS[@]}" \
    -R "${WORK_DIR}/union_GAINED.bed" "${WORK_DIR}/union_LOST.bed" \
    --referencePoint center -a "$FLANK" -b "$FLANK" --binSize "$BIN" \
    --missingDataAsZero -o "${WORK_DIR}/matrix_union.gz" -p "$CORES"
plotHeatmap -m "${WORK_DIR}/matrix_union.gz" -o "${PLOTS_DIR}/A08_union_DARs_merged.pdf" \
    --colorMap "$CMAP" --sortUsing mean --sortUsingSamples 1 --sortRegions descend \
    --zMin "$Z_MIN" --zMax "$Z_MAX" --xAxisLabel "$XLABEL" \
    --samplesLabel "${ALL_MERGED_LABELS[@]}" --regionsLabel "GAINED" "LOST" \
    --heatmapHeight "$HEIGHT" --heatmapWidth "$WIDTH" --whatToShow "heatmap and colorbar"
echo "[A08] Saved: ${PLOTS_DIR}/A08_union_DARs_merged.pdf"

echo "[A08] Union-DAR-by-TSS-distance tornado (proximal <=${PROXIMAL_LABEL} / distal >${PROXIMAL_LABEL})"
computeMatrix reference-point -S "${ALL_MERGED_BWS[@]}" \
    -R "${WORK_DIR}/union_GAINED_proximal.bed" "${WORK_DIR}/union_GAINED_distal.bed" \
       "${WORK_DIR}/union_LOST_proximal.bed" "${WORK_DIR}/union_LOST_distal.bed" \
    --referencePoint center -a "$FLANK" -b "$FLANK" --binSize "$BIN" \
    --missingDataAsZero -o "${WORK_DIR}/matrix_union_by_distance.gz" -p "$CORES"
plotHeatmap -m "${WORK_DIR}/matrix_union_by_distance.gz" -o "$LAST_PLOT" \
    --colorMap "$CMAP" --sortUsing mean --sortUsingSamples 1 --sortRegions descend \
    --zMin "$Z_MIN" --zMax "$Z_MAX" --xAxisLabel "$XLABEL" \
    --samplesLabel "${ALL_MERGED_LABELS[@]}" \
    --regionsLabel "GAINED proximal (<=${PROXIMAL_LABEL})" "GAINED distal (>${PROXIMAL_LABEL})" \
                   "LOST proximal (<=${PROXIMAL_LABEL})" "LOST distal (>${PROXIMAL_LABEL})" \
    --heatmapHeight "$HEIGHT" --heatmapWidth "$WIDTH" --whatToShow "heatmap and colorbar"
echo "[A08] Saved: ${LAST_PLOT}"
echo "[A08] Saved: ${LAST_PLOT}"

conda deactivate
echo ""
echo "[A08] Complete -- $(date). Plots: ${PLOTS_DIR}/A08_*.pdf"
