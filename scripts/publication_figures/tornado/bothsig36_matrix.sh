#!/usr/bin/env bash
set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
AK="$REPO_ROOT/data/external/akiyama_chip"
ED="$REPO_ROOT/data/external/encode_h1hesc"
SCR="$REPO_ROOT/results/atac/tornado_intermediate"
cd "$REPO_ROOT"
source /storage/volume01/emile/miniconda3/etc/profile.d/conda.sh
conda activate viz_env

mapfile -t BEDS < <(cut -f2 "$SCR/bothsig36_order.tsv" | sed "s#^#$SCR/bothsig36_beds/#; s#\$#.bed#")
computeMatrix reference-point --referencePoint center \
  -R "${BEDS[@]}" \
  -S "$AK/MLL4_hES_rep_1.bigWig" "$AK/UTX_hES_avg.bigWig" "$AK/Flag_hES_ctrl.bigWig" \
     "$ED/H3K4me1.bigWig" "$ED/H3K27ac.bigWig" "$ED/H3K27me3.bigWig" \
  --samplesLabel KMT2D KDM6A FLAG H3K4me1 H3K27ac H3K27me3 \
  -b 3000 -a 3000 --binSize 50 --missingDataAsZero \
  -p 12 -o "$SCR/matrix_bothsig36.gz" 2> "$SCR/bothsig36_cm.err"
echo "[bothsig36] DONE rc=$? -> $SCR/matrix_bothsig36.gz"
