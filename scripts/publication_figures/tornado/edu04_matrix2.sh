#!/usr/bin/env bash
set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
AK="$REPO_ROOT/data/external/akiyama_chip"
ED="$REPO_ROOT/data/external/encode_h1hesc"
SCR="$REPO_ROOT/results/atac/tornado_intermediate"
cd "$REPO_ROOT"
source /storage/volume01/emile/miniconda3/etc/profile.d/conda.sh
conda activate viz_env

# bg_regions.bed (from 03_define_background_regions.R) replaces the original,
# unrecoverable bg_fixed.bed -- see that script's header for why.
computeMatrix reference-point --referencePoint center \
  -R "$SCR/bound_highconf_canon.bed" "$SCR/bg_regions.bed" \
  -S "$AK/MLL4_hES_rep_1.bigWig" "$AK/UTX_hES_avg.bigWig" "$AK/Flag_hES_ctrl.bigWig" \
     "$ED/H3K4me1.bigWig" "$ED/H3K27ac.bigWig" "$ED/H3K27me3.bigWig" \
  --samplesLabel KMT2D KDM6A Ctrl H3K4me1 H3K27ac H3K27me3 \
  -b 3000 -a 3000 --binSize 50 --missingDataAsZero \
  -p 12 -o "$SCR/matrix_edu04_6col.gz" 2> "$SCR/edu04_cm2.err"
echo "[edu04] DONE matrix rc=$? -> $SCR/matrix_edu04_6col.gz"
