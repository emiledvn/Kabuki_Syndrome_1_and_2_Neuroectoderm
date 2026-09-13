#!/usr/bin/env bash
# 00_fetch_encode_tracks.sh -- downloads the 3 public ENCODE H1-hESC
# fold-change-over-control bigWigs used as reference tracks in the Figure4
# tornado plot (H3K4me1, H3K27ac, H3K27me3). Same accessions already used by
# scripts/EDU02_Enhancer_Tornado.sh for the enhancer-state illustration --
# reused verbatim here, not re-derived.
#
# Self-checkpointing: skips a file if it already exists.
#
# Run: bash scripts/publication_figures/tornado/00_fetch_encode_tracks.sh
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT_DIR="$REPO_ROOT/data/external/encode_h1hesc"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

declare -A BW_URL=(
  [H3K4me1]="https://encode-public.s3.amazonaws.com/2023/01/15/5d2eef5e-9d65-4455-81eb-c5366d96a3f5/ENCFF396RXV.bigWig"
  [H3K27ac]="https://encode-public.s3.amazonaws.com/2021/02/01/b55106ed-b1c1-44ec-98c1-0b4f6923c7ab/ENCFF919FBG.bigWig"
  [H3K27me3]="https://encode-public.s3.amazonaws.com/2023/01/15/3651400d-7969-470e-8b8b-4a3eb0f9fb43/ENCFF380KPI.bigWig"
)

for mark in H3K4me1 H3K27ac H3K27me3; do
  f="${mark}.bigWig"
  if [ -f "$f" ]; then
    echo "[00] ${f} already present, skipping."
  else
    echo "[00] Fetching ${mark} (${BW_URL[$mark]})"
    curl -sL "${BW_URL[$mark]}" -o "${f}.part" && mv "${f}.part" "$f"
  fi
done
echo "[00] Done -- 3 ENCODE H1-hESC bigWigs in ${OUT_DIR}"
