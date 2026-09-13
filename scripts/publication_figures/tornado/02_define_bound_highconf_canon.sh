#!/usr/bin/env bash
# 02_define_bound_highconf_canon.sh -- canonical-chromosome filter over the
# merged, high-confidence WT-bound TOBIAS footprint set, producing
# bound_highconf_canon.bed for Panel B (footprint-vs-background binding map).
#
# KNOWN GAP: this script expects data/tobias_output/all_WT_bound_highconf.bed
# to already exist, but no script in this repo currently GENERATES that file
# -- it is a merge, across all TFs, of TOBIAS BINDetect's per-TF high-
# confidence WT-bound sites (5-column BED: chrom, start, end, motif_id,
# bound-score), produced ad hoc during exploration and never captured as a
# script. Regenerating it means concatenating each TF's high-confidence
# bound-site BED from data/tobias_output/<contrast>/BINDetect/<TF>/ (see
# scripts/atac/A07b_TOBIAS_DOWNSTREAM.R for how BINDetect output is read and
# what "high-confidence" thresholds it already applies elsewhere) and
# sorting the concatenation. That merge script does not exist yet -- until
# it does, this step needs data/tobias_output/all_WT_bound_highconf.bed
# supplied manually.
#
# What this script itself does is trivial and fully deterministic: keep only
# canonical chromosomes (chr1-22, chrX, chrY), same filter already applied
# throughout this chain (e.g. bothsig36_beds.sh).
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"

IN="data/tobias_output/all_WT_bound_highconf.bed"
OUT_DIR="results/atac/tornado_intermediate"
OUT="${OUT_DIR}/bound_highconf_canon.bed"
mkdir -p "$OUT_DIR"

if [ ! -f "$IN" ]; then
  echo "[02] ERROR: ${IN} not found. See this script's header -- it must be" >&2
  echo "     supplied manually (or regenerated) before this step can run." >&2
  exit 1
fi

awk -F'\t' '$1 ~ /^chr[0-9XY]+$/' "$IN" | sort -k1,1 -k2,2n > "$OUT"
echo "[02] $(wc -l < "$OUT") canonical-chromosome bound sites -> ${OUT}"
