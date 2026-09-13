#!/usr/bin/env bash
# 36 TFs significant in BOTH contrasts (A07b_heatmap_both_sig_with_motifs):
# resolve each symbol to its best JASPAR matrix (most high-conf WT-bound
# sites), then take that motif's top 200 high-conf WT-bound instances by
# TOBIAS bound score. Ordered by TF family, then symbol.
#
# Requires results/atac/tornado_intermediate/bound_highconf_canon.bed --
# unfiltered (all chromosomes) equivalent is fine too, see
# 02_define_bound_highconf_canon.sh for how that file (and its own upstream
# gap) is documented.
set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
HC="$REPO_ROOT/data/tobias_output/all_WT_bound_highconf.bed"
SCR="$REPO_ROOT/results/atac/tornado_intermediate"
OUT="$SCR/bothsig36_beds"
BD="$REPO_ROOT/data/tobias_output/KMT2D_Het_vs_WT/BINDetect"
N=200
mkdir -p "$SCR"
rm -rf "$OUT"; mkdir -p "$OUT"
: > "$SCR/bothsig36_order.tsv"   # family<TAB>symbol<TAB>motif<TAB>nsites

# family : space-separated symbols  (order from A07b both_sig_family)
declare -A FAM=(
  ["AP-1 (Fos/Jun)"]="FOS FOSL1 FOSL2 JDP2 JUNB JUND"
  ["bZIP (BACH/Maf)"]="BACH1 BACH2 NFE2"
  ["FOX (Forkhead)"]="FOXA3 FOXD1 FOXK1 FOXK2 FOXL1 FOXO6 FOXP4"
  ["Nuclear receptor"]="ESRRA ESRRB NR2F2 NR6A1"
  ["RFX"]="RFX1 RFX2 RFX3 RFX5"
  ["SOX (HMG-box)"]="SOX2 SOX4 SOX8 SOX9 SOX10 SOX13 SOX15"
  ["Zinc finger"]="BNC2 KLF15 ZBED4 ZNF610"
  ["Homeodomain"]="PAX3"
)
FAM_ORDER=("AP-1 (Fos/Jun)" "bZIP (BACH/Maf)" "FOX (Forkhead)" "Nuclear receptor" "RFX" "SOX (HMG-box)" "Zinc finger" "Homeodomain")

for fam in "${FAM_ORDER[@]}"; do
  for s in ${FAM[$fam]}; do
    # candidate motifs for this symbol
    mapfile -t cands < <(ls -d "$BD"/*/ 2>/dev/null | xargs -n1 basename \
        | awk -v s="$s" 'BEGIN{IGNORECASE=1}{n=$0; sub(/_MA[0-9].*/,"",n); if (toupper(n)==toupper(s)) print}')
    best=""; bestn=-1
    for m in "${cands[@]}"; do
      c=$(awk -F'\t' -v m="$m" '$4==m && $1 ~ /^chr[0-9XY]+$/' "$HC" | wc -l)
      if [ "$c" -gt "$bestn" ]; then bestn=$c; best=$m; fi
    done
    awk -F'\t' -v m="$best" '$4==m && $1 ~ /^chr[0-9XY]+$/' "$HC" \
      | sort -t$'\t' -k5,5gr -k1,1 -k2,2n | head -n "$N" \
      | cut -f1-3 | sort -k1,1 -k2,2n > "$OUT/${s}.bed"
    kept=$(wc -l < "$OUT/${s}.bed")
    printf '%s\t%s\t%s\t%s\n' "$fam" "$s" "$best" "$kept" >> "$SCR/bothsig36_order.tsv"
  done
done
echo "[bothsig36] $(ls $OUT | wc -l) TF beds, $(cat $OUT/*.bed | wc -l) sites"
column -t -s$'\t' "$SCR/bothsig36_order.tsv"
