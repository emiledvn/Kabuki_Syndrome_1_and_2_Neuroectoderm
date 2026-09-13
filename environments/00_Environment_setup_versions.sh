#!/usr/bin/env bash
# 00_Environment_setup_versions.sh -- checks that all conda envs this pipeline
# needs already exist (nf-run, atac-narrow, viz_env, tobias_env, ks_1_2_r) and
# prints their pinned tool versions. Does NOT auto-create missing envs -- they're
# either shared infrastructure (nf-run, atac-narrow, viz_env, tobias_env, used by
# other analysis/ED2x_* projects too) or a deliberate one-time clone+extend
# (ks_1_2_r, see README "Environments"), so recreating them here would either
# duplicate slow NFS work already done, or silently diverge from the exact pinned
# versions in envs/*.yml.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# Verified-at-runtime versions are only useful for reproducibility if they're
# kept, not just printed to a screen log that gets rotated away -- mirrors the
# R scripts' own results/Session_info/*_session_info.txt convention.
SESSION_DIR="results/Session_info"
mkdir -p "$SESSION_DIR"
OUTFILE="${SESSION_DIR}/00_environment_versions.txt"
exec > >(tee "$OUTFILE") 2>&1

CONDA_BASE="${KS12_CONDA_BASE:-/storage/volume01/emile/miniconda3}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

REQUIRED_ENVS=(nf-run atac-narrow viz_env tobias_env ks_1_2_r)
MISSING=()

echo "[00] Checking required conda environments..."
EXISTING=$(conda env list 2>/dev/null | awk '{print $1}')

for env in "${REQUIRED_ENVS[@]}"; do
    if echo "$EXISTING" | grep -qx "$env"; then
        echo "[00] OK: $env"
    else
        echo "[00] MISSING: $env"
        MISSING+=("$env")
    fi
done

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo ""
    echo "[00] ERROR: ${#MISSING[@]} required env(s) missing: ${MISSING[*]}"
    echo "[00] Not auto-created (see script header). To restore from the pinned spec:"
    for env in "${MISSING[@]}"; do
        echo "[00]   conda env create -n ${env} -f envs/${env}.yml"
    done
    exit 1
fi

echo ""
echo "[00] All required environments present. Key tool versions:"
echo "--- nf-run ---"
conda run -n nf-run nextflow -version 2>/dev/null | grep -i version
conda run -n nf-run apptainer --version 2>/dev/null
echo "--- atac-narrow ---"
conda run -n atac-narrow macs3 --version 2>&1
echo "--- viz_env ---"
conda run -n viz_env bamCoverage --version 2>&1
echo "--- tobias_env ---"
conda run -n tobias_env TOBIAS --version 2>&1 | head -1
echo "--- ks_1_2_r ---"
conda run -n ks_1_2_r Rscript -e 'cat("R", R.version.string, "\n"); cat("DESeq2", as.character(packageVersion("DESeq2")), "\n")' 2>/dev/null

echo ""
echo "[00] Complete -- $(date)"
