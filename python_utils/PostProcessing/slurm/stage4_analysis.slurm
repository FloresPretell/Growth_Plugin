#!/bin/bash
#SBATCH --job-name=Morphoanalysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH -o logs/stage4_%j.out
#SBATCH -e logs/stage4_%j.err

set -euo pipefail

echo "=== Stage 4 Multivariate analysis ==="
echo "Timestamp      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Hostname       : $(hostname)"
echo "Working Dir    : $(pwd)"
echo "Pipeline Root  : ${PIPELINE_ROOT:-<unset>}"
echo "Morpho Root    : ${MORPHO_ROOT:-<unset>}"
echo "Analysis Root  : ${ANALYSIS_ROOT:-<unset>}"

CSV_PATH="${MORPHO_ROOT}/morphometrics.csv"
[[ -f "${CSV_PATH}" ]] || { echo "ERROR: morphometrics.csv not found: ${CSV_PATH}" >&2; exit 1; }
[[ -n "${MORPHO_PYTHON:-}" ]] || { echo "ERROR: MORPHO_PYTHON not set" >&2; exit 1; }

mkdir -p logs "${ANALYSIS_ROOT}"

# === Clean ParaView environment pollution ===
unset PYTHONPATH
source /scratch/flore0a/iops/miniconda3/etc/profile.d/conda.sh
conda activate swc_env

# Verify conda Python works
"${MORPHO_PYTHON}" --version || { echo "ERROR: MORPHO_PYTHON failed: ${MORPHO_PYTHON}" >&2; exit 1; }

cmd=("${MORPHO_PYTHON}" "${PIPELINE_ROOT}/Morphoanalysis_plot/analyze_morphometrics.py" \
     --input "${CSV_PATH}" \
     --output-dir "${ANALYSIS_ROOT}")

printf 'Command        :'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

echo "[success] Multivariate analysis complete: ${ANALYSIS_ROOT}"