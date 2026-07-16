#!/bin/bash
#SBATCH --job-name=Morphometrics
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH -o logs/stage3b_%j.out
#SBATCH -e logs/stage3b_%j.err

set -euo pipefail

echo "=== Stage 3b Morphometrics ==="
echo "Timestamp      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Hostname       : $(hostname)"
echo "Working Dir    : $(pwd)"
echo "Pipeline Root  : ${PIPELINE_ROOT:-<unset>}"
echo "SWC Root       : ${SWC_ROOT:-<unset>}"
echo "Morpho Root    : ${MORPHO_ROOT:-<unset>}"

[[ -d "${SWC_ROOT:-}" ]] || { echo "ERROR: SWC_ROOT does not exist" >&2; exit 1; }
[[ -n "${MORPHO_PYTHON:-}" ]] || { echo "ERROR: MORPHO_PYTHON not set" >&2; exit 1; }

find "${SWC_ROOT}" -type f -name 'clean_skeleton.swc' -print -quit | grep -q . \
    || { echo "ERROR: No clean_skeleton.swc files found under ${SWC_ROOT}" >&2; exit 1; }

mkdir -p logs "${MORPHO_ROOT}"

# === Clean ParaView environment pollution ===
unset PYTHONPATH
source /scratch/flore0a/iops/miniconda3/etc/profile.d/conda.sh
conda activate swc_env

# Verify conda Python works
"${MORPHO_PYTHON}" --version || { echo "ERROR: MORPHO_PYTHON failed: ${MORPHO_PYTHON}" >&2; exit 1; }

cmd=("${MORPHO_PYTHON}" "${PIPELINE_ROOT}/Postanalysis_swc/compute_morphometrics.py" \
     --path "${SWC_ROOT}" \
     --out "${MORPHO_ROOT}/morphometrics.csv")

printf 'Command        :'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

echo "[success] Morphometrics complete: ${MORPHO_ROOT}/morphometrics.csv"