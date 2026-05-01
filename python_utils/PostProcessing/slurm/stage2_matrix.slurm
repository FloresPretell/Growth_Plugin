#!/bin/bash
#SBATCH --job-name=MatrixViz
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH -o logs/stage2_%j.out
#SBATCH -e logs/stage2_%j.err

set -euo pipefail

echo "=== Stage 2 Matrix visualization ==="
echo "Timestamp      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Hostname       : $(hostname)"
echo "Pipeline Root  : ${PIPELINE_ROOT:-<unset>}"
echo "PNG Root       : ${PNG_ROOT:-<unset>}"
echo "Matrix Root    : ${MATRIX_ROOT:-<unset>}"

[[ -d "${PNG_ROOT:-}" ]] || { echo "ERROR: PNG_ROOT does not exist" >&2; exit 1; }
find "${PNG_ROOT}" -type f -name '*.png' -print -quit | grep -q . \
    || { echo "ERROR: No PNG files found under ${PNG_ROOT}" >&2; exit 1; }

[[ -n "${MATRIX_PYTHON:-}" ]] || { echo "ERROR: MATRIX_PYTHON not set" >&2; exit 1; }

mkdir -p logs "${MATRIX_ROOT}"

# === Clean ParaView environment pollution ===
unset PYTHONPATH
source /scratch/flore0a/iops/miniconda3/etc/profile.d/conda.sh
conda activate swc_env

# Verify conda Python works
"${MATRIX_PYTHON}" --version || { echo "ERROR: MATRIX_PYTHON failed: ${MATRIX_PYTHON}" >&2; exit 1; }

cmd=("${MATRIX_PYTHON}" "${PIPELINE_ROOT}/Matrix_viosualization_png/matrix.py" \
     "${PNG_ROOT}" \
     --output "${MATRIX_ROOT}/mosaic_matrix.gif" \
     --png_subdir "u_t" \
     --duration 250 \
     --keep_last_frame)

printf 'Command        :'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

echo "[success] Matrix visualization complete"