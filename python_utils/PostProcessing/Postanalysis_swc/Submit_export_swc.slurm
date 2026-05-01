#!/bin/bash
#SBATCH --job-name=Export_swc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail

echo "=== Stage 3a SWC export ==="
echo "Timestamp      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Hostname       : $(hostname)"
echo "Pipeline Root  : ${PIPELINE_ROOT:-<unset>}"
echo "SWC Root       : ${SWC_ROOT:-<unset>}"
echo "Dirs File      : ${DIRS_FILE:-<unset>}"

[[ -n "${DIRS_FILE:-}" && -f "${DIRS_FILE}" ]] || { echo "ERROR: DIRS_FILE missing" >&2; exit 1; }
[[ -n "${SWC_ROOT:-}" ]] || { echo "ERROR: SWC_ROOT not set" >&2; exit 1; }
[[ -n "${SWC_PYTHON:-}" ]] || { echo "ERROR: SWC_PYTHON not set" >&2; exit 1; }

mkdir -p logs "${SWC_ROOT}"

# === CRITICAL: Clean ParaView environment pollution ===
# These variables were set by ParaView modules and break conda Python
unset PYTHONPATH
source /scratch/flore0a/iops/miniconda3/etc/profile.d/conda.sh
conda activate swc_env
# Don't fully unset LD_LIBRARY_PATH (system needs it), but we'll let conda handle it

# Test conda Python
"${SWC_PYTHON}" --version || { echo "ERROR: SWC_PYTHON failed: ${SWC_PYTHON}" >&2; exit 1; }

SIMDIR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${DIRS_FILE}" | xargs)
[[ -n "${SIMDIR}" ]] || { echo "ERROR: No case found" >&2; exit 1; }
[[ -d "${SIMDIR}" ]] || { echo "ERROR: PNG directory not found: ${SIMDIR}" >&2; exit 1; }

cmd=("${SWC_PYTHON}" "${PIPELINE_ROOT}/Postanalysis_swc/main_export_swc.py" \
     --simdir "${SIMDIR}" \
     --outdir "${SWC_ROOT}")

printf 'Command        :'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

echo "[success] SWC export complete"