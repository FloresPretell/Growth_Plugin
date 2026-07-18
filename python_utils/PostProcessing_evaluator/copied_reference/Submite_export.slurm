#!/bin/bash
#SBATCH --job-name=Export_png
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=23:59:00
#SBATCH --mem=200G
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail

echo "=== Stage 1 PNG export ==="
echo "Timestamp      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "Hostname       : $(hostname)"
echo "Working Dir    : $(pwd)"
echo "Pipeline Root  : ${PIPELINE_ROOT:-<unset>}"
echo "Results Root   : ${RESULTS_ROOT:-<unset>}"
echo "PNG Root       : ${PNG_ROOT:-<unset>}"
echo "Dirs File      : ${DIRS_FILE:-<unset>}"

[[ -n "${DIRS_FILE:-}" && -f "${DIRS_FILE}" ]] || { echo "ERROR: DIRS_FILE is missing or invalid" >&2; exit 1; }
[[ -n "${PNG_ROOT:-}" ]] || { echo "ERROR: PNG_ROOT is not set" >&2; exit 1; }
mkdir -p "${LOGS_DIR:-logs}" "${PNG_ROOT}"


if [[ -n "${PVPYTHON_SOURCE:-}" ]]; then
    [[ -f "${PVPYTHON_SOURCE}" ]] || { echo "ERROR: PVPYTHON_SOURCE not found: ${PVPYTHON_SOURCE}" >&2; exit 1; }
    # shellcheck source=/dev/null
    source "${PVPYTHON_SOURCE}"
fi

SIMDIR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${DIRS_FILE}" | xargs)
[[ -n "${SIMDIR}" ]] || { echo "ERROR: No case found in ${DIRS_FILE} for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
[[ -d "${SIMDIR}" ]] || { echo "ERROR: Simulation directory does not exist: ${SIMDIR}" >&2; exit 1; }

if [[ ! -s "${SIMDIR}/Simulation.vtkhdf" ]]; then
    echo "ERROR: No non-empty Simulation.vtkhdf found in ${SIMDIR}" >&2
    exit 1
fi
cmd=("${PVPYTHON_CMD}")
if [[ -n "${PVPYTHON_ARGS:-}" ]]; then
    # shellcheck disable=SC2206
    extra_args=(${PVPYTHON_ARGS})
    cmd+=("${extra_args[@]}")
fi
cmd+=("${PIPELINE_ROOT}/Postanalysis/main_export.py" --simdir "${SIMDIR}" --outdir "${PNG_ROOT}" --fields_to_export u_u u_p u_t u_b u_ca_cyt --export-mode render )

printf 'Command        :'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"
