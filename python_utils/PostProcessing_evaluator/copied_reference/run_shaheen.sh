#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG="${SCRIPT_DIR}/config_shaheen.sh"

CONFIG_FILE="${DEFAULT_CONFIG}"
DRY_RUN=0
MAX_CASES_OVERRIDE=""

usage() {
    cat <<'EOF'
Usage:
  ./run_shaheen.sh [--config FILE] [--dry-run] [--max-cases N]

Submits the full post-processing pipeline to Slurm on Shaheen.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --max-cases)
            MAX_CASES_OVERRIDE="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

CONFIG_FILE="$(cd "$(dirname "${CONFIG_FILE}")" && pwd)/$(basename "${CONFIG_FILE}")"
[[ -f "${CONFIG_FILE}" ]] || { echo "ERROR: Config file not found: ${CONFIG_FILE}" >&2; exit 1; }

# shellcheck source=/dev/null
source "${CONFIG_FILE}"

require_var() {
    local name="$1"
    [[ -n "${!name:-}" ]] || { echo "ERROR: Required config variable ${name} is empty" >&2; exit 1; }
}

command_preview() {
    local quoted=()
    local arg
    for arg in "$@"; do
        quoted+=("$(printf '%q' "${arg}")")
    done
    printf '%s\n' "${quoted[*]}"
}

log_msg() {
    printf '%s\n' "$*" | tee -a "${RUN_LOG}"
}



log_context() {
    log_msg "Timestamp      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
    log_msg "Hostname       : $(hostname)"
    log_msg "Working Dir    : $(pwd)"
    log_msg "Script Dir     : ${SCRIPT_DIR}"
    log_msg "Config File    : ${CONFIG_FILE}"
    log_msg "Dataset Root   : ${DATASET_ROOT}"
    log_msg "Output Root    : ${OUTPUT_ROOT}"
    log_msg "Results Root   : ${RESULTS_ROOT}"
    log_msg "Max Cases      : ${MAX_CASES:-<all>}"
}

require_var DATASET_ROOT
require_var OUTPUT_ROOT
require_var SLURM_ACCOUNT
require_var SLURM_PARTITION
require_var SLURM_MAIL

SBATCH_DRY_RUN_COUNTER=0

LOCAL_PYTHON="${LOCAL_PYTHON:-python3}"
PVPYTHON_CMD="${PVPYTHON_CMD:-pvbatch}"
PVPYTHON_ARGS="${PVPYTHON_ARGS:---force-offscreen-rendering}"
PVPYTHON_SOURCE="${PVPYTHON_SOURCE:-}"

MATRIX_PYTHON="${MATRIX_PYTHON:-${SWC_PYTHON}}"
SWC_PYTHON="${SWC_PYTHON:-${SWC_PYTHON}}"
MORPHO_PYTHON="${MORPHO_PYTHON:-${SWC_PYTHON}}"

if [[ -n "${MAX_CASES_OVERRIDE}" ]]; then
    MAX_CASES="${MAX_CASES_OVERRIDE}"
else
    MAX_CASES="${MAX_CASES:-}"
fi

[[ -d "${DATASET_ROOT}" ]] || { echo "ERROR: Dataset root does not exist: ${DATASET_ROOT}" >&2; exit 1; }
DATASET_ROOT="$(cd "${DATASET_ROOT}" && pwd)"
OUTPUT_ROOT="${OUTPUT_ROOT%/}"
mkdir -p "${OUTPUT_ROOT}"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
RESULTS_ROOT="${OUTPUT_ROOT}/PostProcessing_${RUN_ID}"
LOGS_DIR="${RESULTS_ROOT}/logs"
DIRS_DIR="${RESULTS_ROOT}/dirs"
PNG_ROOT="${RESULTS_ROOT}/PNG"
MATRIX_ROOT="${RESULTS_ROOT}/MatrixVisualization"
SWC_ROOT="${RESULTS_ROOT}/SWC"
MORPHO_ROOT="${RESULTS_ROOT}/Morphometrics"
ANALYSIS_ROOT="${RESULTS_ROOT}/MultivariateAnalysis"
REPORTS_DIR="${RESULTS_ROOT}/reports"
DIRS_PNG_FILE="${DIRS_DIR}/dirs_png.txt"
DIRS_SWC_FILE="${DIRS_DIR}/dirs_swc.txt"
CONFIG_SNAPSHOT="${RESULTS_ROOT}/config_snapshot.sh"
RUN_LOG="${LOGS_DIR}/run_shaheen.log"

mkdir -p "${LOGS_DIR}" "${DIRS_DIR}" "${PNG_ROOT}" "${MATRIX_ROOT}" "${SWC_ROOT}" "${MORPHO_ROOT}" "${ANALYSIS_ROOT}" "${REPORTS_DIR}"
cp "${CONFIG_FILE}" "${CONFIG_SNAPSHOT}"
{
    echo
    echo "RUN_ID=\"${RUN_ID}\""
    echo "RESULTS_ROOT=\"${RESULTS_ROOT}\""
    echo "MAX_CASES=\"${MAX_CASES}\""
} >> "${CONFIG_SNAPSHOT}"

: > "${RUN_LOG}"

log_msg "=== PostProcessing Shaheen run ==="
log_context
log_msg ""

if [[ "${DRY_RUN}" -eq 0 ]] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch was not found in PATH" >&2
    exit 1
fi

if ! command -v "${LOCAL_PYTHON}" >/dev/null 2>&1; then
    echo "ERROR: Manifest python executable not found in PATH: ${LOCAL_PYTHON}" >&2
    exit 1
fi

log_msg "--- Stage 1 manifest: dataset cases ---"
"${LOCAL_PYTHON}" "${SCRIPT_DIR}/Postanalysis/generate_folder_list_vtkhdf.py" \
    --path "${DATASET_ROOT}" \
    --output "${DIRS_PNG_FILE}" \
    --max-cases "${MAX_CASES}" 2>&1 | tee -a "${RUN_LOG}"

[[ -f "${DIRS_PNG_FILE}" ]] || { echo "ERROR: Stage 1 manifest was not created: ${DIRS_PNG_FILE}" >&2; exit 1; }
N1=$(grep -cve '^[[:space:]]*$' "${DIRS_PNG_FILE}" || true)
[[ "${N1}" -gt 0 ]] || { echo "ERROR: No cases detected in ${DIRS_PNG_FILE}" >&2; exit 1; }
log_msg "Detected cases : ${N1}"
log_msg ""

log_msg "--- Stage 3 manifest: expected PNG to SWC cases ---"
"${LOCAL_PYTHON}" "${SCRIPT_DIR}/Postanalysis_swc/generate_folder_list_swc_vtkhdf.py" \
    --source-dirs-file "${DIRS_PNG_FILE}" \
    --png-root "${PNG_ROOT}" \
    --output "${DIRS_SWC_FILE}" 2>&1 | tee -a "${RUN_LOG}"

[[ -f "${DIRS_SWC_FILE}" ]] || { echo "ERROR: SWC manifest was not created: ${DIRS_SWC_FILE}" >&2; exit 1; }
N3=$(grep -cve '^[[:space:]]*$' "${DIRS_SWC_FILE}" || true)
[[ "${N3}" -gt 0 ]] || { echo "ERROR: No SWC cases detected in ${DIRS_SWC_FILE}" >&2; exit 1; }
log_msg "Detected SWC cases : ${N3}"
log_msg ""

SBATCH_BASE=(
    --account="${SLURM_ACCOUNT}"
    --partition="${SLURM_PARTITION}"
    --mail-user="${SLURM_MAIL}"
    --mail-type=ALL
    --chdir="${RESULTS_ROOT}"
)

############################################################################################
# Stage 1: PNG Export
log_msg "--- Stage 1: PNG Export (${N1} cases) ---"
cmd1=(
    sbatch --parsable
    "${SBATCH_BASE[@]}"
    --array=0-$((N1-1))
    --export=ALL,PIPELINE_ROOT="${SCRIPT_DIR}",RESULTS_ROOT="${RESULTS_ROOT}",LOGS_DIR="${LOGS_DIR}",DIRS_FILE="${DIRS_PNG_FILE}",PNG_ROOT="${PNG_ROOT}",PVPYTHON_CMD="${PVPYTHON_CMD}",PVPYTHON_ARGS="${PVPYTHON_ARGS}",PVPYTHON_SOURCE="${PVPYTHON_SOURCE}"
    "${SCRIPT_DIR}/Postanalysis/Submite_export.slurm"
)
log_msg "Command: $(command_preview "${cmd1[@]}")"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    JOB1="DRYRUN1"
else
    JOB1=$("${cmd1[@]}")
fi
log_msg "Stage 1 job: ${JOB1}"
log_msg ""

# Stage 2: Matrix GIF
log_msg "--- Stage 2: Matrix GIF ---"
cmd2=(
    sbatch --parsable
    "${SBATCH_BASE[@]}"
    --dependency=afterok:${JOB1}
    --export=ALL,PIPELINE_ROOT="${SCRIPT_DIR}",RESULTS_ROOT="${RESULTS_ROOT}",LOGS_DIR="${LOGS_DIR}",PNG_ROOT="${PNG_ROOT}",MATRIX_ROOT="${MATRIX_ROOT}",MATRIX_PYTHON="${MATRIX_PYTHON}"
    "${SCRIPT_DIR}/slurm/stage2_matrix.slurm"
)
log_msg "Command: $(command_preview "${cmd2[@]}")"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    JOB2="DRYRUN2"
else
    JOB2=$("${cmd2[@]}")
fi
log_msg "Stage 2 job: ${JOB2} (after ${JOB1})"
log_msg ""

# Stage 3a: SWC Export
log_msg "--- Stage 3a: SWC Export (${N3} cases) ---"
cmd3a=(
    sbatch --parsable
    "${SBATCH_BASE[@]}"
    --dependency=afterok:${JOB1}
    --array=0-$((N3-1))
    --export=ALL,PIPELINE_ROOT="${SCRIPT_DIR}",RESULTS_ROOT="${RESULTS_ROOT}",LOGS_DIR="${LOGS_DIR}",DIRS_FILE="${DIRS_SWC_FILE}",SWC_ROOT="${SWC_ROOT}",SWC_PYTHON="${SWC_PYTHON}"
    "${SCRIPT_DIR}/Postanalysis_swc/Submit_export_swc.slurm"
)
log_msg "Command: $(command_preview "${cmd3a[@]}")"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    JOB3A="DRYRUN3A"
else
    JOB3A=$("${cmd3a[@]}")
fi
log_msg "Stage 3a job: ${JOB3A} (after ${JOB1})"
log_msg ""

# Stage 3b: Morphometrics
log_msg "--- Stage 3b: Morphometrics ---"
cmd3b=(
    sbatch --parsable
    "${SBATCH_BASE[@]}"
    --dependency=afterok:${JOB3A}
    --export=ALL,PIPELINE_ROOT="${SCRIPT_DIR}",RESULTS_ROOT="${RESULTS_ROOT}",LOGS_DIR="${LOGS_DIR}",SWC_ROOT="${SWC_ROOT}",MORPHO_ROOT="${MORPHO_ROOT}",MORPHO_PYTHON="${MORPHO_PYTHON}"
    "${SCRIPT_DIR}/slurm/stage3_morphometrics.slurm"
)
log_msg "Command: $(command_preview "${cmd3b[@]}")"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    JOB3B="DRYRUN3B"
else
    JOB3B=$("${cmd3b[@]}")
fi
log_msg "Stage 3b job: ${JOB3B} (after ${JOB3A})"
log_msg ""

# Stage 4: Analysis
log_msg "--- Stage 4: Analysis ---"
cmd4=(
    sbatch --parsable
    "${SBATCH_BASE[@]}"
    --dependency=afterok:${JOB3B}
    --export=ALL,PIPELINE_ROOT="${SCRIPT_DIR}",RESULTS_ROOT="${RESULTS_ROOT}",LOGS_DIR="${LOGS_DIR}",MORPHO_ROOT="${MORPHO_ROOT}",ANALYSIS_ROOT="${ANALYSIS_ROOT}",MORPHO_PYTHON="${MORPHO_PYTHON}"
    "${SCRIPT_DIR}/slurm/stage4_analysis.slurm"
)
log_msg "Command: $(command_preview "${cmd4[@]}")"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    JOB4="DRYRUN4"
else
    JOB4=$("${cmd4[@]}")
fi
log_msg "Stage 4 job: ${JOB4} (after ${JOB3B})"
log_msg ""

log_msg "=== Submission Summary ==="
log_msg "Stage 1 PNG export      : ${JOB1}"
log_msg "Stage 2 Matrix GIF      : ${JOB2} (depends on ${JOB1})"
log_msg "Stage 3a SWC export     : ${JOB3A} (depends on ${JOB1})"
log_msg "Stage 3b Morphometrics  : ${JOB3B} (depends on ${JOB3A})"
log_msg "Stage 4 Analysis        : ${JOB4} (depends on ${JOB3B})"
log_msg "Results Root            : ${RESULTS_ROOT}"