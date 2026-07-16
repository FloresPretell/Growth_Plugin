#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG="${SCRIPT_DIR}/config_local.sh"

CONFIG_FILE="${DEFAULT_CONFIG}"
DRY_RUN=0
MAX_CASES_OVERRIDE=""
MAX_TIMESTEPS_OVERRIDE=""

usage() {
    cat <<'EOF'
Usage:
  ./run_local.sh [--config FILE] [--dry-run] [--max-cases N] [--max-timesteps N]

Runs the full post-processing pipeline locally without Slurm.
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
        --max-timesteps)
            MAX_TIMESTEPS_OVERRIDE="$2"
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

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    local label="$1"
    local cmd="$2"
    if command -v "${cmd}" >/dev/null 2>&1; then
        return 0
    fi
    fail "${label} executable not found in PATH: ${cmd}"
}

require_python_module() {
    local python_cmd="$1"
    local module_name="$2"
    local hint="$3"
    if ! "${python_cmd}" -c "import ${module_name}" >/dev/null 2>&1; then
        fail "${module_name} is not available in ${python_cmd}. ${hint}"
    fi
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
run_or_print() {
    local cmd=("$@")
    log_msg "Command: $(command_preview "${cmd[@]}")"
    if [[ "${DRY_RUN}" -eq 0 ]]; then
        "${cmd[@]}" 2>&1 | tee -a "${RUN_LOG}"
    fi
}

run_setup_command() {
    local cmd=("$@")
    log_msg "Command: $(command_preview "${cmd[@]}")"
    "${cmd[@]}" 2>&1 | tee -a "${RUN_LOG}"
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
    log_msg "PNG Mode       : ${PNG_EXPORT_MODE}"
    log_msg "Local Fields   : ${LOCAL_FIELDS_TO_EXPORT:-u_u u_p u_t u_b u_ca_cyt}"
    log_msg "Max Timesteps  : ${LOCAL_MAX_TIMESTEPS:-<all>}"
}

require_var DATASET_ROOT
require_var OUTPUT_ROOT

LOCAL_CONDA_ENV="${LOCAL_CONDA_ENV:-}"
LOCAL_PYTHON="${LOCAL_PYTHON:-python}"
PNG_EXPORT_PYTHON="${PNG_EXPORT_PYTHON:-${LOCAL_PYTHON}}"
PNG_EXPORT_MODE="${PNG_EXPORT_MODE:-local_safe}"
MATRIX_PYTHON="${MATRIX_PYTHON:-${LOCAL_PYTHON}}"
SWC_PYTHON="${SWC_PYTHON:-${LOCAL_PYTHON}}"
MORPHO_PYTHON="${MORPHO_PYTHON:-${LOCAL_PYTHON}}"

if [[ -n "${MAX_CASES_OVERRIDE}" ]]; then
    MAX_CASES="${MAX_CASES_OVERRIDE}"
else
    MAX_CASES="${MAX_CASES:-}"
fi

if [[ -n "${MAX_TIMESTEPS_OVERRIDE}" ]]; then
    LOCAL_MAX_TIMESTEPS="${MAX_TIMESTEPS_OVERRIDE}"
else
    LOCAL_MAX_TIMESTEPS="${LOCAL_MAX_TIMESTEPS:-}"
fi

if [[ -n "${LOCAL_MAX_TIMESTEPS}" ]] && ! [[ "${LOCAL_MAX_TIMESTEPS}" =~ ^[0-9]+$ ]] ; then
    fail "LOCAL_MAX_TIMESTEPS must be a positive integer, got: ${LOCAL_MAX_TIMESTEPS}"
fi

DEFAULT_LOCAL_FIELDS=(u_u u_p u_t u_b u_ca_cyt)
if [[ -n "${LOCAL_FIELDS_TO_EXPORT:-}" ]]; then
    read -r -a LOCAL_FIELDS <<< "${LOCAL_FIELDS_TO_EXPORT}"
else
    LOCAL_FIELDS=("${DEFAULT_LOCAL_FIELDS[@]}")
fi
[[ "${#LOCAL_FIELDS[@]}" -gt 0 ]] || fail "LOCAL_FIELDS_TO_EXPORT resolved to an empty field list"

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
RUN_LOG="${LOGS_DIR}/run_local.log"

mkdir -p "${LOGS_DIR}" "${DIRS_DIR}" "${PNG_ROOT}" "${MATRIX_ROOT}" "${SWC_ROOT}" "${MORPHO_ROOT}" "${ANALYSIS_ROOT}" "${REPORTS_DIR}"
export MPLCONFIGDIR="${RESULTS_ROOT}/.matplotlib"
mkdir -p "${MPLCONFIGDIR}"
export PYTHONUNBUFFERED=1
cp "${CONFIG_FILE}" "${CONFIG_SNAPSHOT}"
{
    echo
    echo "RUN_ID=\"${RUN_ID}\""
    echo "RESULTS_ROOT=\"${RESULTS_ROOT}\""
    echo "MAX_CASES=\"${MAX_CASES}\""
    echo "LOCAL_MAX_TIMESTEPS=\"${LOCAL_MAX_TIMESTEPS}\""
} >> "${CONFIG_SNAPSHOT}"

: > "${RUN_LOG}"

log_msg "=== PostProcessing local run ==="
log_context
log_msg ""

if [[ -n "${LOCAL_CONDA_ENV}" ]]; then
    log_msg "Activating Conda environment: ${LOCAL_CONDA_ENV}"
    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
        conda activate "${LOCAL_CONDA_ENV}"
    else
        echo "ERROR: conda not found, but LOCAL_CONDA_ENV is set to ${LOCAL_CONDA_ENV}" >&2
        exit 1
    fi
fi

require_command "local python" "${LOCAL_PYTHON}"
require_command "PNG export python" "${PNG_EXPORT_PYTHON}"
require_command "matrix python" "${MATRIX_PYTHON}"
require_command "SWC python" "${SWC_PYTHON}"
require_command "morphometrics python" "${MORPHO_PYTHON}"
require_python_module "${PNG_EXPORT_PYTHON}" "paraview.simple" "Activate/install ParaView in the selected local Python environment or update PNG_EXPORT_PYTHON in config_local.sh."
if [[ "${PNG_EXPORT_MODE}" == "local_safe" ]]; then
    require_python_module "${PNG_EXPORT_PYTHON}" "matplotlib" "Install matplotlib in the selected local Python environment for local-safe PNG export."
    require_python_module "${PNG_EXPORT_PYTHON}" "scipy" "Install scipy in the selected local Python environment for local-safe PNG export."
elif [[ "${PNG_EXPORT_MODE}" != "render" ]]; then
    fail "Unsupported PNG_EXPORT_MODE: ${PNG_EXPORT_MODE}"
fi

log_msg "--- Stage 1: dataset manifest ---"
run_setup_command \
    "${LOCAL_PYTHON}" "${SCRIPT_DIR}/Postanalysis/generate_folder_list.py" \
    --path "${DATASET_ROOT}" \
    --output "${DIRS_PNG_FILE}" \
    --max-cases "${MAX_CASES}"

[[ -f "${DIRS_PNG_FILE}" ]] || { echo "ERROR: Stage 1 manifest was not created: ${DIRS_PNG_FILE}" >&2; exit 1; }
N1=$(grep -cve '^[[:space:]]*$' "${DIRS_PNG_FILE}" || true)
[[ "${N1}" -gt 0 ]] || { echo "ERROR: No cases detected in ${DIRS_PNG_FILE}" >&2; exit 1; }
log_msg "Detected cases : ${N1}"
log_msg ""

log_msg "--- Stage 1: PNG export ---"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    while IFS= read -r simdir; do
        [[ -n "${simdir}" ]] || continue
        cmd=("${PNG_EXPORT_PYTHON}" "${SCRIPT_DIR}/Postanalysis/main_export.py" --simdir "${simdir}" --outdir "${PNG_ROOT}" --export-mode "${PNG_EXPORT_MODE}" --fields_to_export "${LOCAL_FIELDS[@]}")
        if [[ -n "${LOCAL_MAX_TIMESTEPS}" ]]; then
            cmd+=(--max-timesteps "${LOCAL_MAX_TIMESTEPS}")
        fi
        log_msg "Command: $(command_preview "${cmd[@]}")"
    done < "${DIRS_PNG_FILE}"
else
    while IFS= read -r simdir; do
        [[ -n "${simdir}" ]] || continue
        [[ -d "${simdir}" ]] || { echo "ERROR: Missing simulation directory: ${simdir}" >&2; exit 1; }
        find "${simdir}" -maxdepth 1 \( -name 'Output_t*.vtu' -o -name 'Output_t*.pvtu' \) -print -quit | grep -q . \
            || { echo "ERROR: No VTU/PVTU files found in ${simdir}" >&2; exit 1; }
        cmd=("${PNG_EXPORT_PYTHON}" "${SCRIPT_DIR}/Postanalysis/main_export.py" --simdir "${simdir}" --outdir "${PNG_ROOT}" --export-mode "${PNG_EXPORT_MODE}" --fields_to_export "${LOCAL_FIELDS[@]}")
        if [[ -n "${LOCAL_MAX_TIMESTEPS}" ]]; then
            cmd+=(--max-timesteps "${LOCAL_MAX_TIMESTEPS}")
        fi
        run_or_print "${cmd[@]}"
    done < "${DIRS_PNG_FILE}"
fi
log_msg ""

log_msg "--- Stage 3 manifest: PNG to SWC cases ---"
run_setup_command \
    "${LOCAL_PYTHON}" "${SCRIPT_DIR}/Postanalysis_swc/generate_folder_list_swc.py" \
    --source-dirs-file "${DIRS_PNG_FILE}" \
    --png-root "${PNG_ROOT}" \
    --output "${DIRS_SWC_FILE}"

[[ -f "${DIRS_SWC_FILE}" ]] || { echo "ERROR: SWC manifest was not created: ${DIRS_SWC_FILE}" >&2; exit 1; }
N3=$(grep -cve '^[[:space:]]*$' "${DIRS_SWC_FILE}" || true)
[[ "${N3}" -gt 0 ]] || { echo "ERROR: No SWC cases detected in ${DIRS_SWC_FILE}" >&2; exit 1; }
log_msg "Detected SWC cases : ${N3}"
log_msg ""

log_msg "--- Stage 2: matrix visualization ---"
run_or_print \
    "${MATRIX_PYTHON}" "${SCRIPT_DIR}/Matrix_viosualization_png/matrix.py" \
    "${PNG_ROOT}" \
    --output "${MATRIX_ROOT}/mosaic_matrix.gif" \
    --png_subdir "u_t" \
    --duration 250 \
    --keep_last_frame
log_msg ""

log_msg "--- Stage 3a: SWC export ---"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    while IFS= read -r png_dir; do
        [[ -n "${png_dir}" ]] || continue
        log_msg "Command: $(command_preview "${SWC_PYTHON}" "${SCRIPT_DIR}/Postanalysis_swc/main_export_swc.py" --simdir "${png_dir}" --outdir "${SWC_ROOT}")"
    done < "${DIRS_SWC_FILE}"
else
    while IFS= read -r png_dir; do
        [[ -n "${png_dir}" ]] || continue
        [[ -d "${png_dir}" ]] || { echo "ERROR: Missing PNG directory for SWC export: ${png_dir}" >&2; exit 1; }
        find "${png_dir}" -maxdepth 1 -name '*.png' -print -quit | grep -q . \
            || { echo "ERROR: No PNG files found for SWC export in ${png_dir}" >&2; exit 1; }
        run_or_print \
            "${SWC_PYTHON}" "${SCRIPT_DIR}/Postanalysis_swc/main_export_swc.py" \
            --simdir "${png_dir}" \
            --outdir "${SWC_ROOT}"
    done < "${DIRS_SWC_FILE}"
fi
log_msg ""

log_msg "--- Stage 3b: morphometrics ---"
run_or_print \
    "${MORPHO_PYTHON}" "${SCRIPT_DIR}/Postanalysis_swc/compute_morphometrics.py" \
    --path "${SWC_ROOT}" \
    --out "${MORPHO_ROOT}/morphometrics.csv"
log_msg ""

log_msg "--- Stage 4: multivariate analysis ---"
run_or_print \
    "${MORPHO_PYTHON}" "${SCRIPT_DIR}/Morphoanalysis_plot/analyze_morphometrics.py" \
    --input "${MORPHO_ROOT}/morphometrics.csv" \
    --output-dir "${ANALYSIS_ROOT}"
log_msg ""

log_msg "=== Local pipeline complete ==="
log_msg "Results Root : ${RESULTS_ROOT}"
