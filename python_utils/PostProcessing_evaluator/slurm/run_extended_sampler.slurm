#!/bin/bash
#SBATCH --job-name=ppe_ext_sampler
#SBATCH --partition=workq
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:30:00
#SBATCH --mem=64G
#SBATCH -o /scratch/flore0a/PostProcessing_evaluator/logs/ppe_ext_sampler_%j.out
#SBATCH -e /scratch/flore0a/PostProcessing_evaluator/logs/ppe_ext_sampler_%j.err
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL

# Run scripts/reproduction/sample_fields_on_swc_extended.py under pvpython.
# Writes new CSVs to /scratch/.../reproduction_runs/<TS>/ — never overwrites
# the existing 20260513_010951 (clean) or 20260513_011220 (pixel) runs.
#
# Required env on submit:
#   CASE_ID     — short case label (e.g. D7p5_V3, D2p5_V3)
#   SIM_FILE    — Simulation.vtkhdf for that case
#   SWC_DIR     — merged/swc directory for that case
#   RUN_TS      — optional fixed timestamp (so both swc-types land in one TS)

set -euo pipefail
mkdir -p /scratch/flore0a/PostProcessing_evaluator/{reproduction_runs,logs}

CASE_ID="${CASE_ID:-D7p5_V3}"
SIM_FILE="${SIM_FILE:-/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/Simulation.vtkhdf}"
SWC_DIR="${SWC_DIR:-/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc}"
RUN_TS="${RUN_TS:-$(date +"%Y%m%d_%H%M%S")}"
RUN_DIR="/scratch/flore0a/PostProcessing_evaluator/reproduction_runs/${RUN_TS}_extended"
mkdir -p "$RUN_DIR/logs"

EVAL_ROOT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator"

set +e
source /scratch/flore0a/Modules.sh 2>&1
set -e
unset PYTHONPATH

PVPYTHON="$(which pvpython 2>/dev/null || true)"
[[ -n "$PVPYTHON" ]] || { echo "ERROR: pvpython missing" >&2; exit 1; }

echo "[$(date '+%F %T')] CASE_ID=${CASE_ID} RUN_DIR=${RUN_DIR}"
echo "  sim_file=${SIM_FILE}"
echo "  swc_dir=${SWC_DIR}"

run_one() {
    local SWC_TYPE="$1"   # clean | pixel
    local OUT_CSV="$2"
    local LOG="$3"
    echo "[$(date '+%F %T')] sampling ${SWC_TYPE} -> ${OUT_CSV}"
    srun --ntasks=1 "${PVPYTHON}" --force-offscreen-rendering \
        "${EVAL_ROOT}/scripts/reproduction/sample_fields_on_swc_extended.py" \
        --sim-file  "${SIM_FILE}" \
        --swc-dir   "${SWC_DIR}" \
        --out-csv   "${OUT_CSV}" \
        --case-id   "${CASE_ID}" \
        --swc-type  "${SWC_TYPE}" \
        >>"${LOG}" 2>&1
}

CLEAN_CSV="${RUN_DIR}/${CASE_ID}_fields.csv"
PIXEL_CSV="${RUN_DIR}/${CASE_ID}_pixel_fields.csv"

run_one clean "${CLEAN_CSV}" "${RUN_DIR}/logs/sample_clean.log"
run_one pixel "${PIXEL_CSV}" "${RUN_DIR}/logs/sample_pixel.log"

echo "[$(date '+%F %T')] done. Run dir: ${RUN_DIR}"
ls -lah "${RUN_DIR}"
