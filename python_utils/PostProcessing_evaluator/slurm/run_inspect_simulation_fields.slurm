#!/bin/bash
#SBATCH --job-name=ppe_inspect
#SBATCH --partition=workq
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH -o /scratch/flore0a/PostProcessing_evaluator/logs/ppe_inspect_%j.out
#SBATCH -e /scratch/flore0a/PostProcessing_evaluator/logs/ppe_inspect_%j.err
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL

# SLURM job: field inventory of one simulation directory under pvpython.
# Inputs are provided via environment variables exported by the submitter:
#   INPUT_DIR  — simulation directory (Simulation.vtkhdf or merged/*.vtu)
#   CASE_ID    — short label embedded in run_metadata.json
# Output: /scratch/flore0a/PostProcessing_evaluator/runs/<TS>/{field_inventory.{csv,md},run_metadata.json,log.txt}
#
# Per /scratch/flore0a/Instruction_every_sessionread.md §17.5, pvpython MUST be
# launched via srun --ntasks=1 on Shaheen III. Direct invocation crashes
# MPI_Init: pals_init2() failed: 2.
#
# Submit from the login node:
#   INPUT_DIR=/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D \
#   CASE_ID=D7p5_V3 \
#   sbatch /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator/slurm/run_inspect_simulation_fields.slurm

set -euo pipefail

EVAL_ROOT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator"
SCRATCH_ROOT="/scratch/flore0a/PostProcessing_evaluator"
mkdir -p "${SCRATCH_ROOT}/runs" "${SCRATCH_ROOT}/logs"

: "${INPUT_DIR:?ERROR: export INPUT_DIR=<path to Simulation2D>}"
: "${CASE_ID:=unspecified}"

PY_SCRIPT="${EVAL_ROOT}/scripts/inspect_simulation_fields.py"

set +e
source /scratch/flore0a/Modules.sh 2>&1
set -e

PVPYTHON="$(which pvpython 2>/dev/null || true)"
if [[ -z "${PVPYTHON}" ]]; then
    echo "ERROR: pvpython not found after sourcing Modules.sh" >&2
    exit 1
fi

echo "[$(date '+%F %T')] starting"
echo "  pvpython : ${PVPYTHON}"
echo "  script   : ${PY_SCRIPT}"
echo "  input    : ${INPUT_DIR}"
echo "  case_id  : ${CASE_ID}"

srun --ntasks=1 "${PVPYTHON}" --force-offscreen-rendering "${PY_SCRIPT}" \
    --input "${INPUT_DIR}" \
    --output-root "${SCRATCH_ROOT}/runs" \
    --case-id "${CASE_ID}"

echo "[$(date '+%F %T')] done."
