#!/bin/bash
# run_inspect_simulation_fields.sh
# Login-node helper: run scripts/inspect_simulation_fields.py on one simulation
# directory and write the timestamped run under scratch.
#
# This script DOES NOT submit a SLURM job. For heavy long runs use the SLURM
# variant in ../slurm/. For the first-pass inspection on one case (44 timesteps,
# fully I/O-bound on the VTKHDF), the login node + plain python3 + vtk is fine.
#
# Default input below points at D7.5_V3 (the AnalysisResults reference case).
# Override on the command line:
#     bash run_inspect_simulation_fields.sh /path/to/Simulation2D

set -euo pipefail

EVAL_ROOT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator"
SCRATCH_ROOT="/scratch/flore0a/PostProcessing_evaluator"

INPUT_DIR="${1:-/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D}"
CASE_ID="${2:-D7p5_V3}"

# Prefer pvpython on Shaheen (vtk is always present in ParaView).
# A plain python3 with vtk also works.
PY_BIN="${PY_BIN:-python3}"

if ! command -v "${PY_BIN}" >/dev/null 2>&1; then
    echo "ERROR: ${PY_BIN} not in PATH. Set PY_BIN, or run the SLURM variant under pvpython." >&2
    exit 1
fi

mkdir -p "${SCRATCH_ROOT}/runs" "${SCRATCH_ROOT}/logs"

echo "[$(date '+%F %T')] inspect_simulation_fields starting"
echo "  python    : $(command -v ${PY_BIN})"
echo "  input     : ${INPUT_DIR}"
echo "  case_id   : ${CASE_ID}"
echo "  scratch   : ${SCRATCH_ROOT}/runs/<timestamp>/"

"${PY_BIN}" "${EVAL_ROOT}/scripts/inspect_simulation_fields.py" \
    --input "${INPUT_DIR}" \
    --output-root "${SCRATCH_ROOT}/runs" \
    --case-id "${CASE_ID}"

echo "[$(date '+%F %T')] done."
