#!/bin/bash
#SBATCH --job-name=ppe_repro_analysisresults
#SBATCH --partition=workq
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH -o /scratch/flore0a/PostProcessing_evaluator/logs/ppe_repro_%j.out
#SBATCH -e /scratch/flore0a/PostProcessing_evaluator/logs/ppe_repro_%j.err
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL

# SLURM job: reproduce the 4-step SWC-sampling chain that produced
# /scratch/flore0a/AnalysisResults, into a fresh timestamped directory under
# /scratch/flore0a/PostProcessing_evaluator/reproduction_runs/.
#
# Inputs via environment:
#   CASE_ID     (default: D7p5_V3)
#   SIM_FILE    (default: D7.5_V3 VTKHDF)
#   SWC_DIR     (default: D7.5_V3 merged/swc)
#   VARIANT     (default: both)   clean | pixel | both
#
# Per /scratch/flore0a/Instruction_every_sessionread.md:
#  - §17.5: pvpython requires srun --ntasks=1 on Shaheen III.
#  - §16:   everything written from the compute node must go to /scratch.
#  - §10:   --mail-user and --mail-type=ALL are required.
#
# Submit from the login node:
#   sbatch /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator/slurm/run_reproduce_analysisresults.slurm

set -euo pipefail

EVAL_ROOT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator"
mkdir -p /scratch/flore0a/PostProcessing_evaluator/{reproduction_runs,logs}

CASE_ID="${CASE_ID:-D7p5_V3}"
SIM_FILE="${SIM_FILE:-/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/Simulation.vtkhdf}"
SWC_DIR="${SWC_DIR:-/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc}"
VARIANT="${VARIANT:-both}"
SWC_PYTHON="${SWC_PYTHON:-/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python}"

set +e
source /scratch/flore0a/Modules.sh 2>&1
set -e

PVPYTHON="$(which pvpython 2>/dev/null || true)"
if [[ -z "${PVPYTHON}" ]]; then
    echo "ERROR: pvpython not found after sourcing Modules.sh" >&2
    exit 1
fi

# Modules.sh sets PYTHONPATH to paraview's Python 3.12 stdlib. That poisons
# swc_env (Python 3.10) and any non-pvpython child process — see the
# `_sre.MAGIC` mismatch traceback in job 12016460. pvpython has its own
# bootstrapped path resolution and does NOT need PYTHONPATH (verified
# empirically: `pvpython -c "import vtk"` works with PYTHONPATH unset).
unset PYTHONPATH

# Determine which python to use for the orchestrator (no vtk needed).
ORCH_PY="${ORCH_PY:-${SWC_PYTHON}}"
if ! command -v "${ORCH_PY}" >/dev/null 2>&1; then
    echo "WARNING: ${ORCH_PY} not found; falling back to system python3" >&2
    ORCH_PY="python3"
fi

# Pre-flight
[[ -f "$SIM_FILE" ]] || { echo "ERROR: SIM_FILE missing: $SIM_FILE" >&2; exit 1; }
[[ -d "$SWC_DIR" ]]  || { echo "ERROR: SWC_DIR missing: $SWC_DIR"  >&2; exit 1; }

echo "[$(date '+%F %T')] starting reproduction"
echo "  evaluator     : ${EVAL_ROOT}"
echo "  case_id       : ${CASE_ID}"
echo "  sim_file      : ${SIM_FILE}"
echo "  swc_dir       : ${SWC_DIR}"
echo "  variant       : ${VARIANT}"
echo "  pvpython      : ${PVPYTHON}"
echo "  swc_python    : ${SWC_PYTHON}"
echo "  orchestrator  : ${ORCH_PY}"

# pvpython call inside run_chain.py is wrapped through Popen, so the
# inner pvpython invocation also needs to go through srun. The orchestrator
# uses PVPYTHON as a literal argv[0] — wrap it with a shim that does srun.
SHIM_DIR="${TMPDIR:-/scratch/flore0a/PostProcessing_evaluator/temporary}/srun_pvpython_$$"
mkdir -p "$SHIM_DIR"
cat > "$SHIM_DIR/pvpython" <<EOF
#!/bin/bash
exec srun --ntasks=1 "${PVPYTHON}" --force-offscreen-rendering "\$@"
EOF
chmod +x "$SHIM_DIR/pvpython"
export PATH="$SHIM_DIR:$PATH"
export PVPYTHON="$SHIM_DIR/pvpython"

"${ORCH_PY}" "${EVAL_ROOT}/scripts/reproduction/run_chain.py" \
    --case-id "${CASE_ID}" \
    --sim-file "${SIM_FILE}" \
    --swc-dir  "${SWC_DIR}" \
    --skeleton-variant "${VARIANT}" \
    --pvpython "${PVPYTHON}" \
    --swc-python "${SWC_PYTHON}"

rc=$?
echo "[$(date '+%F %T')] done. rc=$rc"
exit $rc
