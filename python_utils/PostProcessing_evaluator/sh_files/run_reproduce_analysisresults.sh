#!/bin/bash
# run_reproduce_analysisresults.sh
# Login-node orchestrator: re-run the 4-step SWC-sampling chain that produced
# /scratch/flore0a/AnalysisResults, into a fresh timestamped directory under
# /scratch/flore0a/PostProcessing_evaluator/reproduction_runs/.
#
# Outputs:
#   /scratch/flore0a/PostProcessing_evaluator/reproduction_runs/<TS>/
#       <case>_fields.csv, <case>_fields_enriched.csv (clean)
#       <case>_pixel_fields.csv, <case>_pixel_fields_enriched.csv (pixel)
#       figures/<case>/fig{1..4}_*.{png,svg}
#       figures/<case>/front_aware_analysis/{figA..figF,summary_report.txt}
#       figures/<case>/front_aware_pixel/{figA..figF,summary_report.txt}
#       logs/<step>.{out,err}
#       execution_manifest.json
#
# Note: this is for the *login node*. On Shaheen III, pvpython requires
# srun + PALS context (see /scratch/flore0a/Instruction_every_sessionread.md
# §17.5). The pvpython sampling step (especially --skeleton-variant pixel)
# is heavy and should normally go through the SLURM variant in ../slurm/.
# This script will still work on the login node for the clean-skeleton
# variant, which is fast (984 nodes total).

set -euo pipefail

EVAL_ROOT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator"

CASE_ID="${CASE_ID:-D7p5_V3}"
SIM_FILE="${SIM_FILE:-/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/Simulation.vtkhdf}"
SWC_DIR="${SWC_DIR:-/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc}"
VARIANT="${VARIANT:-clean}"   # clean | pixel | both
PVPYTHON_BIN="${PVPYTHON_BIN:-pvpython}"
SWC_PYTHON="${SWC_PYTHON:-/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python}"

# Pre-flight check
[[ -f "$SIM_FILE" ]] || { echo "ERROR: SIM_FILE not found: $SIM_FILE" >&2; exit 1; }
[[ -d "$SWC_DIR" ]]  || { echo "ERROR: SWC_DIR not found: $SWC_DIR"  >&2; exit 1; }

# Use the orchestrator python (no vtk needed; subprocess wrapper only)
python3 "${EVAL_ROOT}/scripts/reproduction/run_chain.py" \
    --case-id "${CASE_ID}" \
    --sim-file "${SIM_FILE}" \
    --swc-dir  "${SWC_DIR}" \
    --skeleton-variant "${VARIANT}" \
    --pvpython "${PVPYTHON_BIN}" \
    --swc-python "${SWC_PYTHON}"
