#!/bin/bash
# run_figure_design.sh
# Login-node runner for the second-generation figure-design iteration.
# Reads one enriched CSV; writes a fresh timestamped directory under
# /scratch/flore0a/PostProcessing_evaluator/figure_design_runs/<TS>/<variant>/.
#
# These are exploratory candidate figures for discussion with Nicole.
# Old reproduction outputs are not modified or overwritten.

set -euo pipefail

EVAL_ROOT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator"
SWC_PY="${SWC_PY:-/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python}"

CSV="${1:?usage: $0 <enriched CSV> [variant] [run_ts]}"
VARIANT="${2:-clean}"
RUN_TS="${3:-}"

cmd=("${SWC_PY}" "${EVAL_ROOT}/scripts/figure_design/run_figure_design.py"
     --csv "${CSV}" --variant "${VARIANT}")
if [[ -n "$RUN_TS" ]]; then
    cmd+=(--run-ts "$RUN_TS")
fi

echo "[$(date '+%F %T')] starting figure_design (${VARIANT})"
echo "  csv : ${CSV}"
"${cmd[@]}"
echo "[$(date '+%F %T')] done."
