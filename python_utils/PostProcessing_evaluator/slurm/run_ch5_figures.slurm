#!/bin/bash
#SBATCH --job-name=ppe_ch5
#SBATCH --partition=workq
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH -o /scratch/flore0a/PostProcessing_evaluator/logs/ppe_ch5_%j.out
#SBATCH -e /scratch/flore0a/PostProcessing_evaluator/logs/ppe_ch5_%j.err
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL

# Lightweight matplotlib-only job (no pvpython, no MPI). Runs the 4 CH5
# modules on one enriched CSV. Inputs via environment:
#   CSV     (default: D7p5_V3 clean enriched)
#   VARIANT (default: clean)
#
# Submission:
#   CSV=... sbatch /project/.../PostProcessing_evaluator/slurm/run_ch5_figures.slurm

set -euo pipefail
mkdir -p /scratch/flore0a/PostProcessing_evaluator/{ch5_figures,logs}

CSV="${CSV:-/scratch/flore0a/PostProcessing_evaluator/reproduction_runs/20260513_010951/D7p5_V3_fields_enriched.csv}"
VARIANT="${VARIANT:-clean}"
EVAL_ROOT="/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator"
SWC_PY="${SWC_PY:-/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python}"

echo "[$(date '+%F %T')] CH5 SLURM start; CSV=${CSV}  VARIANT=${VARIANT}"
"${SWC_PY}" "${EVAL_ROOT}/scripts/ch5/main_ch5.py" \
    --csv "${CSV}" --variant "${VARIANT}"
echo "[$(date '+%F %T')] CH5 SLURM done."
