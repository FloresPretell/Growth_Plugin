#!/bin/bash
#SBATCH --job-name=Morphometrics_swc
#SBATCH --partition=workq
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail
PYTHON_BIN="/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python"

mkdir -p logs

BASE_PATH=/scratch/flore0a/Dataset_test2/Cases

echo "Computing morphometrics from ${BASE_PATH}"

time "${PYTHON_BIN}" ./compute_morphometrics.py \
    --path "${BASE_PATH}" \
    --out morphometrics.csv

echo "Done. CSV written to morphometrics.csv"
