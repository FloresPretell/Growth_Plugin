#!/usr/bin/env bash
# config_dataset1_shaheen.sh
# PostProcessing pipeline config for Dataset_test1
# D = 3,6,9,12,15,18,21,24,27,30  x  V = 0.5,1,1.5,2,2.5,3,3.5,4,4.5,5
# 99 cases (vtkhdf format; D12_V0p5 has merged/ only — pipeline skips it gracefully)
# Created: 2026-06-29

DATASET_ROOT="/scratch/flore0a/Dataset_test1/Cases_dataset1"
OUTPUT_ROOT="/scratch/flore0a/Dataset_test1_Results"
MAX_CASES=""

LOCAL_PYTHON="python3"
PVPYTHON_CMD="pvbatch"
PVPYTHON_ARGS="--force-offscreen-rendering"
PNG_EXPORT_MODE="render"
PVPYTHON_SOURCE="/scratch/flore0a/Modules.sh"

# === Use conda environment for all non-ParaView tasks ===
CONDA_PYTHON="/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python"

MATRIX_PYTHON="${CONDA_PYTHON}"      # ✅ Use conda
SWC_PYTHON="${CONDA_PYTHON}"         # ✅ Use conda
MORPHO_PYTHON="${CONDA_PYTHON}"      # ✅ Use conda

SLURM_ACCOUNT="k10070"
SLURM_PARTITION="workq"
SLURM_MAIL="elsanicole.florespretell@kaust.edu.sa"
