#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# value_extraction_pipeline.slurm
# Full SWC-field value-extraction pipeline (single case or array over cases)
#
# USAGE — single case:
#   sbatch --export=ALL,CASE_ID=D7p5_V3 value_extraction_pipeline.slurm
#
# USAGE — array over a list of case IDs in cases.txt (one ID per line):
#   sbatch --array=0-<N-1> value_extraction_pipeline.slurm
#   where N = number of lines in cases.txt
#
# PIPELINE STAGES
#   1. pvpython sample_fields_on_swc.py   — probe simulation fields at SWC nodes
#   2. python   enrich_swc_samples.py     — add morpho-dynamic columns (L0, GC, …)
#   3. python   plot_swc_field_samples.py — 4-panel overview figures
#   4. python   plot_front_aware_analysis.py — front-aware kymographs + A-F figures
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=ValueExtract
#SBATCH --partition=workq
#SBATCH -A k10070
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --mail-user=elsanicole.florespretell@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

set -euo pipefail
mkdir -p logs

# ─── USER-CONFIGURABLE PATHS ─────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES_ROOT="${CASES_ROOT:-/scratch/flore0a/Dataset_test2/Cases}"
OUT_ROOT="${OUT_ROOT:-/scratch/flore0a/AnalysisResults}"
CONDA_SH="${CONDA_SH:-/scratch/flore0a/iops/miniconda3/etc/profile.d/conda.sh}"
SWC_TYPE="${SWC_TYPE:-pixel}"          # pixel | clean
SMOOTH_SIGMA="${SMOOTH_SIGMA:-1.5}"    # Gaussian sigma for pixel-skeleton smoothing
R_GC="${R_GC:-0.1}"                   # Growth-cone radius (physical units)

# ─── CASE ID RESOLUTION ──────────────────────────────────────────────────────
# If CASE_ID is set (--export=ALL,CASE_ID=...) use it directly.
# Otherwise read from cases.txt using the array task index.
if [[ -z "${CASE_ID:-}" ]]; then
    CASES_FILE="${SCRIPT_DIR}/cases.txt"
    if [[ ! -f "$CASES_FILE" ]]; then
        echo "[ERROR] cases.txt not found and CASE_ID not set." >&2
        exit 1
    fi
    CASE_ID=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CASES_FILE")
fi

if [[ -z "$CASE_ID" ]]; then
    echo "[ERROR] Could not determine CASE_ID." >&2
    exit 1
fi

echo "========================================"
echo " Value Extraction Pipeline"
echo " Case ID : $CASE_ID"
echo " Node    : $(hostname)"
echo " Date    : $(date)"
echo "========================================"

# ─── DERIVED PATHS ───────────────────────────────────────────────────────────
SIM_FILE="${CASES_ROOT}/${CASE_ID}/Simulation2D/Simulation.vtkhdf"
SWC_DIR="${CASES_ROOT}/${CASE_ID}/Simulation2D/merged/swc"
RAW_CSV="${OUT_ROOT}/${CASE_ID}_fields.csv"
ENRICHED_CSV="${OUT_ROOT}/${CASE_ID}_fields_enriched.csv"
FIG_DIR="${OUT_ROOT}/figures/${CASE_ID}"
FA_FIG_DIR="${FIG_DIR}/front_aware_analysis"

mkdir -p "${OUT_ROOT}" "${FIG_DIR}" "${FA_FIG_DIR}"

if [[ ! -f "$SIM_FILE" ]]; then
    echo "[ERROR] Simulation file not found: $SIM_FILE" >&2
    exit 1
fi
if [[ ! -d "$SWC_DIR" ]]; then
    echo "[ERROR] SWC directory not found: $SWC_DIR" >&2
    exit 1
fi

# ─── STAGE 1 — Field sampling (pvpython) ─────────────────────────────────────
echo ""
echo "[stage 1] Sampling simulation fields on SWC skeleton ..."
module load paraview/6.1.0-mesa

pvpython "${SCRIPT_DIR}/sample_fields_on_swc.py" \
    --sim-file   "$SIM_FILE" \
    --swc-dir    "$SWC_DIR" \
    --out-csv    "$RAW_CSV" \
    --case-id    "$CASE_ID" \
    --swc-type   "$SWC_TYPE" \
    --smooth-sigma "$SMOOTH_SIGMA"

echo "[stage 1] Done → $RAW_CSV  ($(wc -l < "$RAW_CSV") lines)"

# ─── STAGES 2-4 — Enrichment + plotting (swc_env conda) ──────────────────────
module unload paraview/6.1.0-mesa 2>/dev/null || true

source "$CONDA_SH"
conda activate swc_env

echo ""
echo "[stage 2] Enriching CSV with morpho-dynamic columns ..."
python "${SCRIPT_DIR}/enrich_swc_samples.py" \
    --csv         "$RAW_CSV" \
    --out         "$ENRICHED_CSV" \
    --swc-dir     "$SWC_DIR" \
    --swc-type    "$SWC_TYPE" \
    --r-gc        "$R_GC"

echo "[stage 2] Done → $ENRICHED_CSV  ($(wc -l < "$ENRICHED_CSV") lines)"

echo ""
echo "[stage 3] Generating overview figures ..."
python "${SCRIPT_DIR}/plot_swc_field_samples.py" \
    --csv      "$RAW_CSV" \
    --out-dir  "$FIG_DIR" \
    --case-id  "$CASE_ID"

echo "[stage 3] Done → $FIG_DIR"

echo ""
echo "[stage 4] Generating front-aware analysis figures ..."
python "${SCRIPT_DIR}/plot_front_aware_analysis.py" \
    --csv     "$ENRICHED_CSV" \
    --out-dir "$FA_FIG_DIR" \
    --case-id "$CASE_ID" \
    --r-gc    "$R_GC"

echo "[stage 4] Done → $FA_FIG_DIR"

echo ""
echo "========================================"
echo " Pipeline complete for $CASE_ID"
echo " $(date)"
echo "========================================"
