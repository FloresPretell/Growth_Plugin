cd /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing

RESULT_ROOT=/scratch/flore0a/Dataset_test3_Results/PostProcessing_20260427_131826
PIPELINE_ROOT=/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing

mkdir -p "${RESULT_ROOT}/logs"
mkdir -p "${RESULT_ROOT}/MatrixVisualization"

echo "Check PNGs:"
find "${RESULT_ROOT}/PNG" -type f -path '*/u_t/*.png' | head -20
find "${RESULT_ROOT}/PNG" -type f -path '*/u_t/*.png' | wc -l

mv "${RESULT_ROOT}/MatrixVisualization/mosaic_matrix.gif" \
   "${RESULT_ROOT}/MatrixVisualization/mosaic_matrix_BAD_previous.gif" 2>/dev/null || true

sbatch \
  --chdir="${RESULT_ROOT}" \
  --output="${RESULT_ROOT}/logs/stage2_rerun_%j.out" \
  --error="${RESULT_ROOT}/logs/stage2_rerun_%j.err" \
  --export=ALL,PIPELINE_ROOT="${PIPELINE_ROOT}",PNG_ROOT="${RESULT_ROOT}/PNG",MATRIX_ROOT="${RESULT_ROOT}/MatrixVisualization",MATRIX_PYTHON=/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python \
  "${PIPELINE_ROOT}/slurm/stage2_matrix.slurm"

## resubmit the 3part of case 3


cd /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing

RESULT_ROOT=/scratch/flore0a/Dataset_test4_Results/PostProcessing_20260427_133226
PIPELINE_ROOT=/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing

mkdir -p "${RESULT_ROOT}/logs"
mkdir -p "${RESULT_ROOT}/Morphometrics"

sbatch \
  --chdir="${RESULT_ROOT}" \
  --output="${RESULT_ROOT}/logs/stage3b_rerun_%j.out" \
  --error="${RESULT_ROOT}/logs/stage3b_rerun_%j.err" \
  --export=ALL,PIPELINE_ROOT="${PIPELINE_ROOT}",SWC_ROOT="${RESULT_ROOT}/SWC",MORPHO_ROOT="${RESULT_ROOT}/Morphometrics",MORPHO_PYTHON=/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python \
  "${PIPELINE_ROOT}/slurm/stage3_morphometrics.slurm"

