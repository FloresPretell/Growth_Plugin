# copied_reference/ — Provenance map

All files in this folder are **frozen** byte-for-byte copies of the official PostProcessing pipeline. They must not be edited. To adapt a script for evaluator use, copy it again into `../scripts/reproduction/` (for AnalysisResults reproduction) or `../scripts/` (for new evaluator tools), and document the change in the script header (original path, copy date, change list, reason).

Source root: `/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing/`

## Files (this directory)

| Frozen copy here | Original source |
|---|---|
| `main_export.py` | `Postanalysis/main_export.py` |
| `util_postanalisis_export_png.py` | `Postanalysis/util_postanalisis_export_png.py` |
| `Submite_export.slurm` | `Postanalysis/Submite_export.slurm` |
| `main_export_swc.py` | `Postanalysis_swc/main_export_swc.py` |
| `util_postanalisis_export_swc.py` | `Postanalysis_swc/util_postanalisis_export_swc.py` |
| `compute_morphometrics.py` | `Postanalysis_swc/compute_morphometrics.py` |
| `compute_morphometrics_morphoanalysis_plot_variant.py` | `Morphoanalysis_plot/compute_morphometrics.py` (renamed: same name as above in source, different content) |
| `Submit_export_swc.slurm` | `Postanalysis_swc/Submit_export_swc.slurm` |
| `analyze_morphometrics.py` | `Morphoanalysis_plot/analyze_morphometrics.py` |
| `matrix.py` | `Matrix_viosualization_png/matrix.py` |
| `stage2_matrix.slurm` | `slurm/stage2_matrix.slurm` |
| `stage3_morphometrics.slurm` | `slurm/stage3_morphometrics.slurm` |
| `stage4_analysis.slurm` | `slurm/stage4_analysis.slurm` |
| `run_shaheen.sh` | `run_shaheen.sh` |
| `config_shaheen.sh` | `config_shaheen.sh` |

## Files (old_analysisresults_scripts/)

These four scripts produced everything currently sitting in `/scratch/flore0a/AnalysisResults/` (see `/scratch/flore0a/IA_report/IA_session_report_20260503_153249.md`). Treat them as the reproduction reference.

| Frozen copy | Original source |
|---|---|
| `old_analysisresults_scripts/sample_fields_on_swc.py` | `Postanalysis_swc/sample_fields_on_swc.py` |
| `old_analysisresults_scripts/enrich_swc_samples.py` | `Postanalysis_swc/enrich_swc_samples.py` |
| `old_analysisresults_scripts/plot_swc_field_samples.py` | `Postanalysis_swc/plot_swc_field_samples.py` |
| `old_analysisresults_scripts/plot_front_aware_analysis.py` | `Postanalysis_swc/plot_front_aware_analysis.py` |
| `old_analysisresults_scripts/FIELD_SAMPLING_REPORT.md` | `Postanalysis_swc/FIELD_SAMPLING_REPORT.md` |

Note: the byte-identical duplicates under `PostProcessing/Value_Extraction/` were not re-copied (same content as `Postanalysis_swc/`).

## Rule

Files here are `chmod 0444`. Do not chmod or modify. If a change is needed, copy out, do not edit in place.
