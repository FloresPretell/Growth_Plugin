# PostProcessing_evaluator

Evaluator-only workspace for developing and validating postprocessing analysis tools that extract time-dependent values along simulated neuronal geometries. The long-term goal is to identify robust quantities and figures useful for publication.

The **official postprocessing pipeline must not be modified** until the evaluator version is validated.

## What is protected (must not be modified)

- `/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing/` — the official Python postprocessing pipeline (this evaluator's reference).
- `/project/k10070/Nicole/UG4/ug4/apps/Neuronal_Process/PostProcessing/` — the older app-side mirror (treated as read-only too).
- `/project/k10070/Nicole/UG4/ug4/apps/Neuronal_Process/Scripts_lua/Model_only_final_steps.lua` and other production Lua.
- `/scratch/flore0a/AnalysisResults/` — historical postprocessing outputs (read-only reference for reproduction validation).

## Directory layout

```
PostProcessing_evaluator/
├── README.md
├── copied_reference/                    # frozen byte-for-byte copies of official scripts
│   ├── README_PROVENANCE.md             #   provenance map (original path of each file)
│   └── old_analysisresults_scripts/     #   scripts that produced /scratch/flore0a/AnalysisResults
├── scripts/
│   ├── inspect_simulation_fields.py     # Stage 1 — autodetect & inventory all fields in a sim dir
│   └── reproduction/                    # evaluator-adapted variants (see header for changes)
├── configs/
├── sh_files/
│   ├── run_inspect_simulation_fields.sh
│   └── run_reproduce_analysisresults.sh
├── slurm/
│   ├── run_inspect_simulation_fields.slurm
│   └── run_reproduce_analysisresults.slurm
├── reports/
├── presentation/
│   ├── presentation.tex
│   └── compile_presentation.sh
├── tests/
└── logs/
```

Heavy outputs (CSV, figures, intermediate data, run logs) are written to scratch:

```
/scratch/flore0a/PostProcessing_evaluator/
├── inputs/
├── outputs/
├── reports/         # timestamped inspection / inventory / provenance / comparison reports
├── logs/            # SLURM and login-node logs
├── figures/
├── tables/
├── runs/            # /<TS>/  one inspection run per timestamp
├── reproduction_runs/  # /<TS>/  one full reproduction per timestamp
└── temporary/
```

## How to run

### Stage 1 — field inventory of one simulation directory

Login node (uses plain python3; vtk must be installed):
```bash
bash /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator/sh_files/run_inspect_simulation_fields.sh \
     /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D D7p5_V3
```

SLURM (uses pvpython under `srun --ntasks=1`, which is required on Shaheen III — see `/scratch/flore0a/Instruction_every_sessionread.md` §17.5):
```bash
INPUT_DIR=/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D \
CASE_ID=D7p5_V3 \
sbatch /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator/slurm/run_inspect_simulation_fields.slurm
```

Output goes to `/scratch/flore0a/PostProcessing_evaluator/runs/<TS>/`:
- `field_inventory.csv` (long format, one row per timestep × field × component)
- `field_inventory.md` (per-step summary table)
- `run_metadata.json` (command, executable, hostname, vtk version, detected source)
- `log.txt`

### AnalysisResults reproduction (4-stage chain)

Login node — the script orchestrates the 4 SWC-sampling steps end-to-end into a fresh timestamped directory:
```bash
bash /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator/sh_files/run_reproduce_analysisresults.sh
```

SLURM (for the heavy pixel-skeleton sampling):
```bash
sbatch /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/PostProcessing_evaluator/slurm/run_reproduce_analysisresults.slurm
```

Output goes to `/scratch/flore0a/PostProcessing_evaluator/reproduction_runs/<TS>/` with an `execution_manifest.json` capturing the command, inputs, outputs, executables, hostname, and date. The historical `/scratch/flore0a/AnalysisResults/` is **never** touched.

## Recovered workflow — pipeline stages

See [`/scratch/flore0a/PostProcessing_evaluator/reports/recovered_workflow_map_20260513_003846.md`](file:///scratch/flore0a/PostProcessing_evaluator/reports/recovered_workflow_map_20260513_003846.md) for the full classification. Summary:

1. **Field inventory** — `scripts/inspect_simulation_fields.py` (NEW, evaluator-only). Autodetects fields; no hard-coded names.
2. **Field sampling on geometry** — `copied_reference/old_analysisresults_scripts/sample_fields_on_swc.py` (frozen). pvpython + vtkProbeFilter at SWC node coordinates.
3. **Front-aware morphological enrichment** — `copied_reference/old_analysisresults_scripts/enrich_swc_samples.py` (frozen).
4. **Overview visualization** — `copied_reference/old_analysisresults_scripts/plot_swc_field_samples.py` (frozen). Produces fig1–fig4.
5. **Front-aware publication visualization** — `copied_reference/old_analysisresults_scripts/plot_front_aware_analysis.py` (frozen). Produces figA–figF + summary_report.txt.
6. **Morphometric extraction** — `copied_reference/compute_morphometrics.py` (frozen, not exercised yet).
7. **Multivariate / statistical analysis** — `copied_reference/analyze_morphometrics.py` (frozen, not exercised yet; needs multi-case dataset).
8. **Mosaic / matrix visualization** — `copied_reference/matrix.py` (frozen, not exercised yet).
9. **Publication figures** — TBD, will derive from Stage 5 outputs as evaluator iterates.

## How this evaluator differs from the official pipeline

- **No PNG-stage dependency for the first pass.** The official pipeline begins with a PNG export (Stage 1 of `run_shaheen.sh`). For the reproduction path we use the already-existing `Simulation.vtkhdf` directly, skipping PNG export entirely.
- **All outputs are timestamped.** The official pipeline does timestamp its run root (`PostProcessing_${RUN_ID}/`), but the AnalysisResults artifacts under `/scratch/flore0a/AnalysisResults/` are flat (no timestamp). The evaluator's reproduction runs are always under `/scratch/flore0a/PostProcessing_evaluator/reproduction_runs/<TS>/`.
- **`copied_reference/` files are `chmod 0444`.** Source scripts are not in PATH; adapted variants live in `scripts/reproduction/` with explicit header documentation.
- **Field discovery first.** A first-pass field inventory is the prerequisite for any field-specific analysis; the official pipeline assumes a fixed list (`u_u, u_p, u_t, u_b, u_ca_cyt`).

## Current limitations

- The reproduction script targets only `D7p5_V3` because that is the only case present in `AnalysisResults`. Multi-case reproduction is deferred to the next evaluator iteration.
- The pixel-skeleton sampling step is the heavy I/O step (~92,125 vtkProbeFilter calls × 11 fields). On Shaheen this must run under SLURM/pvpython; on the login node it will be I/O-bound and slow.
- No automated pixel-level diff against historical figures yet. The reproduction comparison report compares filenames, sizes, and shapes — visual equivalence is not asserted.
- Plain login-node `python3` does not have `vtk`; the inspect script needs `pvpython` for ParaView's VTK or a conda env with `vtk`. The SLURM script handles this.

## Next planned analyses

- Run inspection on every case in `/scratch/flore0a/Dataset_test{2,3,4}/Cases/*/Simulation2D` to build a cross-case field consistency report.
- Extend `plot_front_aware_analysis.py` with a Lmax(t) plot showing the elongation-arrest plateau visible in `summary_report.txt` (D7.5_V3 plateaus at L_max ≈ 3.05 from frame ~5 onward).
- Cross-validate the SWC-pixel and clean-skeleton variants on the same enriched fields (already done visually in the two front_aware_* subfolders; need a quantitative comparison).
- Add Stage 6 + 7 (morphometrics + multivariate) when more than one case is reproduced.

## Safety boundary

- Official PostProcessing (`/project/.../python_utils/PostProcessing/`): never modified.
- Production Lua: never modified.
- `/scratch/flore0a/AnalysisResults/`: never modified or moved.
- All heavy outputs go to `/scratch/flore0a/PostProcessing_evaluator/`.
- All SLURM scripts include `--mail-user=elsanicole.florespretell@kaust.edu.sa` and `--mail-type=ALL`.

See `/scratch/flore0a/Instruction_every_sessionread.md` for full operational rules on Shaheen III.
