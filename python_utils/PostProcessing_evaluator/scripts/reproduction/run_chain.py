#!/usr/bin/env python3
"""
run_chain.py
============
Evaluator orchestrator for the 4-step SWC-sampling chain that produced
/scratch/flore0a/AnalysisResults.

Original sources (frozen reference, do NOT edit there):
    copied_reference/old_analysisresults_scripts/sample_fields_on_swc.py
    copied_reference/old_analysisresults_scripts/enrich_swc_samples.py
    copied_reference/old_analysisresults_scripts/plot_swc_field_samples.py
    copied_reference/old_analysisresults_scripts/plot_front_aware_analysis.py

This orchestrator was created on 2026-05-13 for the PostProcessing_evaluator.
Changes vs the original 2026-05-03 ad-hoc invocations:
    - Output goes under a timestamped reproduction_runs/<TS>/ directory
      instead of the flat /scratch/flore0a/AnalysisResults/.
    - An execution_manifest.json captures every command, executable,
      hostname, input path, and output path.
    - Stdout/stderr of each step is logged to <TS>/logs/<step>.{out,err}.
    - The four frozen scripts are invoked directly from copied_reference/
      via -m runpy.run_path semantics with sys.argv overridden — they are
      NOT edited. This keeps copied_reference/ truly frozen.

Why the change is evaluator-only:
    The original AnalysisResults outputs are valuable historical reference;
    we must not overwrite them. Reproductions go into a fresh timestamped
    directory so the historical snapshot remains intact and the evaluator
    can compare new vs old.

The orchestrator does not invent any analysis — it only reschedules the
existing scripts and captures provenance.
"""

import argparse
import json
import os
import shlex
import socket
import subprocess
import sys
import platform
from datetime import datetime
from pathlib import Path


EVAL_ROOT = Path(__file__).resolve().parents[2]
FROZEN_DIR = EVAL_ROOT / 'copied_reference' / 'old_analysisresults_scripts'

SAMPLE       = FROZEN_DIR / 'sample_fields_on_swc.py'
ENRICH       = FROZEN_DIR / 'enrich_swc_samples.py'
PLOT_OVERV   = FROZEN_DIR / 'plot_swc_field_samples.py'
PLOT_FRONTA  = FROZEN_DIR / 'plot_front_aware_analysis.py'


def parse_args():
    ap = argparse.ArgumentParser(
        description="Reproduce the AnalysisResults SWC-sampling chain into a "
                    "timestamped reproduction_runs/<TS>/ directory.")
    ap.add_argument('--case-id', default='D7p5_V3')
    ap.add_argument('--sim-file',
                    default='/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/Simulation.vtkhdf')
    ap.add_argument('--swc-dir',
                    default='/scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc')
    ap.add_argument('--repro-root',
                    default='/scratch/flore0a/PostProcessing_evaluator/reproduction_runs')
    ap.add_argument('--pvpython',
                    default=os.environ.get('PVPYTHON', 'pvpython'),
                    help="Path to pvpython used for the sampling step.")
    ap.add_argument('--swc-python',
                    default=os.environ.get('SWC_PYTHON',
                        '/scratch/flore0a/iops/miniconda3/envs/swc_env/bin/python'),
                    help="Conda Python used for enrichment and plotting.")
    ap.add_argument('--skeleton-variant',
                    choices=['clean', 'pixel', 'both'], default='both',
                    help="Which skeleton variant(s) to reproduce.")
    ap.add_argument('--smooth-sigma', default='1.5',
                    help="Pass-through to sample_fields_on_swc.py for pixel variant.")
    ap.add_argument('--dry-run', action='store_true',
                    help="Print the commands; do not run them.")
    return ap.parse_args()


def run_step(label, cmd, logs_dir, manifest):
    out_path = logs_dir / f'{label}.out'
    err_path = logs_dir / f'{label}.err'
    started = datetime.now().isoformat(timespec='seconds')
    cmd_str = ' '.join(shlex.quote(str(c)) for c in cmd)
    print(f"\n=== {label} ===\n  $ {cmd_str}\n  log → {out_path}", flush=True)
    rc = -1
    with open(out_path, 'w') as oh, open(err_path, 'w') as eh:
        proc = subprocess.run(cmd, stdout=oh, stderr=eh)
        rc = proc.returncode
    finished = datetime.now().isoformat(timespec='seconds')
    manifest['steps'].append({
        'label': label,
        'command': cmd_str,
        'argv': [str(c) for c in cmd],
        'returncode': rc,
        'started_at': started,
        'finished_at': finished,
        'stdout_log': str(out_path),
        'stderr_log': str(err_path),
    })
    if rc != 0:
        print(f"  ERROR: step {label!r} returned {rc}. See {err_path}.", flush=True)
    return rc


def main():
    args = parse_args()
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir = Path(args.repro_root).resolve() / ts
    run_dir.mkdir(parents=True, exist_ok=False)
    logs_dir = run_dir / 'logs'
    logs_dir.mkdir()

    manifest = {
        'evaluator': 'PostProcessing_evaluator',
        'run_dir': str(run_dir),
        'started_at': datetime.now().isoformat(timespec='seconds'),
        'hostname': socket.gethostname(),
        'platform': platform.platform(),
        'python': sys.executable,
        'argv': sys.argv,
        'case_id': args.case_id,
        'inputs': {
            'sim_file': args.sim_file,
            'swc_dir': args.swc_dir,
        },
        'executables': {
            'pvpython': args.pvpython,
            'swc_python': args.swc_python,
        },
        'skeleton_variants_requested': args.skeleton_variant,
        'frozen_scripts': {
            'sample_fields_on_swc.py': str(SAMPLE),
            'enrich_swc_samples.py':   str(ENRICH),
            'plot_swc_field_samples.py':     str(PLOT_OVERV),
            'plot_front_aware_analysis.py':  str(PLOT_FRONTA),
        },
        'steps': [],
        'outputs': {},
    }

    for name in ('sample_fields_on_swc.py', 'enrich_swc_samples.py',
                 'plot_swc_field_samples.py', 'plot_front_aware_analysis.py'):
        p = FROZEN_DIR / name
        if not p.is_file():
            print(f"FATAL: frozen reference missing: {p}", file=sys.stderr)
            sys.exit(2)

    variants = (['clean', 'pixel'] if args.skeleton_variant == 'both'
                else [args.skeleton_variant])

    rc_any_fail = 0
    overview_done = False

    for variant in variants:
        tag = 'clean' if variant == 'clean' else 'pixel'
        raw_csv = run_dir / f'{args.case_id}_{tag}_fields.csv' if variant == 'pixel' \
                  else run_dir / f'{args.case_id}_fields.csv'
        enriched = run_dir / f'{args.case_id}_{tag}_fields_enriched.csv' if variant == 'pixel' \
                   else run_dir / f'{args.case_id}_fields_enriched.csv'
        out_subdir = run_dir / 'figures' / args.case_id / ('front_aware_pixel' if variant == 'pixel' else 'front_aware_analysis')
        out_subdir.mkdir(parents=True, exist_ok=True)
        overview_dir = run_dir / 'figures' / args.case_id
        overview_dir.mkdir(parents=True, exist_ok=True)

        # Step 1 — sample
        sample_cmd = [
            args.pvpython, str(SAMPLE),
            '--sim-file', args.sim_file,
            '--swc-dir',  args.swc_dir,
            '--out-csv',  str(raw_csv),
            '--case-id',  args.case_id,
            '--swc-type', variant,
        ]
        if variant == 'pixel':
            sample_cmd += ['--smooth-sigma', args.smooth_sigma]
        if args.dry_run:
            manifest['steps'].append({'label': f'sample_{variant}', 'command': ' '.join(map(shlex.quote, sample_cmd)), 'dry_run': True})
        else:
            rc = run_step(f'sample_{variant}', sample_cmd, logs_dir, manifest)
            if rc != 0:
                rc_any_fail = rc; continue

        # Step 2 — enrich
        enrich_cmd = [
            args.swc_python, str(ENRICH),
            '--csv', str(raw_csv),
            '--out', str(enriched),
            '--swc-dir', args.swc_dir,
            '--r-gc', '0.1',
        ]
        if args.dry_run:
            manifest['steps'].append({'label': f'enrich_{variant}', 'command': ' '.join(map(shlex.quote, enrich_cmd)), 'dry_run': True})
        else:
            rc = run_step(f'enrich_{variant}', enrich_cmd, logs_dir, manifest)
            if rc != 0:
                rc_any_fail = rc; continue

        # Step 3 — overview figures (clean only; same as the original session)
        if variant == 'clean' and not overview_done:
            overview_cmd = [
                args.swc_python, str(PLOT_OVERV),
                '--csv', str(raw_csv),
                '--out-dir', str(overview_dir),
                '--case-id', args.case_id,
            ]
            if args.dry_run:
                manifest['steps'].append({'label': 'plot_overview', 'command': ' '.join(map(shlex.quote, overview_cmd)), 'dry_run': True})
            else:
                rc = run_step('plot_overview', overview_cmd, logs_dir, manifest)
                if rc != 0:
                    rc_any_fail = rc
            overview_done = True

        # Step 4 — front-aware figures
        fa_cmd = [
            args.swc_python, str(PLOT_FRONTA),
            '--csv', str(enriched),
            '--out-dir', str(out_subdir),
        ]
        if args.dry_run:
            manifest['steps'].append({'label': f'plot_frontaware_{variant}', 'command': ' '.join(map(shlex.quote, fa_cmd)), 'dry_run': True})
        else:
            rc = run_step(f'plot_frontaware_{variant}', fa_cmd, logs_dir, manifest)
            if rc != 0:
                rc_any_fail = rc

        manifest['outputs'][variant] = {
            'raw_csv': str(raw_csv),
            'enriched_csv': str(enriched),
            'figures_dir': str(out_subdir),
        }

    manifest['finished_at'] = datetime.now().isoformat(timespec='seconds')
    manifest['final_returncode'] = rc_any_fail
    manifest['status'] = 'ok' if rc_any_fail == 0 else 'errors'

    with open(run_dir / 'execution_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nreproduction_run: {run_dir}")
    print(f"manifest:         {run_dir/'execution_manifest.json'}")
    if rc_any_fail != 0:
        print("Some step(s) failed; see logs and execution_manifest.json.")
    sys.exit(rc_any_fail)


if __name__ == '__main__':
    main()
