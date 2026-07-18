#!/usr/bin/env python3
"""
run_figure_design.py — driver for the second-generation figure-design iteration.

Reads one enriched CSV and produces a fresh timestamped directory under
/scratch/flore0a/PostProcessing_evaluator/figure_design_runs/<TS>/<variant>/
containing the A–E candidate figures and an execution_manifest.json.

This is exploratory; outputs are clearly marked. No publication claims.

Usage:
    python3 run_figure_design.py \
        --csv /scratch/flore0a/PostProcessing_evaluator/reproduction_runs/<TS>/D7p5_V3_fields_enriched.csv \
        --variant clean \
        --root /scratch/flore0a/PostProcessing_evaluator/figure_design_runs
"""

import argparse
import json
import os
import socket
import subprocess
import sys
from datetime import datetime
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCRIPTS = {
    'A_kymograph':       HERE / 'figdesign_A_kymograph.py',
    'B_gc_timecourse':   HERE / 'figdesign_B_gc_timecourse.py',
    'C_lag_corr':        HERE / 'figdesign_C_lag_corr.py',
    'D_new_region':      HERE / 'figdesign_D_new_region.py',
    'E_branch_summary':  HERE / 'figdesign_E_branch_summary.py',
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True,
                    help="Enriched CSV (produced by enrich_swc_samples.py).")
    ap.add_argument('--variant', default='clean',
                    help="Label for the variant (clean, pixel, ...).")
    ap.add_argument('--root',
                    default='/scratch/flore0a/PostProcessing_evaluator/figure_design_runs')
    ap.add_argument('--run-ts', default=None,
                    help="Optional fixed timestamp (so multiple variants land in one run).")
    args = ap.parse_args()

    csv = Path(args.csv).resolve()
    if not csv.is_file():
        print(f'ERROR: CSV not found: {csv}', file=sys.stderr)
        sys.exit(2)

    ts = args.run_ts or datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir = Path(args.root).resolve() / ts
    out_dir = run_dir / args.variant / 'figures'
    out_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = run_dir / args.variant / 'logs'
    logs_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        'evaluator': 'PostProcessing_evaluator',
        'iteration': 'figure_design',
        'run_dir': str(run_dir),
        'variant': args.variant,
        'csv': str(csv),
        'started_at': datetime.now().isoformat(timespec='seconds'),
        'hostname': socket.gethostname(),
        'python': sys.executable,
        'argv': sys.argv,
        'steps': [],
    }

    rc_final = 0
    env = dict(os.environ)
    env['PYTHONPATH'] = str(HERE) + os.pathsep + env.get('PYTHONPATH', '')

    for label, script in SCRIPTS.items():
        out_log = logs_dir / f'{label}.out'
        err_log = logs_dir / f'{label}.err'
        cmd = [sys.executable, str(script),
               '--csv', str(csv),
               '--out-dir', str(out_dir)]
        started = datetime.now().isoformat(timespec='seconds')
        print(f'\n=== {label} ===', flush=True)
        print(f'  $ {" ".join(cmd)}', flush=True)
        with open(out_log, 'w') as oh, open(err_log, 'w') as eh:
            proc = subprocess.run(cmd, stdout=oh, stderr=eh, env=env)
            rc = proc.returncode
        finished = datetime.now().isoformat(timespec='seconds')
        if rc != 0:
            rc_final = rc
            print(f'  ERROR: rc={rc}; see {err_log}')
        manifest['steps'].append({
            'label': label, 'command': cmd,
            'returncode': rc,
            'started_at': started, 'finished_at': finished,
            'stdout_log': str(out_log), 'stderr_log': str(err_log),
        })

    manifest['finished_at'] = datetime.now().isoformat(timespec='seconds')
    manifest['status'] = 'ok' if rc_final == 0 else 'errors'
    with open(run_dir / args.variant / 'execution_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)

    print(f'\nfigure_design run: {run_dir / args.variant}')
    print(f'manifest:          {run_dir / args.variant / "execution_manifest.json"}')
    sys.exit(rc_final)


if __name__ == '__main__':
    main()
