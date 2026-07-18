#!/usr/bin/env python3
"""
main_ch5.py — run all CH5 figure modules for one enriched CSV.

Usage:
    python main_ch5.py \
        --csv  /scratch/flore0a/.../D7p5_V3_fields_enriched.csv \
        --variant clean \
        --out-root /scratch/flore0a/PostProcessing_evaluator/ch5_figures

Stops on the first module that crashes only if --strict is passed; by default
it logs the error and continues with the next module.
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
MODULES = [
    ('fig_A', HERE / 'fig_A_elongation.py'),
    ('fig_B', HERE / 'fig_B_avoidance.py'),
    ('fig_C', HERE / 'fig_C_pruning.py'),
    ('fig_D', HERE / 'fig_D_global.py'),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--variant', default='clean')
    ap.add_argument('--out-root',
                    default='/scratch/flore0a/PostProcessing_evaluator/ch5_figures')
    ap.add_argument('--run-ts', default=None)
    ap.add_argument('--strict', action='store_true')
    args = ap.parse_args()

    csv = Path(args.csv).resolve()
    if not csv.is_file():
        print(f'ERROR: CSV not found: {csv}', file=sys.stderr)
        sys.exit(2)

    ts = args.run_ts or datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir = Path(args.out_root).resolve() / ts / args.variant
    out_dir = run_dir / 'figures'
    logs_dir = run_dir / 'logs'
    out_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        'evaluator': 'PostProcessing_evaluator',
        'iteration': 'CH5',
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

    for label, script in MODULES:
        out_log = logs_dir / f'{label}.out'
        err_log = logs_dir / f'{label}.err'
        cmd = [sys.executable, str(script),
               '--csv', str(csv), '--out-dir', str(out_dir)]
        print(f'\n=== {label} ===\n  $ {" ".join(cmd)}', flush=True)
        started = datetime.now().isoformat(timespec='seconds')
        with open(out_log, 'w') as oh, open(err_log, 'w') as eh:
            proc = subprocess.run(cmd, stdout=oh, stderr=eh, env=env)
            rc = proc.returncode
        finished = datetime.now().isoformat(timespec='seconds')
        manifest['steps'].append({
            'label': label, 'command': cmd, 'returncode': rc,
            'started_at': started, 'finished_at': finished,
            'stdout_log': str(out_log), 'stderr_log': str(err_log),
        })
        if rc != 0:
            rc_final = rc
            print(f'  ERROR: rc={rc} — see {err_log}')
            if args.strict:
                break

    manifest['finished_at'] = datetime.now().isoformat(timespec='seconds')
    manifest['status'] = 'ok' if rc_final == 0 else 'errors'
    with open(run_dir / 'ch5_execution_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print(f'\nCH5 run: {run_dir}')
    print(f'manifest: {run_dir / "ch5_execution_manifest.json"}')
    sys.exit(rc_final)


if __name__ == '__main__':
    main()
