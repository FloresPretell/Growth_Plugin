#!/usr/bin/env python3
"""
inspect_simulation_fields.py
============================
Read-only inspection of one simulation directory. Detects the input format
(VTKHDF preferred, VTU/PVTU fallback), enumerates timesteps, lists every
point-data and cell-data array, and reports per-field min/max/mean/std/N/NaN/Inf.

Designed to run under either:
    - pvpython (ParaView's bundled Python; VTK is always present)
    - a plain python3 that has the `vtk` package installed

No field names are hard-coded. Whatever the file actually contains is what
gets reported.

Outputs (all timestamped) under <output_root>/<TS>/:
    field_inventory.csv     long format: one row per (timestep_index, field, kind)
    field_inventory.md      human-readable summary
    run_metadata.json       command, executable, hostname, paths, detected source
    log.txt                 stdout/stderr mirror

Usage:
    pvpython inspect_simulation_fields.py \\
        --input /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D \\
        --output-root /scratch/flore0a/PostProcessing_evaluator/runs

    python3 inspect_simulation_fields.py \\
        --input /scratch/.../Simulation2D \\
        --output-root /scratch/.../runs \\
        --max-timesteps 3
"""

import argparse
import csv
import datetime as dt
import glob
import json
import math
import os
import platform
import re
import socket
import sys
import traceback
from pathlib import Path


# ---------------------------------------------------------------------------
# VTK import (works both in pvpython and plain python with `vtk` installed)
# ---------------------------------------------------------------------------

def import_vtk():
    """
    Import VTK. numpy is intentionally NOT required: the current Shaheen III
    ParaView build (pvpython 6.1.0-mesa, May 2026) ships without numpy in its
    bundled Python, and the swc_env conda has numpy but no vtk. To work in
    both environments without environment changes, this script computes
    statistics in pure Python (Welford one-pass) plus native VTK
    `GetRange(component)` for min/max.
    """
    try:
        import vtk
        return vtk
    except Exception as exc:
        sys.stderr.write(
            "ERROR: Could not import vtk. Run via pvpython, or pip install vtk.\n"
            f"  underlying error: {exc!r}\n"
        )
        sys.exit(2)


# ---------------------------------------------------------------------------
# Source detection
# ---------------------------------------------------------------------------

def detect_source(input_dir: Path):
    """
    Return a dict describing what was found in `input_dir`:
        kind: 'vtkhdf' | 'vtu_series' | 'pvtu_series' | 'unknown'
        primary_path:   path to canonical input
        time_files:     list of (timestep_index_int_or_None, path) ordered by index
        merged_dir:     optional /merged/ subdir
    """
    info = {
        'kind': 'unknown',
        'primary_path': None,
        'time_files': [],
        'merged_dir': None,
        'notes': [],
    }

    # 1) Direct VTKHDF
    vtkhdf = input_dir / 'Simulation.vtkhdf'
    if vtkhdf.is_file() and vtkhdf.stat().st_size > 0:
        info['kind'] = 'vtkhdf'
        info['primary_path'] = str(vtkhdf)
        info['notes'].append('Found Simulation.vtkhdf')
        return info

    # 2) Merged/ subfolder
    merged = input_dir / 'merged'
    if merged.is_dir():
        info['merged_dir'] = str(merged)
        info['notes'].append(f'Found merged/ directory: {merged}')
        # Prefer pvtu series, fallback to vtu
        pvtu = sorted(merged.glob('Output_t*.pvtu'))
        if pvtu:
            info['kind'] = 'pvtu_series'
            info['primary_path'] = str(merged)
            info['time_files'] = [(_extract_step(p), str(p)) for p in pvtu]
            return info
        vtu = sorted(merged.glob('Output_t*.vtu'))
        if vtu:
            info['kind'] = 'vtu_series'
            info['primary_path'] = str(merged)
            info['time_files'] = [(_extract_step(p), str(p)) for p in vtu]
            return info

    # 3) Bare Output_t*.{vtu,pvtu} directly in input_dir
    pvtu = sorted(input_dir.glob('Output_t*.pvtu'))
    if pvtu:
        info['kind'] = 'pvtu_series'
        info['primary_path'] = str(input_dir)
        info['time_files'] = [(_extract_step(p), str(p)) for p in pvtu]
        return info
    vtu = sorted(input_dir.glob('Output_t*.vtu'))
    if vtu:
        info['kind'] = 'vtu_series'
        info['primary_path'] = str(input_dir)
        info['time_files'] = [(_extract_step(p), str(p)) for p in vtu]
        return info

    # 4) Generic .vtu/.pvtu anywhere in input_dir
    pvtu = sorted(input_dir.glob('*.pvtu'))
    if pvtu:
        info['kind'] = 'pvtu_series'
        info['primary_path'] = str(input_dir)
        info['time_files'] = [(_extract_step(p), str(p)) for p in pvtu]
        info['notes'].append('No Output_t*.pvtu found; falling back to *.pvtu glob')
        return info
    vtu = sorted(input_dir.glob('*.vtu'))
    if vtu:
        info['kind'] = 'vtu_series'
        info['primary_path'] = str(input_dir)
        info['time_files'] = [(_extract_step(p), str(p)) for p in vtu]
        info['notes'].append('No Output_t*.vtu found; falling back to *.vtu glob')
        return info

    return info


_step_re = re.compile(r'_t(\d+)')


def _extract_step(path):
    m = _step_re.search(Path(path).name)
    return int(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# VTK loaders
# ---------------------------------------------------------------------------

def open_vtkhdf(vtk_mod, vtkhdf_path):
    """
    Return (reader, n_steps, timestep_values).
    Uses vtkHDFReader (VTK >= 9.2, ParaView >= 5.11).
    """
    rdr = vtk_mod.vtkHDFReader()
    rdr.SetFileName(str(vtkhdf_path))
    rdr.UpdateInformation()
    info = rdr.GetOutputInformation(0)
    n_steps = info.Length(vtk_mod.vtkStreamingDemandDrivenPipeline.TIME_STEPS()) \
        if info.Has(vtk_mod.vtkStreamingDemandDrivenPipeline.TIME_STEPS()) else 0
    if n_steps > 0:
        time_values = [
            info.Get(vtk_mod.vtkStreamingDemandDrivenPipeline.TIME_STEPS(), i)
            for i in range(n_steps)
        ]
    else:
        time_values = [0.0]
        n_steps = 1
    return rdr, n_steps, time_values


def open_one_vtu(vtk_mod, path):
    """Open one VTU or PVTU file and return its dataset (after Update())."""
    p = str(path)
    if p.endswith('.pvtu'):
        rdr = vtk_mod.vtkXMLPUnstructuredGridReader()
    else:
        rdr = vtk_mod.vtkXMLUnstructuredGridReader()
    rdr.SetFileName(p)
    rdr.Update()
    return rdr.GetOutput()


def get_vtkhdf_dataset_at_step(vtk_mod, rdr, time_value):
    rdr.UpdateTimeStep(time_value)
    rdr.Update()
    return rdr.GetOutput()


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def _scan_component(arr, n_tuples, comp_index, stride=1):
    """
    One-pass Welford scan over component `comp_index` of a vtkDataArray.
    Iterates with `stride` (1 = exact; >1 = subsample). Pure Python; numpy-free.
    Returns (mean, var, n_finite, n_nan, n_inf, min_finite, max_finite).
    """
    m = 0.0; s = 0.0
    n_fin = n_nan = n_inf = 0
    mn = float('inf'); mx = float('-inf')
    gc = arr.GetComponent
    for i in range(0, n_tuples, stride):
        v = gc(i, comp_index)
        if v != v: n_nan += 1; continue
        if v == float('inf') or v == float('-inf'): n_inf += 1; continue
        n_fin += 1
        if v < mn: mn = v
        if v > mx: mx = v
        d = v - m; m += d / n_fin; s += d * (v - m)
    var = (s / (n_fin - 1)) if n_fin > 1 else 0.0
    if n_fin == 0: mn = None; mx = None
    return m, var, n_fin, n_nan, n_inf, mn, mx


def _scan_magnitude(arr, n_tuples, n_comp, stride=1):
    """Welford on magnitude (subsampled by `stride`)."""
    m = 0.0; s = 0.0
    n_fin = n_nan = n_inf = 0
    mn = float('inf'); mx = float('-inf')
    gt = arr.GetTuple
    for i in range(0, n_tuples, stride):
        t = gt(i)
        bad_nan = bad_inf = False; sq = 0.0
        for v in t:
            if v != v: bad_nan = True; break
            if v == float('inf') or v == float('-inf'): bad_inf = True; break
            sq += v * v
        if bad_nan: n_nan += 1; continue
        if bad_inf: n_inf += 1; continue
        v = sq ** 0.5
        n_fin += 1
        if v < mn: mn = v
        if v > mx: mx = v
        d = v - m; m += d / n_fin; s += d * (v - m)
    var = (s / (n_fin - 1)) if n_fin > 1 else 0.0
    if n_fin == 0: mn = None; mx = None
    return m, var, n_fin, n_nan, n_inf, mn, mx


def stats_for_array(arr, mode='fast', stride=1):
    """
    Return list of stat-dicts for a single VTK array.

    mode='fast' (default): use native VTK GetRange() for (min, max) per
        component, no mean/std/NaN/Inf scan. Suitable for field-presence checks.
        Time: O(1) per array.

    mode='full': Welford one-pass over every n-th value (stride>=1). Produces
        mean/std/n_nan/n_inf. Default stride=1 is exact but slow in pure Python.

    Vector / tensor arrays are reported as one magnitude row plus one row per
    component.
    """
    name = arr.GetName()
    n_tuples = arr.GetNumberOfTuples()
    n_comp = arr.GetNumberOfComponents()
    rows = []
    if n_comp == 1:
        if mode == 'fast':
            mn, mx = arr.GetRange(0)
            rows.append(_pack_row(name, 'scalar', n_tuples, n_comp,
                                   None, None, None, None, None, mn, mx,
                                   note='fast (GetRange only)'))
        else:
            r = _scan_component(arr, n_tuples, 0, stride)
            rows.append(_pack_row(name, 'scalar', n_tuples, n_comp, *r))
    else:
        if mode == 'fast':
            # magnitude min/max via vtkDataArray::GetRange(-1) (the -1 component
            # means "vector magnitude" in VTK)
            try:
                mn, mx = arr.GetRange(-1)
            except Exception:
                mn = mx = None
            rows.append(_pack_row(name, 'magnitude', n_tuples, n_comp,
                                   None, None, None, None, None, mn, mx,
                                   note='fast (GetRange -1 = magnitude)'))
            for c in range(n_comp):
                mn, mx = arr.GetRange(c)
                rows.append(_pack_row(name, f'c{c}', n_tuples, n_comp,
                                       None, None, None, None, None, mn, mx,
                                       note='fast (GetRange only)'))
        else:
            r = _scan_magnitude(arr, n_tuples, n_comp, stride)
            rows.append(_pack_row(name, 'magnitude', n_tuples, n_comp, *r))
            for c in range(n_comp):
                r = _scan_component(arr, n_tuples, c, stride)
                rows.append(_pack_row(name, f'c{c}', n_tuples, n_comp, *r))
    return rows


def _pack_row(name, component, n_tuples, n_comp,
              mean, var, n_fin, n_nan, n_inf, mn, mx, note=None):
    row = {
        'name': name, 'component': component,
        'n_tuples': int(n_tuples), 'n_components': int(n_comp),
        'min':  (float(mn) if mn is not None else None),
        'max':  (float(mx) if mx is not None else None),
        'mean': (float(mean) if mean is not None else None),
        'std':  (float(var ** 0.5) if (var is not None and var >= 0) else None),
        'n_nan': (int(n_nan) if n_nan is not None else None),
        'n_inf': (int(n_inf) if n_inf is not None else None),
        'note':  note if note is not None else (
            '' if (n_nan == 0 and n_inf == 0)
            else 'contains non-finite'),
    }
    if n_fin is not None and n_fin == 0 and note is None:
        row['note'] = 'all values non-finite'
    return row


def inventory_dataset(dataset, mode='fast', stride=1):
    """Return list of dicts: rows for point-data + cell-data + field-data."""
    rows = []
    pd_ = dataset.GetPointData()
    cd_ = dataset.GetCellData()
    fd_ = dataset.GetFieldData()
    n_points = dataset.GetNumberOfPoints()
    n_cells = dataset.GetNumberOfCells()

    def _emit(kind_label, container):
        for i in range(container.GetNumberOfArrays()):
            arr = container.GetArray(i)
            if arr is None:
                continue
            for r in stats_for_array(arr, mode=mode, stride=stride):
                r['kind'] = kind_label
                r['n_points'] = n_points
                r['n_cells'] = n_cells
                rows.append(r)

    _emit('point_data', pd_)
    _emit('cell_data', cd_)
    _emit('field_data', fd_)
    return rows


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

CSV_COLUMNS = [
    'time_step_index', 'time_value', 'kind', 'name', 'component',
    'n_points', 'n_cells', 'n_tuples', 'n_components',
    'min', 'max', 'mean', 'std', 'n_nan', 'n_inf', 'note',
]


def write_csv(rows, csv_path):
    with open(csv_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_md(rows, md_path, source_info, time_values, started_at, finished_at):
    by_step = {}
    for r in rows:
        by_step.setdefault(r['time_step_index'], []).append(r)
    with open(md_path, 'w') as fh:
        fh.write(f"# Field inventory — generated {finished_at}\n\n")
        fh.write(f"- Input:        `{source_info.get('input_dir')}`\n")
        fh.write(f"- Detected:     `{source_info.get('kind')}`\n")
        fh.write(f"- Primary path: `{source_info.get('primary_path')}`\n")
        fh.write(f"- Timesteps:    {len(time_values)} "
                 f"(values: {time_values[:5]}{' …' if len(time_values) > 5 else ''})\n")
        fh.write(f"- Started:      {started_at}\n")
        fh.write(f"- Finished:     {finished_at}\n\n")

        # Cross-step summary table: fields observed
        seen = {}
        for r in rows:
            k = (r['kind'], r['name'], r.get('n_components'))
            if k not in seen:
                seen[k] = 0
            seen[k] += 1
        fh.write("## Fields observed (across all timesteps)\n\n")
        fh.write("| kind | field | n_components | timesteps_seen |\n")
        fh.write("|------|-------|-------------:|---------------:|\n")
        for (kind, name, ncomp), count in sorted(seen.items()):
            fh.write(f"| {kind} | `{name}` | {ncomp} | {count} |\n")
        fh.write("\n")

        # Per-timestep summary (compact)
        fh.write("## Per-timestep summary (scalar / magnitude rows only)\n\n")
        for step in sorted(by_step.keys(),
                           key=lambda x: (x is None, x if x is not None else -1)):
            entries = by_step[step]
            t_vals = sorted({e.get('time_value') for e in entries
                             if e.get('time_value') is not None})
            t_str = f", t={t_vals[0]}" if t_vals else ''
            n_points = entries[0].get('n_points')
            n_cells = entries[0].get('n_cells')
            fh.write(f"### timestep_index={step}{t_str}  "
                     f"(n_points={n_points}, n_cells={n_cells})\n\n")
            fh.write("| kind | field | comp | min | max | mean | std | nan | inf |\n")
            fh.write("|------|-------|------|---:|---:|---:|---:|---:|---:|\n")
            for e in entries:
                if e['component'] not in ('scalar', 'magnitude'):
                    continue
                fh.write(f"| {e['kind']} | `{e['name']}` | {e['component']} "
                         f"| {_fmt(e['min'])} | {_fmt(e['max'])} "
                         f"| {_fmt(e['mean'])} | {_fmt(e['std'])} "
                         f"| {e['n_nan']} | {e['n_inf']} |\n")
            fh.write("\n")


def _fmt(x):
    if x is None:
        return '—'
    if isinstance(x, float):
        if not math.isfinite(x):
            return str(x)
        return f"{x:.6g}"
    return str(x)


def write_metadata(meta, json_path):
    with open(json_path, 'w') as fh:
        json.dump(meta, fh, indent=2, default=str)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Field inventory of one simulation directory.")
    ap.add_argument('--input', required=True,
                    help="Simulation directory (containing Simulation.vtkhdf "
                         "or merged/Output_t*.{vtu,pvtu}).")
    ap.add_argument('--output-root', required=True,
                    help="Root under which a timestamped run directory is created.")
    ap.add_argument('--case-id', default='',
                    help="Optional case identifier embedded in run metadata.")
    ap.add_argument('--max-timesteps', type=int, default=None,
                    help="Limit the number of timesteps inspected.")
    ap.add_argument('--timestep-stride', type=int, default=1,
                    help="Inspect every Nth timestep (default 1).")
    ap.add_argument('--mode', choices=['fast', 'full'], default='fast',
                    help="'fast' (default): native VTK GetRange() only, "
                         "no mean/std/NaN/Inf scan. 'full': pure-Python "
                         "Welford scan (slow without numpy).")
    ap.add_argument('--stats-stride', type=int, default=1,
                    help="Welford subsample stride when --mode=full "
                         "(default 1 = exact).")
    args = ap.parse_args()

    vtk = import_vtk()

    started_at = dt.datetime.now().isoformat(timespec='seconds')
    ts_dir_name = dt.datetime.now().strftime('%Y%m%d_%H%M%S')
    output_root = Path(args.output_root).resolve()
    run_dir = output_root / ts_dir_name
    run_dir.mkdir(parents=True, exist_ok=False)

    log_path = run_dir / 'log.txt'
    log_fh = open(log_path, 'w')
    class Tee:
        def __init__(self, *streams):
            self.streams = streams
        def write(self, s):
            for st in self.streams:
                st.write(s); st.flush()
        def flush(self):
            for st in self.streams:
                st.flush()
    sys.stdout = Tee(sys.__stdout__, log_fh)
    sys.stderr = Tee(sys.__stderr__, log_fh)

    print(f"[{started_at}] inspect_simulation_fields starting")
    print(f"  python      : {sys.executable}")
    print(f"  vtk         : {getattr(vtk, 'VTK_VERSION', 'unknown')}")
    print(f"  input       : {args.input}")
    print(f"  output_root : {output_root}")
    print(f"  run_dir     : {run_dir}")

    input_dir = Path(args.input).resolve()
    if not input_dir.is_dir():
        print(f"ERROR: --input must be a directory: {input_dir}")
        sys.exit(2)
    source = detect_source(input_dir)
    source['input_dir'] = str(input_dir)
    print(f"  detected    : kind={source['kind']}, primary={source['primary_path']}")
    for note in source['notes']:
        print(f"               note: {note}")

    if source['kind'] == 'unknown':
        print("ERROR: No recognized simulation output found in input directory.")
        print("       Looked for Simulation.vtkhdf, merged/Output_t*.{vtu,pvtu}, and *.{vtu,pvtu}.")
        meta = {
            'started_at': started_at,
            'finished_at': dt.datetime.now().isoformat(timespec='seconds'),
            'command': sys.argv,
            'python': sys.executable,
            'vtk_version': getattr(vtk, 'VTK_VERSION', 'unknown'),
            'hostname': socket.gethostname(),
            'platform': platform.platform(),
            'case_id': args.case_id,
            'input': str(input_dir),
            'detected': source,
            'status': 'no_input_found',
        }
        write_metadata(meta, run_dir / 'run_metadata.json')
        sys.exit(3)

    # Inspect
    all_rows = []
    time_values_used = []
    if source['kind'] == 'vtkhdf':
        rdr, n_steps, time_values = open_vtkhdf(vtk, source['primary_path'])
        indices = list(range(n_steps))
    else:
        # VTU/PVTU series: one file per step
        files = source['time_files']
        indices = list(range(len(files)))
        time_values = [step for (step, _) in files]
        n_steps = len(files)

    if args.timestep_stride > 1:
        indices = indices[::args.timestep_stride]
    if args.max_timesteps is not None:
        indices = indices[:args.max_timesteps]

    print(f"  timesteps   : total={n_steps}, will inspect={len(indices)}")

    for i_step in indices:
        if source['kind'] == 'vtkhdf':
            t_val = time_values[i_step]
            ds = get_vtkhdf_dataset_at_step(vtk, rdr, t_val)
            step_idx = i_step
        else:
            step_idx, fpath = source['time_files'][i_step]
            t_val = step_idx if step_idx is not None else float(i_step)
            ds = open_one_vtu(vtk, fpath)
        rows = inventory_dataset(ds, mode=args.mode, stride=args.stats_stride)
        for r in rows:
            r['time_step_index'] = step_idx
            r['time_value'] = t_val
        all_rows.extend(rows)
        time_values_used.append(t_val)
        scalars = [r for r in rows if r['component'] in ('scalar', 'magnitude')]
        print(f"   step {step_idx} (t={t_val}): "
              f"n_points={rows[0]['n_points'] if rows else '?'}, "
              f"n_cells={rows[0]['n_cells'] if rows else '?'}, "
              f"unique_fields={len(set((r['kind'], r['name']) for r in scalars))}")

    finished_at = dt.datetime.now().isoformat(timespec='seconds')

    csv_path = run_dir / 'field_inventory.csv'
    md_path = run_dir / 'field_inventory.md'
    json_path = run_dir / 'run_metadata.json'
    write_csv(all_rows, csv_path)
    write_md(all_rows, md_path, source, time_values_used, started_at, finished_at)

    meta = {
        'started_at': started_at,
        'finished_at': finished_at,
        'command': sys.argv,
        'python': sys.executable,
        'vtk_version': getattr(vtk, 'VTK_VERSION', 'unknown'),
        'hostname': socket.gethostname(),
        'platform': platform.platform(),
        'case_id': args.case_id,
        'input': str(input_dir),
        'detected': source,
        'timesteps_total': n_steps,
        'timesteps_inspected': len(indices),
        'time_values_inspected': time_values_used,
        'rows_written': len(all_rows),
        'csv': str(csv_path),
        'md': str(md_path),
        'log': str(log_path),
        'status': 'ok',
    }
    write_metadata(meta, json_path)

    print(f"[{finished_at}] DONE.")
    print(f"  csv  : {csv_path}")
    print(f"  md   : {md_path}")
    print(f"  json : {json_path}")
    print(f"  log  : {log_path}")


if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception:
        sys.stderr.write("FATAL:\n")
        traceback.print_exc()
        sys.exit(1)
