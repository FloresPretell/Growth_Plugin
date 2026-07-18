#!/usr/bin/env pvpython
"""
sample_fields_on_swc_extended.py
=================================
Adapted from:
  /project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/
  PostProcessing_evaluator/copied_reference/old_analysisresults_scripts/
  sample_fields_on_swc.py    (frozen, chmod 0444 — never edited)

Date copied: 2026-05-13
Author of adaptation: PostProcessing_evaluator session

Changes vs frozen reference (exhaustive list):

  1. DEFAULT_FIELDS extended with three new fields:
        'norm_vel_int'   (point_data scalar) — interface normal velocity, v_n
        'curvature'      (point_data scalar) — level-set curvature, κ
        'GrowthVel'      (cell_data  vec3 )  — cell-data growth velocity vector

  2. vtkProbeFilter call: PassCellArraysOn() so cell-data is interpolated to
     probe points alongside point-data. This makes 'GrowthVel' accessible
     by name on the probed output the same way as point arrays.

  3. Per-node scalar derived column added:
        growth_vel_n = GrowthVel  ·  t̂_node
     where t̂_node is the unit SWC tangent — direction from the parent node
     to the current node (in physical µm), normalised. For the root, the
     direction toward the first child is used. For the rare case of a node
     with no parent and no children, growth_vel_n is set to NaN.

  4. BASE_COLUMNS extended with 'tangent_x', 'tangent_y' (the per-node unit
     tangent components, in physical µm), purely for downstream sanity check.

Nothing else is changed. The pixel→physical transform, probe-filter logic,
SWC parser, smoothing, frame discovery, output schema (CSV writing), and CLI
flags are reused verbatim by importing the helpers from the frozen reference.

Why evaluator-only:
  The official pipeline is validated and must not be changed. This variant
  is required by CH5 §5.3 (velocity-law test with the real v_n instead of
  the dLmax/dt proxy) and §5.6 (curvature-vs-inhibitor test with the real κ
  from the level-set instead of the SWC-triplet Menger curvature).

Outputs:
  /scratch/flore0a/PostProcessing_evaluator/reproduction_runs/<TS>/
  (never to AnalysisResults or official PostProcessing)
"""

# pvpython runs Python 3.12. importlib.util is in the stdlib.
import argparse
import csv
import importlib.util
import math
import sys
from pathlib import Path

import vtk

# ---- Load frozen helpers as a module -----------------------------------------
FROZEN = (
    '/project/k10070/Nicole/UG4/ug4/plugins/Growth_Plugin/python_utils/'
    'PostProcessing_evaluator/copied_reference/old_analysisresults_scripts/'
    'sample_fields_on_swc.py'
)
_spec = importlib.util.spec_from_file_location('frozen_sampler', FROZEN)
frozen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(frozen)

parse_swc            = frozen.parse_swc
compute_topology     = frozen.compute_topology
smooth_pixel_swc     = frozen.smooth_pixel_swc
build_pixel_to_physical = frozen.build_pixel_to_physical
open_vtkhdf          = frozen.open_vtkhdf
open_vtu_series      = frozen.open_vtu_series
get_mesh_vtkhdf      = frozen.get_mesh_vtkhdf
get_mesh_vtu         = frozen.get_mesh_vtu
discover_swc_frames  = frozen.discover_swc_frames


# ---- Extended field list -----------------------------------------------------
DEFAULT_FIELDS_EXT = [
    'u_t', 'u_u', 'u_b', 'u_p', 'u_ca_cyt', 'inhibitor', 'lsf',
    'norm_vel_int', 'curvature',   # newly propagated point-data scalars
    'GrowthVel',                   # cell-data vec3 → scalar via tangent
]

BASE_COLUMNS_EXT = [
    'case_id', 'time_step', 'time', 'swc_frame',
    'swc_node_id', 'parent_id', 'node_type',
    'x_px', 'y_px',
    'x_phys', 'y_phys', 'z_phys',
    'distance_from_soma', 'path_length_from_soma', 'branch_order', 'is_tip',
    'tangent_x', 'tangent_y',
]


# ---- Helper: per-node unit tangent from SWC topology ------------------------
def compute_tangents(nodes):
    """
    For each node in `nodes`, add 'tangent_x', 'tangent_y' (unit vector in
    physical µm). If the node has a parent, the tangent is from parent →
    current; for the root, from current → first child; for isolated nodes,
    NaN.

    Requires nodes already have 'x_phys', 'y_phys' set (after the transform)
    and 'parent', 'id' from the SWC. Adds entries in-place.
    """
    by_id = {n['id']: n for n in nodes}
    children = {n['id']: [] for n in nodes}
    for n in nodes:
        if n['parent'] in children:
            children[n['parent']].append(n['id'])
    NAN = float('nan')
    for n in nodes:
        pid = n['parent']
        if pid in by_id and pid != n['id']:
            par = by_id[pid]
            dx = n['x_phys'] - par['x_phys']
            dy = n['y_phys'] - par['y_phys']
        else:
            childs = children.get(n['id'], [])
            if not childs:
                n['tangent_x'] = NAN
                n['tangent_y'] = NAN
                continue
            ch = by_id[childs[0]]
            dx = ch['x_phys'] - n['x_phys']
            dy = ch['y_phys'] - n['y_phys']
        mag = math.hypot(dx, dy)
        if mag < 1e-12:
            n['tangent_x'] = NAN
            n['tangent_y'] = NAN
        else:
            n['tangent_x'] = dx / mag
            n['tangent_y'] = dy / mag


# ---- Extended probe: point-data + GrowthVel cell-data + dot product ---------
def probe_fields_extended(mesh, phys_coords, scalar_fields, vector_field,
                            tangents):
    """
    Like the frozen probe_fields, but also samples a cell-data 3-vector field
    `vector_field` (e.g. 'GrowthVel') alongside scalar point-data, then
    returns a `growth_vel_n` scalar per node via:
        growth_vel_n = vector · t̂_node      (in 2D: vx·tx + vy·ty)

    scalar_fields:  list of strings to probe as scalar (point or cell data).
    vector_field:   name of the 3-vector cell_data array, or None to skip.
    tangents:       list of (tx, ty, tz=0) per node, parallel to phys_coords.

    Returns list of dicts (one per probe point):
        {<field>: '...', 'valid_point': 1|0, 'growth_vel_n': '...'}
    """
    pts = vtk.vtkPoints()
    for x, y, z in phys_coords:
        pts.InsertNextPoint(x, y, z)
    probe_poly = vtk.vtkPolyData()
    probe_poly.SetPoints(pts)
    probe = vtk.vtkProbeFilter()
    probe.SetSourceData(mesh)
    probe.SetInputData(probe_poly)
    probe.PassPointArraysOn()
    probe.PassCellArraysOn()
    probe.Update()

    out = probe.GetOutput()
    out_pd = out.GetPointData()
    valid_arr = out_pd.GetArray('vtkValidPointMask')

    # Validate that the requested fields exist (warn once per missing field)
    missing_scalar = []
    for fname in scalar_fields:
        if out_pd.GetArray(fname) is None:
            missing_scalar.append(fname)
    if missing_scalar:
        print(f'[warn] probe_extended: missing scalar fields '
              f'(probe returned NaN): {missing_scalar}')
    vec_arr = None
    if vector_field:
        vec_arr = out_pd.GetArray(vector_field)
        if vec_arr is None:
            print(f'[warn] probe_extended: vector field {vector_field!r} '
                  f'not found; growth_vel_n will be empty.')

    n = len(phys_coords)
    rows = []
    for i in range(n):
        valid = bool(int(valid_arr.GetTuple1(i))) if valid_arr is not None else True
        row = {'valid_point': int(valid)}
        # scalar fields
        for fname in scalar_fields:
            arr = out_pd.GetArray(fname)
            if arr is None or not valid:
                row[fname] = ''
            else:
                try:
                    row[fname] = f"{arr.GetValue(i):.8g}"
                except Exception:
                    row[fname] = ''
        # vector → scalar dot with tangent
        if vec_arr is not None and valid:
            try:
                vx, vy, vz = vec_arr.GetTuple(i)
            except Exception:
                vx = vy = vz = float('nan')
            tx, ty, tz = tangents[i]
            if (tx != tx) or (ty != ty):  # NaN check
                row['growth_vel_n'] = ''
            else:
                vn = vx * tx + vy * ty
                row['growth_vel_n'] = f"{vn:.8g}"
        else:
            row['growth_vel_n'] = ''
        # Pass-through the full vector components (vector_field) — useful for QA
        if vec_arr is not None and valid:
            try:
                vx, vy, vz = vec_arr.GetTuple(i)
                row[vector_field] = f"{vx:.6g};{vy:.6g};{vz:.6g}"
            except Exception:
                row[vector_field] = ''
        else:
            row[vector_field] = ''
        rows.append(row)
    return rows


# ---- Main processing (replaces frozen.process_case for the extended path) ---
def process_case(args):
    case_id = args.case_id
    dry_run = args.dry_run
    swc_type = args.swc_type
    image_size = args.image_size
    requested_fields = [f.strip() for f in args.fields.split(',') if f.strip()]

    # Separate scalar vs vector (GrowthVel is a vector; everything else scalar)
    vector_name = 'GrowthVel' if 'GrowthVel' in requested_fields else None
    scalar_fields = [f for f in requested_fields if f != 'GrowthVel']

    # Open simulation reader (delegate to frozen helpers)
    if args.sim_file:
        reader, times, t0_bounds, available_fields = open_vtkhdf(args.sim_file)
        reader_mode = 'vtkhdf'
        vtu_files = None
        vtu_reader = None
    else:
        vtu_files, vtu_reader, t0_bounds, available_fields = open_vtu_series(
            args.vtu_dir, args.vtu_pattern)
        times = list(range(len(vtu_files)))
        reader_mode = 'vtu'
        reader = None

    # Inform about missing requested scalars
    missing = [f for f in scalar_fields if f not in available_fields]
    if missing:
        print(f'[warn] Requested point-data fields not found in simulation: '
              f'{missing}')
        print(f'[warn] Available point-data: {available_fields}')

    transform, S, px_scale = build_pixel_to_physical(t0_bounds, image_size)
    print(f'[info] Pixel scale: {px_scale:.6f} phys_units/pixel  (S = {S:.4f})')

    all_frames = discover_swc_frames(args.swc_dir, swc_type)
    print(f'[info] Found {len(all_frames)} SWC frames (type: {swc_type})')
    if not all_frames:
        sys.exit('[ERROR] No SWC frames found.')

    frames_to_process = all_frames[:1] if dry_run else all_frames

    out_rows = []
    warn_count = 0

    for frame_str, swc_path in frames_to_process:
        time_step_idx = int(frame_str)
        if time_step_idx >= len(times):
            print(f'[warn] Frame {frame_str}: index {time_step_idx} out of '
                  f'range ({len(times)} time steps), skipping.')
            warn_count += 1
            continue
        time_val = times[time_step_idx]

        nodes = parse_swc(swc_path)
        if not nodes:
            print(f'[warn] Frame {frame_str}: empty SWC, skipping.')
            warn_count += 1
            continue
        if swc_type == 'pixel' and args.smooth_sigma > 0:
            nodes = smooth_pixel_swc(nodes, sigma=args.smooth_sigma)
        nodes = compute_topology(nodes)

        if dry_run:
            nodes = nodes[:10]
            print(f'\n[DRY RUN] Frame {frame_str}, time={time_val:.4g}, '
                  f'sampling {len(nodes)} nodes:')

        phys_coords = []
        for n in nodes:
            xp, yp = transform(n['x'], n['y'])
            n['x_phys'], n['y_phys'], n['z_phys'] = xp, yp, 0.0
            n['distance_from_soma'] = n['distance_px'] * px_scale
            n['path_length_from_soma'] = n['path_length_px'] * px_scale
            phys_coords.append((xp, yp, 0.0))

        # Per-node unit tangent (after the transform so we work in µm)
        compute_tangents(nodes)
        tangents = [(n.get('tangent_x', float('nan')),
                      n.get('tangent_y', float('nan')), 0.0)
                     for n in nodes]

        # Get mesh at this time step
        if reader_mode == 'vtkhdf':
            mesh = get_mesh_vtkhdf(reader, time_val)
        else:
            mesh = get_mesh_vtu(vtu_reader, vtu_files[time_step_idx])

        # Probe (extended)
        field_rows = probe_fields_extended(mesh, phys_coords,
                                             scalar_fields, vector_name,
                                             tangents)

        n_invalid = sum(1 for fr in field_rows if not fr['valid_point'])
        if n_invalid > 0:
            print(f'[warn] Frame {frame_str}: {n_invalid}/{len(nodes)} nodes '
                  'outside mesh domain (valid_point=0)')
            warn_count += 1

        for n, fr in zip(nodes, field_rows):
            row = {
                'case_id':               case_id,
                'time_step':             time_step_idx,
                'time':                  f'{time_val:.6g}',
                'swc_frame':             frame_str,
                'swc_node_id':           n['id'],
                'parent_id':             n['parent'],
                'node_type':             n['type'],
                'x_px':                  f"{n['x']:.3f}",
                'y_px':                  f"{n['y']:.3f}",
                'x_phys':                f"{n['x_phys']:.6f}",
                'y_phys':                f"{n['y_phys']:.6f}",
                'z_phys':                f"{n['z_phys']:.6f}",
                'distance_from_soma':    f"{n['distance_from_soma']:.6f}",
                'path_length_from_soma': f"{n['path_length_from_soma']:.6f}",
                'branch_order':          n['branch_order'],
                'is_tip':                n['is_tip'],
                'tangent_x':             ('' if n.get('tangent_x') != n.get('tangent_x')
                                           else f"{n.get('tangent_x', 0):.6f}"),
                'tangent_y':             ('' if n.get('tangent_y') != n.get('tangent_y')
                                           else f"{n.get('tangent_y', 0):.6f}"),
            }
            row.update(fr)
            out_rows.append(row)

            if dry_run:
                vals = '  '.join(
                    f'{fname}={fr.get(fname, "N/A")}' for fname in requested_fields
                )
                flag = '' if fr['valid_point'] else ' [INVALID]'
                print(f"  node {n['id']:3d}: "
                      f"phys=({n['x_phys']:.4f},{n['y_phys']:.4f})  "
                      f"tangent=({n.get('tangent_x', 0):.3f},"
                      f"{n.get('tangent_y', 0):.3f})  "
                      f"path={n['path_length_from_soma']:.4f}  {vals}{flag}")

    print(f'\n[info] Total rows assembled: {len(out_rows)}')
    if warn_count:
        print(f'[info] Warnings during processing: {warn_count}')

    if dry_run:
        print('\n[DRY RUN] No CSV written. Remove --dry-run to produce the full CSV.')
        return
    if not out_rows:
        sys.exit('[ERROR] No rows to write. Check SWC files and VTK time steps.')

    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    all_columns = (BASE_COLUMNS_EXT + scalar_fields
                   + ([vector_name] if vector_name else [])
                   + (['growth_vel_n'] if vector_name else [])
                   + ['valid_point'])
    with open(out_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=all_columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(out_rows)
    print(f'[done] Written {len(out_rows)} rows → {out_path}')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sim_group = parser.add_mutually_exclusive_group(required=True)
    sim_group.add_argument('--sim-file',
                            help='Path to Simulation.vtkhdf (preferred)')
    sim_group.add_argument('--vtu-dir',
                            help='Directory containing Output_t*.vtu files')
    parser.add_argument('--vtu-pattern', default='Output_t*.vtu')
    parser.add_argument('--swc-dir', required=True,
                         help='Directory of per-frame SWC subdirectories')
    parser.add_argument('--out-csv', default='swc_field_samples_ext.csv')
    parser.add_argument('--case-id', default='unknown')
    parser.add_argument('--fields', default=','.join(DEFAULT_FIELDS_EXT),
                         help='Comma-separated list of fields to sample '
                              '(includes norm_vel_int, curvature, GrowthVel '
                              'by default)')
    parser.add_argument('--swc-type', choices=['clean', 'pixel'],
                         default='clean')
    parser.add_argument('--image-size', type=int, default=2048)
    parser.add_argument('--smooth-sigma', type=float, default=1.5)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    process_case(args)


if __name__ == '__main__':
    main()
