#!/usr/bin/env pvpython
"""
sample_fields_on_swc.py
=======================
Sample simulation field values at SWC skeleton node coordinates and write
a long-format CSV for quantitative morphological analysis.

Requires pvpython (VTK is bundled with ParaView):
    module load paraview/6.1.0-mesa
    pvpython sample_fields_on_swc.py [options]

Also works with any Python that has vtk >= 9.1 installed.

────────────────────────────────────────────────────────────────────────────
USAGE EXAMPLES
────────────────────────────────────────────────────────────────────────────

Dry-run / inspect (no CSV written, samples 10 nodes from first SWC frame):
    pvpython sample_fields_on_swc.py \\
        --sim-file /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/Simulation.vtkhdf \\
        --swc-dir  /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc \\
        --out-csv  /scratch/flore0a/AnalysisResults/D7p5_V3_fields.csv \\
        --case-id  D7p5_V3 \\
        --dry-run

Full run (VTKHDF, preferred):
    pvpython sample_fields_on_swc.py \\
        --sim-file /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/Simulation.vtkhdf \\
        --swc-dir  /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc \\
        --out-csv  /scratch/flore0a/AnalysisResults/D7p5_V3_fields.csv \\
        --case-id  D7p5_V3

VTU series fallback:
    pvpython sample_fields_on_swc.py \\
        --vtu-dir /scratch/.../Cases/D7.5_V3/Simulation2D/merged \\
        --vtu-pattern "Output_t*.vtu" \\
        --swc-dir  /scratch/.../merged/swc \\
        --out-csv  output.csv \\
        --case-id  D7p5_V3

────────────────────────────────────────────────────────────────────────────
BATCH EXAMPLE (bash loop over Dataset_test2):
────────────────────────────────────────────────────────────────────────────

    CASES_ROOT=/scratch/flore0a/Dataset_test2/Cases
    OUT_ROOT=/scratch/flore0a/AnalysisResults/field_samples

    module load paraview/6.1.0-mesa
    for case_dir in $CASES_ROOT/*/; do
        case_id=$(basename "$case_dir")
        sim_file="$case_dir/Simulation2D/Simulation.vtkhdf"
        swc_dir="$case_dir/Simulation2D/merged/swc"
        out_csv="$OUT_ROOT/${case_id}_fields.csv"
        if [ -f "$sim_file" ] && [ -d "$swc_dir" ]; then
            pvpython sample_fields_on_swc.py \\
                --sim-file "$sim_file" \\
                --swc-dir  "$swc_dir" \\
                --out-csv  "$out_csv" \\
                --case-id  "$case_id"
        fi
    done

────────────────────────────────────────────────────────────────────────────
COORDINATE TRANSFORM (see FIELD_SAMPLING_REPORT.md for derivation)
────────────────────────────────────────────────────────────────────────────

SWC coordinates are pixel coordinates from the 2048×2048 rendered PNG:
    x_swc = pixel column (0 = left)
    y_swc = -(pixel row)  (0 = top, so y_swc is negative for most nodes)

The transform to physical simulation coordinates is reconstructed from the
VTK domain bounds at t=0 (same formula used by the PNG export pipeline):

    S         = max(domain_width, domain_height) / 2
    x_phys    = center_x + (x_swc / 2048 - 0.5) * 2 * S
    y_phys    = center_y + (0.5 + y_swc / 2048) * 2 * S
    scale     = 2 * S / 2048  [physical units per pixel, uniform]

Verified on D7.5_V3: soma node (x_swc=1032, y_swc=-1907) → (0.023, -1.437),
lsf = -0.018 (inside cell ✓), u_t = 0.249, u_ca_cyt = 0.251 at t=0.
"""

import argparse
import csv
import glob
import math
import os
import sys
from collections import deque
from pathlib import Path

import vtk


# ─────────────────────────────────────────────────────────────────────────────
# SWC parsing and topology
# ─────────────────────────────────────────────────────────────────────────────

def parse_swc(swc_path):
    """Parse SWC file, return list of node dicts."""
    nodes = []
    with open(swc_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            nodes.append({
                'id':     int(parts[0]),
                'type':   int(parts[1]),
                'x':      float(parts[2]),  # x_swc = pixel column
                'y':      float(parts[3]),  # y_swc = -(pixel row)
                'z':      float(parts[4]),
                'radius': float(parts[5]),
                'parent': int(parts[6]),
            })
    return nodes


def compute_topology(nodes):
    """
    Add topological metrics to each node (in-place):
        path_length_px   — cumulative edge length along tree path from soma (pixels)
        distance_px      — Euclidean distance from soma (pixels)
        branch_order     — number of bifurcation points on path from soma
        is_tip           — 1 if leaf node (no children), else 0
    """
    by_id = {n['id']: n for n in nodes}
    children = {n['id']: [] for n in nodes}
    root_id = None

    for n in nodes:
        if n['parent'] == -1:
            root_id = n['id']
        else:
            if n['parent'] in children:
                children[n['parent']].append(n['id'])

    if root_id is None:
        # Fallback: treat first node as root
        root_id = nodes[0]['id']

    root = by_id[root_id]
    path_len = {root_id: 0.0}
    branch_ord = {root_id: 0}
    queue = deque([root_id])

    while queue:
        nid = queue.popleft()
        n = by_id[nid]
        is_bifurcation = len(children[nid]) >= 2
        for cid in children[nid]:
            c = by_id[cid]
            edge = math.hypot(c['x'] - n['x'], c['y'] - n['y'])
            path_len[cid] = path_len[nid] + edge
            branch_ord[cid] = branch_ord[nid] + (1 if is_bifurcation else 0)
            queue.append(cid)

    for n in nodes:
        nid = n['id']
        n['path_length_px'] = path_len.get(nid, 0.0)
        n['distance_px'] = math.hypot(n['x'] - root['x'], n['y'] - root['y'])
        n['branch_order'] = branch_ord.get(nid, 0)
        n['is_tip'] = int(len(children.get(nid, [])) == 0 and nid != root_id)

    return nodes


# ─────────────────────────────────────────────────────────────────────────────
# Pixel-skeleton coordinate smoothing
# ─────────────────────────────────────────────────────────────────────────────

def smooth_pixel_swc(nodes, sigma=1.5):
    """
    Gaussian-smooth (x, y) pixel coordinates along each tree branch to remove
    the single-pixel staircase oscillations typical of pixel-wise skeletons.

    The tree is decomposed into linear segments (root/junction → next
    junction/tip).  A 1-D Gaussian with the given sigma (pixels) is applied
    to x and y independently along each segment.  The global root node is
    always pinned to its original position.

    Implemented in pure Python (no numpy) so it runs inside pvpython.
    """
    if len(nodes) <= 2:
        return nodes

    by_id = {n['id']: n for n in nodes}
    children = {n['id']: [] for n in nodes}
    root_id = None
    for n in nodes:
        if n['parent'] == -1:
            root_id = n['id']
        elif n['parent'] in children:
            children[n['parent']].append(n['id'])
    if root_id is None:
        return nodes

    # Build Gaussian kernel (pure Python, reflect-pad for boundaries)
    ksize = max(3, int(math.ceil(6 * sigma)) | 1)   # odd integer
    half  = ksize // 2
    raw_k = [math.exp(-0.5 * (i - half) ** 2 / sigma ** 2) for i in range(ksize)]
    total = sum(raw_k)
    kernel = [v / total for v in raw_k]

    def _gauss1d(arr):
        n = len(arr)
        if n < 3:
            return list(arr)
        # reflect-pad: left = arr[half:0:-1], right = arr[-2:-2-half:-1]
        left  = [arr[min(half - j, n - 1)] for j in range(half, 0, -1)]
        right = [arr[max(n - 2 - j, 0)]   for j in range(half)]
        padded = left + list(arr) + right
        # 1-D convolution (valid mode)
        out = []
        for i in range(n):
            val = sum(padded[i + j] * kernel[j] for j in range(ksize))
            out.append(val)
        return out

    # DFS: extract linear segments between junctions/tips
    def _segments(nid, seg):
        seg.append(nid)
        c = children[nid]
        if len(c) == 0:
            yield list(seg)
        elif len(c) == 1:
            yield from _segments(c[0], seg)
        else:
            yield list(seg)              # segment ends at junction
            for cid in c:
                yield from _segments(cid, [nid])   # next segment from junction

    for seg in _segments(root_id, []):
        if len(seg) < 3:
            continue
        xs = [by_id[nid]['x'] for nid in seg]
        ys = [by_id[nid]['y'] for nid in seg]
        xs_s = _gauss1d(xs)
        ys_s = _gauss1d(ys)
        for i, nid in enumerate(seg):
            if nid == root_id:
                continue    # pin soma
            by_id[nid]['x'] = xs_s[i]
            by_id[nid]['y'] = ys_s[i]

    return nodes


# ─────────────────────────────────────────────────────────────────────────────
# Coordinate transform
# ─────────────────────────────────────────────────────────────────────────────

def build_pixel_to_physical(t0_bounds, image_size=2048):
    """
    Reconstruct pixel→physical coordinate transform from VTK domain bounds at t=0.

    Uses the same formula as util_postanalisis_export_png.py:
        - camera centered at domain center
        - parallel scale S = max(domain_width, domain_height) / 2
        - square image W = H = image_size

    Returns (transform_fn, S, px_scale) where:
        transform_fn(x_swc, y_swc) -> (x_phys, y_phys)
        S          = parallel scale (half-width of the view in physical units)
        px_scale   = physical units per pixel (isotropic)
    """
    xmin, xmax, ymin, ymax = t0_bounds[0], t0_bounds[1], t0_bounds[2], t0_bounds[3]
    cx = (xmin + xmax) / 2.0
    cy = (ymin + ymax) / 2.0
    dw = xmax - xmin
    dh = ymax - ymin
    S = max(dw, dh) / 2.0
    W = H = float(image_size)
    px_scale = 2.0 * S / W  # physical units per pixel

    def transform(x_swc, y_swc):
        x_phys = cx + (x_swc / W - 0.5) * 2.0 * S
        y_phys = cy + (0.5 + y_swc / H) * 2.0 * S
        return x_phys, y_phys

    return transform, S, px_scale


# ─────────────────────────────────────────────────────────────────────────────
# VTK readers
# ─────────────────────────────────────────────────────────────────────────────

def open_vtkhdf(path):
    """
    Open a VTKHDF time-series file.
    Returns (reader, sorted_times, t0_bounds, point_field_names).
    """
    reader = vtk.vtkHDFReader()
    reader.SetFileName(str(path))
    reader.UpdateInformation()

    info = reader.GetOutputInformation(0)
    TS_KEY = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()
    n = info.Length(TS_KEY)
    if n == 0:
        sys.exit(f"[ERROR] No time steps found in {path}")
    times = sorted([info.Get(TS_KEY, i) for i in range(n)])

    # Update at t=0 to get bounds and field names
    info.Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), times[0])
    reader.Update()
    mesh = reader.GetOutput()
    t0_bounds = mesh.GetBounds()

    field_names = _get_point_field_names(mesh)
    print(f"[info] VTKHDF opened: {n} time steps, "
          f"t0_bounds = {[round(x, 4) for x in t0_bounds]}")
    print(f"[info] Point arrays: {field_names}")

    return reader, times, t0_bounds, field_names


def open_vtu_series(vtu_dir, pattern="Output_t*.vtu"):
    """
    Open a VTU file series.
    Returns (vtu_files, reader, t0_bounds, point_field_names).
    """
    vtu_files = sorted(glob.glob(os.path.join(str(vtu_dir), pattern)))
    if not vtu_files:
        sys.exit(f"[ERROR] No VTU files matching '{pattern}' in {vtu_dir}")

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_files[0])
    reader.Update()
    mesh = reader.GetOutput()
    t0_bounds = mesh.GetBounds()

    field_names = _get_point_field_names(mesh)
    print(f"[info] VTU series: {len(vtu_files)} files, "
          f"t0_bounds = {[round(x, 4) for x in t0_bounds]}")
    print(f"[info] Point arrays: {field_names}")

    return vtu_files, reader, t0_bounds, field_names


def _get_point_field_names(ugrid):
    pd = ugrid.GetPointData()
    return [pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())]


def get_mesh_vtkhdf(reader, time_val):
    """Update VTKHDF reader to the requested time, return the mesh."""
    info = reader.GetOutputInformation(0)
    info.Set(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP(), time_val)
    reader.Update()
    return reader.GetOutput()


def get_mesh_vtu(reader, vtu_file):
    """Read a single VTU file, return the mesh."""
    reader.SetFileName(str(vtu_file))
    reader.Update()
    return reader.GetOutput()


# ─────────────────────────────────────────────────────────────────────────────
# Field probing
# ─────────────────────────────────────────────────────────────────────────────

def probe_fields(mesh, phys_coords, requested_fields):
    """
    Sample field values at physical coordinates from an unstructured mesh.

    phys_coords: list of (x, y, z) tuples
    requested_fields: list of field names to extract

    Returns list of dicts, one per point:
        {field_name: value_str, ..., 'valid_point': 1|0}
    Invalid points (outside mesh domain) get empty strings for field values.
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
    probe.Update()

    out_pd = probe.GetOutput().GetPointData()
    valid_arr = out_pd.GetArray("vtkValidPointMask")

    rows = []
    n = len(phys_coords)
    for i in range(n):
        valid = bool(valid_arr.GetValue(i)) if valid_arr else True
        row = {'valid_point': int(valid)}
        for fname in requested_fields:
            arr = out_pd.GetArray(fname)
            if arr is None:
                row[fname] = ''
            elif valid:
                row[fname] = f"{arr.GetValue(i):.8g}"
            else:
                row[fname] = ''  # outside domain — flagged by valid_point=0
        rows.append(row)
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# SWC frame discovery
# ─────────────────────────────────────────────────────────────────────────────

def discover_swc_frames(swc_dir, swc_type="clean"):
    """
    Find all per-frame SWC files in swc_dir.
    Expects structure: swc_dir/<frame_idx>/{clean,pixel}_skeleton.swc

    Returns sorted list of (frame_idx_str, swc_path).
    """
    swc_name = "clean_skeleton.swc" if swc_type == "clean" else "pixel_skeleton.swc"
    frames = []
    swc_path_obj = Path(swc_dir)
    if not swc_path_obj.is_dir():
        sys.exit(f"[ERROR] swc-dir does not exist: {swc_dir}")

    for entry in sorted(swc_path_obj.iterdir()):
        if entry.is_dir():
            candidate = entry / swc_name
            if candidate.exists():
                frames.append((entry.name, candidate))

    return frames


# ─────────────────────────────────────────────────────────────────────────────
# CSV column definitions
# ─────────────────────────────────────────────────────────────────────────────

BASE_COLUMNS = [
    'case_id', 'time_step', 'time', 'swc_frame',
    'swc_node_id', 'parent_id', 'node_type',
    'x_px', 'y_px',
    'x_phys', 'y_phys', 'z_phys',
    'distance_from_soma', 'path_length_from_soma', 'branch_order', 'is_tip',
]

DEFAULT_FIELDS = ['u_t', 'u_u', 'u_b', 'u_p', 'u_ca_cyt', 'inhibitor', 'lsf']


# ─────────────────────────────────────────────────────────────────────────────
# Main processing
# ─────────────────────────────────────────────────────────────────────────────

def process_case(args):
    case_id = args.case_id
    dry_run = args.dry_run
    swc_type = args.swc_type
    image_size = args.image_size
    requested_fields = [f.strip() for f in args.fields.split(',') if f.strip()]

    # ── Open simulation reader ────────────────────────────────────────────────
    if args.sim_file:
        reader, times, t0_bounds, available_fields = open_vtkhdf(args.sim_file)
        reader_mode = "vtkhdf"
        vtu_files = None
        vtu_reader = None
    else:
        vtu_files, vtu_reader, t0_bounds, available_fields = open_vtu_series(
            args.vtu_dir, args.vtu_pattern)
        times = list(range(len(vtu_files)))
        reader_mode = "vtu"
        reader = None

    # Warn about requested fields not present in the file
    missing = [f for f in requested_fields if f not in available_fields]
    if missing:
        print(f"[warn] Requested fields not found in simulation file: {missing}")
        print(f"[warn] Available: {available_fields}")

    # ── Coordinate transform ──────────────────────────────────────────────────
    transform, S, px_scale = build_pixel_to_physical(t0_bounds, image_size)
    print(f"[info] Pixel scale: {px_scale:.6f} phys_units/pixel  (S = {S:.4f})")

    # ── Discover SWC frames ───────────────────────────────────────────────────
    all_frames = discover_swc_frames(args.swc_dir, swc_type)
    print(f"[info] Found {len(all_frames)} SWC frames (type: {swc_type})")
    if not all_frames:
        sys.exit("[ERROR] No SWC frames found.")

    frames_to_process = all_frames[:1] if dry_run else all_frames

    # ── Process each frame ────────────────────────────────────────────────────
    out_rows = []
    warn_count = 0

    for frame_str, swc_path in frames_to_process:
        time_step_idx = int(frame_str)
        if time_step_idx >= len(times):
            print(f"[warn] Frame {frame_str}: index {time_step_idx} out of range "
                  f"({len(times)} time steps), skipping.")
            warn_count += 1
            continue
        time_val = times[time_step_idx]

        # Load SWC
        nodes = parse_swc(swc_path)
        if not nodes:
            print(f"[warn] Frame {frame_str}: empty SWC, skipping.")
            warn_count += 1
            continue
        if swc_type == 'pixel' and args.smooth_sigma > 0:
            nodes = smooth_pixel_swc(nodes, sigma=args.smooth_sigma)
        nodes = compute_topology(nodes)

        if dry_run:
            nodes = nodes[:10]
            print(f"\n[DRY RUN] Frame {frame_str}, time={time_val:.4g}, "
                  f"sampling {len(nodes)} nodes:")

        # Transform pixel → physical
        phys_coords = []
        for n in nodes:
            xp, yp = transform(n['x'], n['y'])
            n['x_phys'], n['y_phys'], n['z_phys'] = xp, yp, 0.0
            n['distance_from_soma'] = n['distance_px'] * px_scale
            n['path_length_from_soma'] = n['path_length_px'] * px_scale
            phys_coords.append((xp, yp, 0.0))

        # Get mesh at this time step
        if reader_mode == "vtkhdf":
            mesh = get_mesh_vtkhdf(reader, time_val)
        else:
            mesh = get_mesh_vtu(vtu_reader, vtu_files[time_step_idx])

        # Probe fields
        field_rows = probe_fields(mesh, phys_coords, requested_fields)

        # Check for invalid points
        n_invalid = sum(1 for fr in field_rows if not fr['valid_point'])
        if n_invalid > 0:
            print(f"[warn] Frame {frame_str}: {n_invalid}/{len(nodes)} nodes "
                  "outside mesh domain (valid_point=0)")
            warn_count += 1

        # Assemble output rows
        for n, fr in zip(nodes, field_rows):
            row = {
                'case_id':               case_id,
                'time_step':             time_step_idx,
                'time':                  f"{time_val:.6g}",
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
            }
            row.update(fr)
            out_rows.append(row)

            if dry_run:
                vals = "  ".join(
                    f"{fname}={fr.get(fname, 'N/A')}"
                    for fname in requested_fields
                )
                flag = "" if fr['valid_point'] else " [INVALID]"
                print(f"  node {n['id']:3d}: "
                      f"px=({n['x']:.0f},{n['y']:.0f})  "
                      f"phys=({n['x_phys']:.4f},{n['y_phys']:.4f})  "
                      f"path={n['path_length_from_soma']:.4f}  "
                      f"bo={n['branch_order']}  tip={n['is_tip']}  "
                      f"{vals}{flag}")

    print(f"\n[info] Total rows assembled: {len(out_rows)}")
    if warn_count:
        print(f"[info] Warnings during processing: {warn_count}")

    if dry_run:
        print("\n[DRY RUN] No CSV written. Remove --dry-run to produce the full CSV.")
        return

    if not out_rows:
        sys.exit("[ERROR] No rows to write. Check that SWC files and VTK time steps align.")

    # ── Write CSV ─────────────────────────────────────────────────────────────
    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    all_columns = BASE_COLUMNS + requested_fields + ['valid_point']

    with open(out_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=all_columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"[done] Written {len(out_rows)} rows → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    sim_group = parser.add_mutually_exclusive_group(required=True)
    sim_group.add_argument(
        '--sim-file',
        help='Path to Simulation.vtkhdf (preferred)',
    )
    sim_group.add_argument(
        '--vtu-dir',
        help='Directory containing Output_t*.vtu files (fallback)',
    )

    parser.add_argument(
        '--vtu-pattern', default='Output_t*.vtu',
        help='Glob pattern for VTU files when using --vtu-dir (default: Output_t*.vtu)',
    )
    parser.add_argument(
        '--swc-dir', required=True,
        help='Directory containing per-frame SWC subdirectories '
             '(e.g. .../merged/swc with sub-dirs 000/, 001/, ...)',
    )
    parser.add_argument(
        '--out-csv', default='swc_field_samples.csv',
        help='Output CSV file path (default: swc_field_samples.csv)',
    )
    parser.add_argument(
        '--case-id', default='unknown',
        help='Case identifier written to the case_id CSV column (e.g. D7p5_V3)',
    )
    parser.add_argument(
        '--fields', default=','.join(DEFAULT_FIELDS),
        help=f'Comma-separated list of simulation fields to sample '
             f'(default: {",".join(DEFAULT_FIELDS)})',
    )
    parser.add_argument(
        '--swc-type', choices=['clean', 'pixel'], default='clean',
        help='SWC variant to use: clean (junction nodes only, default) '
             'or pixel (every skeleton pixel)',
    )
    parser.add_argument(
        '--image-size', type=int, default=2048,
        help='Side length of the square PNG used for SWC extraction '
             '(default: 2048, matches the rendering pipeline)',
    )
    parser.add_argument(
        '--smooth-sigma', type=float, default=1.5,
        help='Gaussian smoothing sigma (pixels) applied to pixel-skeleton '
             'coordinates before sampling (default: 1.5). '
             'Set to 0 to disable. Only active when --swc-type pixel.',
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='Sample only 10 nodes from the first SWC frame, '
             'print results, and exit without writing CSV',
    )

    args = parser.parse_args()
    process_case(args)


if __name__ == '__main__':
    main()
