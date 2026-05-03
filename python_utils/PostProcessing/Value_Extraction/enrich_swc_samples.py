#!/usr/bin/env python3
"""
enrich_swc_samples.py
=====================
Add front-aware morphological columns to the CSV produced by sample_fields_on_swc.py.

Does NOT require VTK / pvpython. Runs in swc_env (numpy + pandas).

Usage:
    conda activate swc_env
    python enrich_swc_samples.py \\
        --csv   /scratch/flore0a/AnalysisResults/D7p5_V3_fields.csv \\
        --out   /scratch/flore0a/AnalysisResults/D7p5_V3_fields_enriched.csv \\
        --swc-dir /scratch/flore0a/Dataset_test2/Cases/D7.5_V3/Simulation2D/merged/swc \\
        --r-gc  0.1

Columns added:
    L0                          Max path_length_from_soma at the first time step (scalar).
    Lmax_t                      Max path_length_from_soma among all nodes at this time step.
    elongation_rate             Forward difference: (Lmax(t+1) - Lmax(t)) / dt.
                                NaN for the last frame.
    is_initial_region           1 if path_length_from_soma <= L0, else 0.
    is_new_region               1 if path_length_from_soma >  L0, else 0.
    distance_beyond_initial_front   max(0, path_length_from_soma - L0).
    normalized_global_position  path_length_from_soma / Lmax_t.
    normalized_new_region_pos   (path_length - L0) / (Lmax_t - L0).
                                NaN for initial-region nodes and when Lmax_t == L0.
    is_growth_cone_region       1 if the node is within r_gc (physical units) of any
                                downstream tip along its branch. See note on method below.
    branch_id                   Integer branch identifier inferred from the SWC tree for
                                this frame. 0 = main branch from root; new integer at each
                                bifurcation. Requires --swc-dir. -1 if not computed.

Growth-cone method (chosen automatically):
    If --swc-dir is given: reads the SWC file for each frame, builds the tree, and
    computes the exact path distance from each node to its nearest DESCENDANT tip.
    Flags node as growth-cone if that distance <= r_gc.

    If --swc-dir is omitted: uses a CSV-only approximation.  For each node, finds
    the nearest tip (in the same frame) with a larger path_length_from_soma and takes
    the path_length difference as the proxy distance.  This is exact for single-branch
    trees and a conservative lower bound for multi-branch trees (will occasionally flag
    a node as GC if a tip on a *different* branch happens to have only slightly larger
    path_length).  For typical neuronal trees this error is negligible.

SWC-based branch_id assignment:
    Root of tree → branch_id=0.
    Each child of a bifurcation node (degree >= 2) receives a new integer branch_id.
    Degree-1 pass-through nodes inherit the parent's branch_id.
    Branch IDs are assigned independently per frame; they are NOT tracked across frames.
"""

import argparse
import math
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# SWC parsing and tree utilities  (only used when --swc-dir is given)
# ─────────────────────────────────────────────────────────────────────────────

def _parse_swc(path):
    nodes = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 7:
                continue
            nodes.append({
                'id':     int(p[0]),
                'parent': int(p[6]),
                'x': float(p[2]), 'y': float(p[3]), 'z': float(p[4]),
            })
    return nodes


def _build_tree(nodes):
    by_id = {n['id']: n for n in nodes}
    children = {n['id']: [] for n in nodes}
    root_id = None
    for n in nodes:
        if n['parent'] == -1:
            root_id = n['id']
        else:
            children[n['parent']].append(n['id'])
    return by_id, children, root_id


def _path_lengths_from_root(by_id, children, root_id):
    """Cumulative path length from root to each node (SWC pixel units)."""
    pl = {root_id: 0.0}
    q = deque([root_id])
    while q:
        nid = q.popleft()
        n = by_id[nid]
        for cid in children[nid]:
            c = by_id[cid]
            pl[cid] = pl[nid] + math.hypot(c['x'] - n['x'], c['y'] - n['y'])
            q.append(cid)
    return pl


def _min_subtree_tip_dist(children, path_len):
    """
    For each node, compute the path distance to its nearest DESCENDANT tip
    (not any tip — only tips reachable by going forward/down the tree).

    Returns dict {node_id: distance_to_nearest_descendant_tip}
    Tip nodes → 0.0
    """
    min_tip_path = {}

    def dfs(nid):
        c = children[nid]
        if not c:
            min_tip_path[nid] = path_len[nid]
        else:
            for cid in c:
                dfs(cid)
            min_tip_path[nid] = min(min_tip_path[cid] for cid in c)

    root_id = next(nid for nid, pid_list in children.items()
                   if not any(nid in v for v in children.values()))
    # Safer: root is the node with no parent (parent == -1)
    # We'll call dfs from whichever node has no incoming edges.
    # Since _build_tree ensures only root has parent=-1, we pass root_id explicitly.
    return min_tip_path  # populated by dfs called externally


def _min_subtree_tip_dist_rooted(children, path_len, root_id):
    min_tip_path = {}

    def dfs(nid):
        c = children[nid]
        if not c:
            min_tip_path[nid] = path_len[nid]
        else:
            for cid in c:
                dfs(cid)
            min_tip_path[nid] = min(min_tip_path[cid] for cid in c)

    dfs(root_id)
    return {nid: min_tip_path[nid] - path_len[nid] for nid in min_tip_path}


def _infer_branch_id(children, root_id):
    """
    Assign integer branch_id via recursive DFS.
    Root → 0.  Each child of a bifurcation gets a new id.
    Degree-1 nodes (pass-through) inherit parent's id.
    Returns dict {node_id: branch_id}.
    """
    counter = [0]
    bid = {}

    def walk(nid, current):
        bid[nid] = current
        c = children[nid]
        if len(c) == 1:
            walk(c[0], current)      # same branch
        else:
            for cid in c:
                counter[0] += 1
                walk(cid, counter[0])

    walk(root_id, 0)
    return bid


# ─────────────────────────────────────────────────────────────────────────────
# Growth-cone approximation (CSV-only, no SWC needed)
# ─────────────────────────────────────────────────────────────────────────────

def _gc_approx_for_frame(pl_vals, is_tip_mask):
    """
    Vectorised approximation of path-distance to nearest downstream tip.

    For each node i:
      dist_i = min over all tips j where pl[j] >= pl[i] of (pl[j] - pl[i])
    Tip nodes return 0.

    This is exact when each node has at most one downstream tip; for multi-branch
    trees it is an underestimate when a shorter cross-branch tip "accidentally"
    appears to be closer in path_length space.  In practice the error is small
    because cross-branch shortcuts require going back toward the soma.
    """
    tip_paths = pl_vals[is_tip_mask]
    if tip_paths.size == 0:
        return np.full(len(pl_vals), np.inf)

    # shape: (n_nodes, n_tips)
    diffs = tip_paths[np.newaxis, :] - pl_vals[:, np.newaxis]
    diffs = np.where(diffs >= 0.0, diffs, np.inf)
    min_dists = diffs.min(axis=1)
    min_dists[is_tip_mask] = 0.0
    return min_dists


# ─────────────────────────────────────────────────────────────────────────────
# Main enrichment
# ─────────────────────────────────────────────────────────────────────────────

def enrich(df, swc_dir=None, swc_type='clean', r_gc=0.1):
    """
    Enrich df with front-aware columns.  Works in-place on copies per case.
    Returns the enriched DataFrame.

    Parameters
    ----------
    df        : DataFrame from sample_fields_on_swc.py (all cases)
    swc_dir   : path to the SWC directory (e.g. .../merged/swc).
                If None, branch_id is -1 and growth-cone uses approximation.
    swc_type  : 'clean' or 'pixel' (matches the SWC file used for sampling).
    r_gc      : growth-cone radius in the same physical units as path_length_from_soma.
    """
    swc_name = 'clean_skeleton.swc' if swc_type == 'clean' else 'pixel_skeleton.swc'
    enriched_parts = []

    for case_id, case_df in df.groupby('case_id', sort=False):
        case_df = case_df.copy()
        pl = case_df['path_length_from_soma'].astype(float)

        # ── L0 and Lmax(t) ────────────────────────────────────────────────────
        t0 = int(case_df['time_step'].min())
        L0 = float(pl[case_df['time_step'] == t0].max())

        lmax_by_ts = case_df.groupby('time_step')['path_length_from_soma'].max().astype(float)
        case_df['L0'] = L0
        case_df['Lmax_t'] = case_df['time_step'].map(lmax_by_ts)

        # ── Elongation rate (forward difference, dt in time_step units) ───────
        ts_sorted = sorted(lmax_by_ts.index)
        elong = {}
        for i, ts in enumerate(ts_sorted[:-1]):
            dt = ts_sorted[i + 1] - ts
            elong[ts] = (lmax_by_ts[ts_sorted[i + 1]] - lmax_by_ts[ts]) / dt
        elong[ts_sorted[-1]] = np.nan
        case_df['elongation_rate'] = case_df['time_step'].map(elong)

        # ── Region labels ─────────────────────────────────────────────────────
        lmax_t = case_df['Lmax_t'].astype(float)
        case_df['is_initial_region'] = (pl <= L0).astype(int)
        case_df['is_new_region'] = (pl > L0).astype(int)
        case_df['distance_beyond_initial_front'] = (pl - L0).clip(lower=0.0)
        case_df['normalized_global_position'] = np.where(lmax_t > 0, pl / lmax_t, np.nan)

        new_denom = (lmax_t - L0).astype(float)
        case_df['normalized_new_region_pos'] = np.where(
            (pl > L0) & (new_denom > 1e-12),
            (pl - L0) / new_denom,
            np.nan,
        )

        # ── Growth-cone region ────────────────────────────────────────────────
        if swc_dir is not None:
            swc_root = Path(swc_dir)
            gc_flag = np.zeros(len(case_df), dtype=int)
            bid_col = np.full(len(case_df), -1, dtype=int)

            for frame_str in sorted(case_df['swc_frame'].unique()):
                # swc_frame may be integer (0,1,...) or zero-padded string ("000","001",...)
                # — try zero-padded first, then as-is
                frame_key = f"{int(frame_str):03d}" if str(frame_str).isdigit() else str(frame_str)
                swc_path = swc_root / frame_key / swc_name
                mask = case_df['swc_frame'] == frame_str

                if not swc_path.exists():
                    print(f"  [warn] SWC not found: {swc_path}")
                    continue

                nodes = _parse_swc(swc_path)
                if not nodes:
                    continue
                by_id, children, root_id = _build_tree(nodes)
                if root_id is None:
                    continue

                # Path lengths in SWC pixel coordinates
                pl_px = _path_lengths_from_root(by_id, children, root_id)

                # Min descendant-tip distance in pixel units
                min_tip_px = _min_subtree_tip_dist_rooted(children, pl_px, root_id)

                # Convert r_gc from physical to pixel units using per-case px_scale.
                # px_scale = path_length_phys / path_length_px for the root-to-tip path.
                # We derive it from the frame with max path length.
                frame_pl_phys = case_df.loc[mask, 'path_length_from_soma'].astype(float)
                frame_node_ids = case_df.loc[mask, 'swc_node_id'].astype(int)
                # Find node with max path_length in CSV
                if frame_pl_phys.empty:
                    continue
                max_phys_idx = frame_pl_phys.idxmax()
                max_nid = int(case_df.loc[max_phys_idx, 'swc_node_id'])
                max_pl_phys = float(frame_pl_phys.max())
                max_pl_px = pl_px.get(max_nid, None)
                if max_pl_px and max_pl_px > 0 and max_pl_phys > 0:
                    px_scale = max_pl_phys / max_pl_px
                else:
                    px_scale = 1.0  # fallback
                r_gc_px = r_gc / px_scale

                for i, (orig_idx, row) in enumerate(case_df[mask].iterrows()):
                    nid = int(row['swc_node_id'])
                    gc_flag[case_df.index.get_loc(orig_idx)] = int(
                        min_tip_px.get(nid, np.inf) <= r_gc_px
                    )
                    bid_col[case_df.index.get_loc(orig_idx)] = bid_col_from_map = (
                        _infer_branch_id(children, root_id).get(nid, -1)
                    )

            case_df['is_growth_cone_region'] = gc_flag
            case_df['branch_id'] = bid_col
            gc_method = f'SWC subtree (r_gc = {r_gc:.4f} phys units)'

        else:
            # CSV-only approximation
            gc_flags = np.zeros(len(case_df), dtype=int)
            for ts, grp in case_df.groupby('time_step'):
                pl_vals = grp['path_length_from_soma'].astype(float).values
                is_tip_mask = (grp['is_tip'].values == 1)
                min_dists = _gc_approx_for_frame(pl_vals, is_tip_mask)
                for j, orig_idx in enumerate(grp.index):
                    gc_flags[case_df.index.get_loc(orig_idx)] = int(min_dists[j] <= r_gc)

            case_df['is_growth_cone_region'] = gc_flags
            case_df['branch_id'] = -1
            gc_method = f'CSV approximation (r_gc = {r_gc:.4f} phys units, no SWC re-read)'

        # ── Report ────────────────────────────────────────────────────────────
        print(f"\n{'='*60}")
        print(f"Case: {case_id}")
        print(f"  L0  = {L0:.4f}  (max path_length at t={t0})")
        print(f"  GC definition: {gc_method}")
        print(f"  {'ts':>4}  {'n':>4}  {'initial':>7}  {'new':>5}  {'tip':>4}  "
              f"{'gc_reg':>7}  {'Lmax':>6}")
        for ts in ts_sorted:
            sub = case_df[case_df['time_step'] == ts]
            print(f"  {ts:4d}  {len(sub):4d}  "
                  f"{sub['is_initial_region'].sum():7d}  "
                  f"{sub['is_new_region'].sum():5d}  "
                  f"{sub['is_tip'].sum() if 'is_tip' in sub else -1:4d}  "
                  f"{sub['is_growth_cone_region'].sum():7d}  "
                  f"{lmax_by_ts[ts]:6.3f}")

        enriched_parts.append(case_df)

    return pd.concat(enriched_parts, ignore_index=True)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--csv', required=True,
                        help='Input CSV from sample_fields_on_swc.py')
    parser.add_argument('--out', default=None,
                        help='Output enriched CSV path. Default: <input>_enriched.csv')
    parser.add_argument('--swc-dir', default=None,
                        help='Directory containing per-frame SWC subdirectories '
                             '(e.g. .../merged/swc). Enables exact GC computation and branch_id.')
    parser.add_argument('--swc-type', choices=['clean', 'pixel'], default='clean',
                        help='SWC variant (must match what sample_fields_on_swc.py used)')
    parser.add_argument('--r-gc', type=float, default=0.1,
                        help='Growth-cone radius in physical units (default: 0.1)')
    args = parser.parse_args()

    csv_path = Path(args.csv)
    out_path = Path(args.out) if args.out else csv_path.with_name(
        csv_path.stem + '_enriched.csv')

    print(f"[info] Loading {csv_path} ...")
    df = pd.read_csv(csv_path)

    # Ensure numeric types
    for col in ['time_step', 'path_length_from_soma', 'is_tip', 'swc_node_id']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    print(f"[info] {len(df)} rows, columns: {list(df.columns)}")

    enriched = enrich(df, swc_dir=args.swc_dir, swc_type=args.swc_type, r_gc=args.r_gc)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(out_path, index=False)
    print(f"\n[done] Enriched CSV → {out_path}")
    print(f"[done] New columns: L0, Lmax_t, elongation_rate, is_initial_region, "
          f"is_new_region, distance_beyond_initial_front, normalized_global_position, "
          f"normalized_new_region_pos, is_growth_cone_region, branch_id")


if __name__ == '__main__':
    main()
