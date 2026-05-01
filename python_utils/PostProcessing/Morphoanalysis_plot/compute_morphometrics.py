#!/usr/bin/env python3
"""
Phase 2: Scan all clean_skeleton.swc files, compute morphology metrics,
and write one aggregated CSV dataset.

Expected path pattern:
  {base_path}/{case}/Simulation2D/merged/swc/{timestep}/clean_skeleton.swc

Usage:
  python3 compute_morphometrics.py --path /scratch/flore0a/Dataset_test2/Cases --out morphometrics.csv
"""
import argparse
import csv
import glob
import math
import os
import statistics
from collections import deque

try:
    import numpy as np
    from scipy.spatial import ConvexHull
    _SCIPY = True
except ImportError:
    _SCIPY = False


def parse_swc(swc_path):
    """Parse SWC file. Returns {node_id: (x, y, z, parent_id)}."""
    nodes = {}
    with open(swc_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            nid = int(parts[0])
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            pid = int(parts[6])
            nodes[nid] = (x, y, z, pid)
    return nodes


def _asymmetry_metrics(nodes, children, root_id, root_x, bfs_order):
    """
    Compute left/right asymmetry and topological subtree asymmetry.

    Left/right split: vertical line x = root_x.
      Each non-root node is assigned 'left' (x < root_x) or 'right' (x > root_x).
      Nodes on the axis (x == root_x) are excluded from the L/R totals.

    Subtree asymmetry at binary branch points:
      |n1 - n2| / (n1 + n2)  where n1, n2 = tip counts in each child subtree.
      Computed via reverse-BFS (post-order) to avoid recursion.
    """
    # Left / right asymmetry
    L_length = R_length = 0.0
    L_tips = R_tips = 0

    for nid, (x, y, z, pid) in nodes.items():
        if nid == root_id:
            continue
        if x < root_x:
            side = "L"
        elif x > root_x:
            side = "R"
        else:
            continue                       # on the vertical axis — exclude

        if pid in nodes:
            px, py, pz, _ = nodes[pid]
            edge_len = math.sqrt((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2)
            if side == "L":
                L_length += edge_len
            else:
                R_length += edge_len

        if len(children[nid]) == 0:        # tip
            if side == "L":
                L_tips += 1
            else:
                R_tips += 1

    total_lr_len  = L_length + R_length
    total_lr_tips = L_tips + R_tips

    lr_asym_length = (
        abs(L_length - R_length) / total_lr_len if total_lr_len > 0.0 else float("nan")
    )
    lr_asym_tips = (
        abs(L_tips - R_tips) / total_lr_tips if total_lr_tips > 0 else float("nan")
    )

    # Subtree tip counts via reverse-BFS (post-order)
    subtree_tips = {}
    for nid in reversed(bfs_order):
        if len(children[nid]) == 0:
            subtree_tips[nid] = 1
        else:
            subtree_tips[nid] = sum(subtree_tips[c] for c in children[nid])

    # Partition asymmetry at each binary branch point
    asym_values = []
    for nid in nodes:
        if len(children[nid]) == 2:
            n1 = subtree_tips[children[nid][0]]
            n2 = subtree_tips[children[nid][1]]
            asym_values.append(abs(n1 - n2) / (n1 + n2))
        # branch points with >2 children are skipped

    subtree_asym_mean = (
        round(statistics.mean(asym_values), 4) if asym_values else float("nan")
    )

    return {
        "left_right_asymmetry_length": round(lr_asym_length, 4)
                                       if not math.isnan(lr_asym_length) else float("nan"),
        "left_right_asymmetry_tips":   round(lr_asym_tips, 4)
                                       if not math.isnan(lr_asym_tips) else float("nan"),
        "subtree_asymmetry_mean":      subtree_asym_mean,
    }


def _branching_geometry_metrics(nodes, children):
    """
    Branch angles and section-length statistics from the clean SWC.

    In the clean SWC every edge is a morphological section (no degree-2 nodes).
    Terminal sections are edges whose child node is a tip (degree 0).

    Branch angle at a binary branch point:
        v1 = child1 - node,  v2 = child2 - node
        angle = arccos( clamp(v1·v2 / (|v1|·|v2|), -1, 1) )  [degrees]
    Multifurcations (>2 children) are skipped.
    """
    tip_set = {nid for nid in nodes if len(children[nid]) == 0}

    # Collect section lengths (all edges) and terminal section lengths
    section_lengths = []
    terminal_lengths = []
    for nid, (x, y, z, pid) in nodes.items():
        if pid == -1 or pid not in nodes:
            continue
        px, py, pz, _ = nodes[pid]
        elen = math.sqrt((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2)
        section_lengths.append(elen)
        if nid in tip_set:
            terminal_lengths.append(elen)

    term_mean = (
        round(sum(terminal_lengths) / len(terminal_lengths), 2)
        if terminal_lengths else float("nan")
    )

    if section_lengths:
        sl_mean = sum(section_lengths) / len(section_lengths)
        if sl_mean > 0.0:
            sl_std = math.sqrt(
                sum((s - sl_mean) ** 2 for s in section_lengths) / len(section_lengths)
            )
            section_cv = round(sl_std / sl_mean, 4)
        else:
            section_cv = float("nan")
    else:
        section_cv = float("nan")

    # Branch angles at binary branch points
    branch_angles = []
    for nid in nodes:
        if len(children[nid]) != 2:
            continue
        nx, ny, nz, _ = nodes[nid]
        c1, c2 = children[nid]
        c1x, c1y, c1z, _ = nodes[c1]
        c2x, c2y, c2z, _ = nodes[c2]
        v1 = (c1x - nx, c1y - ny, c1z - nz)
        v2 = (c2x - nx, c2y - ny, c2z - nz)
        len1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
        len2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
        if len1 == 0.0 or len2 == 0.0:
            continue
        cos_a = (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) / (len1 * len2)
        cos_a = max(-1.0, min(1.0, cos_a))        # clamp for numerical safety
        branch_angles.append(math.degrees(math.acos(cos_a)))

    if branch_angles:
        ang_mean = sum(branch_angles) / len(branch_angles)
        ang_std = math.sqrt(
            sum((a - ang_mean) ** 2 for a in branch_angles) / len(branch_angles)
        )
    else:
        ang_mean = ang_std = float("nan")

    return {
        "mean_branch_angle_deg":          round(ang_mean, 2)
                                          if not math.isnan(ang_mean) else float("nan"),
        "std_branch_angle_deg":           round(ang_std, 2)
                                          if not math.isnan(ang_std) else float("nan"),
        "terminal_section_length_mean_px": term_mean,
        "section_length_cv":              section_cv,
    }


def compute_metrics(swc_path):
    """Compute morphology metrics from a single SWC file."""
    nodes = parse_swc(swc_path)
    if not nodes:
        return None

    children = {nid: [] for nid in nodes}
    root_id = None
    for nid, (x, y, z, pid) in nodes.items():
        if pid == -1:
            root_id = nid
        elif pid in children:
            children[pid].append(nid)

    if root_id is None:
        return None

    n_nodes = len(nodes)
    n_tips        = sum(1 for nid in nodes if len(children[nid]) == 0 and nid != root_id)
    n_branchpoints = sum(1 for nid in nodes if len(children[nid]) >= 2)
    n_sections    = n_nodes - 1  # tree has exactly n_nodes-1 edges

    # Total Euclidean length
    total_length = 0.0
    for nid, (x, y, z, pid) in nodes.items():
        if pid != -1 and pid in nodes:
            px, py, pz, _ = nodes[pid]
            total_length += math.sqrt((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2)

    mean_section_length = total_length / n_sections if n_sections > 0 else 0.0

    # BFS: accumulate root-to-node path distances
    root_x, root_y, root_z, _ = nodes[root_id]
    dist_from_root = {root_id: 0.0}
    bfs_order = []
    q = deque([root_id])
    max_path = 0.0
    while q:
        nid = q.popleft()
        bfs_order.append(nid)
        x, y, z, _ = nodes[nid]
        for child in children[nid]:
            cx, cy, cz, _ = nodes[child]
            d = math.sqrt((cx - x) ** 2 + (cy - y) ** 2 + (cz - z) ** 2)
            dist_from_root[child] = dist_from_root[nid] + d
            if dist_from_root[child] > max_path:
                max_path = dist_from_root[child]
            q.append(child)

    # --- new metrics ---

    # Spatial extents: width, height, aspect ratio
    all_x = [v[0] for v in nodes.values()]
    all_y = [v[1] for v in nodes.values()]
    width_px  = max(all_x) - min(all_x)
    height_px = max(all_y) - min(all_y)
    aspect_ratio = (width_px / height_px) if height_px > 0.0 else float("nan")

    # Convex hull area (2D, x-y plane)
    if _SCIPY and n_nodes >= 3:
        try:
            pts = np.array([[v[0], v[1]] for v in nodes.values()])
            hull = ConvexHull(pts)
            convex_hull_area = round(hull.volume, 2)   # hull.volume = area in 2D
        except Exception:
            convex_hull_area = float("nan")
    else:
        convex_hull_area = float("nan")

    # Per-tip: path length and tortuosity
    tip_ids = [nid for nid in nodes if len(children[nid]) == 0 and nid != root_id]
    tip_path_lengths = [dist_from_root[t] for t in tip_ids]

    tortuosities = []
    for t in tip_ids:
        tx, ty, tz, _ = nodes[t]
        straight = math.sqrt(
            (tx - root_x) ** 2 + (ty - root_y) ** 2 + (tz - root_z) ** 2
        )
        if straight > 0.0:
            tortuosities.append(dist_from_root[t] / straight)

    mean_path_length = (
        sum(tip_path_lengths) / len(tip_path_lengths) if tip_path_lengths else 0.0
    )
    mean_tortuosity = (
        sum(tortuosities) / len(tortuosities) if tortuosities else float("nan")
    )
    max_tortuosity = max(tortuosities) if tortuosities else float("nan")

    return {
        "n_nodes":                n_nodes,
        "n_tips":                 n_tips,
        "n_branchpoints":         n_branchpoints,
        "n_sections":             n_sections,
        "total_length_px":        round(total_length, 2),
        "mean_section_length_px": round(mean_section_length, 2),
        "max_path_length_px":     round(max_path, 2),
        "mean_path_length_px":    round(mean_path_length, 2),
        "width_px":               round(width_px, 2),
        "height_px":              round(height_px, 2),
        "aspect_ratio":           round(aspect_ratio, 4) if not math.isnan(aspect_ratio) else float("nan"),
        "convex_hull_area_px2":   convex_hull_area,
        "mean_tortuosity":        round(mean_tortuosity, 4) if not math.isnan(mean_tortuosity) else float("nan"),
        "max_tortuosity":         round(max_tortuosity, 4) if not math.isnan(max_tortuosity) else float("nan"),
        # --- refinement batch ---
        "occupancy_density":      round(total_length / convex_hull_area, 6)
                                  if (not math.isnan(convex_hull_area) and convex_hull_area > 0.0)
                                  else float("nan"),
        "median_path_length_px":  round(statistics.median(tip_path_lengths), 2)
                                  if tip_path_lengths else float("nan"),
        "path_length_cv":         round(
                                      math.sqrt(sum((p - mean_path_length) ** 2 for p in tip_path_lengths)
                                                / len(tip_path_lengths)) / mean_path_length, 4
                                  ) if (len(tip_path_lengths) > 0 and mean_path_length > 0.0)
                                  else float("nan"),
        "tip_branch_ratio":       round(n_tips / n_branchpoints, 4)
                                  if n_branchpoints > 0 else float("nan"),
        # --- asymmetry batch ---
        **_asymmetry_metrics(nodes, children, root_id, root_x, bfs_order),
        # --- branching geometry batch ---
        **_branching_geometry_metrics(nodes, children),
    }


def parse_identifiers(swc_path):
    """
    Extract case name and timestep from a path like:
      .../Cases/{case}/Simulation2D/merged/swc/{timestep}/clean_skeleton.swc
    """
    parts = os.path.normpath(swc_path).split(os.sep)
    try:
        idx = parts.index("Simulation2D")
        case_name = parts[idx - 1]
        timestep  = parts[idx + 3]   # Simulation2D / merged / swc / {timestep}
    except (ValueError, IndexError):
        case_name = "unknown"
        timestep  = "unknown"
    return case_name, timestep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--path", required=True, help="Base directory containing all cases")
    ap.add_argument("--out", default="morphometrics.csv", help="Output CSV path")
    args = ap.parse_args()

    pattern = os.path.join(args.path, "**", "clean_skeleton.swc")
    swc_files = sorted(glob.glob(pattern, recursive=True))

    if not swc_files:
        print(f"WARNING: no clean_skeleton.swc files found under {args.path}")
        return

    print(f"Found {len(swc_files)} SWC files")

    rows = []
    for swc_path in swc_files:
        case_name, timestep = parse_identifiers(swc_path)
        print(f"  {case_name} / t={timestep} ...", end="", flush=True)
        metrics = compute_metrics(swc_path)
        if metrics is None:
            print(" SKIP (empty or invalid)")
            continue
        row = {"case": case_name, "timestep": timestep, "swc_file": swc_path}
        row.update(metrics)
        rows.append(row)
        print(f" nodes={metrics['n_nodes']}  tips={metrics['n_tips']}  "
              f"BP={metrics['n_branchpoints']}  L={metrics['total_length_px']}")

    if not rows:
        print("No valid SWC files processed.")
        return

    fieldnames = list(rows[0].keys())
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
