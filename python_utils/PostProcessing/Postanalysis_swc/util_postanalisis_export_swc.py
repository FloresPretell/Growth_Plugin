import glob
import math
import os
from collections import deque
from pathlib import Path

import numpy as np
from PIL import Image
from skimage.color import rgb2hsv
from skimage.measure import label, regionprops
from skimage.morphology import (
    binary_closing,
    disk,
    remove_small_holes,
    remove_small_objects,
    skeletonize,
)
from skimage.segmentation import clear_border


NEIGHBORS_8 = [
    (-1, -1), (-1, 0), (-1, 1),
    ( 0, -1),           ( 0, 1),
    ( 1, -1), ( 1, 0),  ( 1, 1),
]


def load_rgba_or_rgb(image_path):
    img = Image.open(image_path)
    if img.mode == "RGBA":
        rgba = np.array(img)
        rgb = rgba[..., :3].astype(np.float32)
        alpha = rgba[..., 3:4].astype(np.float32) / 255.0
        black = np.zeros_like(rgb)
        composited = rgb * alpha + black * (1.0 - alpha)
        return composited.astype(np.uint8)
    return np.array(img.convert("RGB"))


def save_bool_image(arr, path):
    Image.fromarray((arr.astype(np.uint8) * 255)).save(path)


def segment_arbor(
    img_rgb,
    value_threshold=0.15,
    saturation_threshold=0.10,
    min_object_size=200,
    closing_radius=2,
    clear_border_objects=False,
):
    hsv = rgb2hsv(img_rgb)
    s = hsv[..., 1]
    v = hsv[..., 2]

    mask = (v > value_threshold) & ((s > saturation_threshold) | (v > 0.30))
    mask = binary_closing(mask, disk(closing_radius))
    mask = remove_small_objects(mask, min_size=min_object_size)
    mask = remove_small_holes(mask, area_threshold=200)

    if clear_border_objects:
        mask = clear_border(mask)

    lab = label(mask)
    if lab.max() == 0:
        raise RuntimeError("No foreground object detected. Lower thresholds.")

    regions = regionprops(lab)
    largest = max(regions, key=lambda r: r.area)
    return lab == largest.label


def skeleton_degree(skel, y, x):
    h, w = skel.shape
    deg = 0
    for dy, dx in NEIGHBORS_8:
        yy, xx = y + dy, x + dx
        if 0 <= yy < h and 0 <= xx < w and skel[yy, xx]:
            deg += 1
    return deg


def compute_skeleton(mask):
    return skeletonize(mask)


def classify_pixels(skel):
    endpoints = []
    branchpoints = []
    regular = []
    ys, xs = np.nonzero(skel)
    for y, x in zip(ys, xs):
        d = skeleton_degree(skel, y, x)
        if d == 1:
            endpoints.append((y, x))
        elif d >= 3:
            branchpoints.append((y, x))
        else:
            regular.append((y, x))
    return endpoints, branchpoints, regular


def choose_root(mask, skel):
    ys, xs = np.nonzero(mask)
    y_bottom = int(ys.max())
    x_center = int(xs.mean())
    pts = np.column_stack(np.nonzero(skel))
    d2 = (pts[:, 0] - y_bottom) ** 2 + (pts[:, 1] - x_center) ** 2
    idx = int(np.argmin(d2))
    return int(pts[idx, 0]), int(pts[idx, 1])


def raw_pixel_tree_from_skeleton(skel, root):
    ys, xs = np.nonzero(skel)
    pixels = [(int(y), int(x)) for y, x in zip(ys, xs)]
    pixset = set(pixels)

    parent = {root: None}
    order = [root]
    q = deque([root])

    while q:
        y, x = q.popleft()
        for dy, dx in NEIGHBORS_8:
            yy, xx = y + dy, x + dx
            nxt = (yy, xx)
            if nxt in pixset and nxt not in parent:
                parent[nxt] = (y, x)
                order.append(nxt)
                q.append(nxt)

    return order, parent


# ---------------------------------------------------------------------------
# Graph-cleaning pipeline (Phase 1 quality improvement)
# ---------------------------------------------------------------------------

def build_graph_from_skeleton(skel):
    """Build undirected weighted graph {node: {neighbor: dist}} from skeleton pixels."""
    ys, xs = np.nonzero(skel)
    nodes = set(zip(ys.tolist(), xs.tolist()))
    graph = {n: {} for n in nodes}
    for y, x in nodes:
        for dy, dx in NEIGHBORS_8:
            yy, xx = y + dy, x + dx
            if (yy, xx) in nodes:
                graph[(y, x)][(yy, xx)] = math.hypot(dy, dx)
    return graph


def merge_branchpoints(graph):
    """
    Merge clusters of adjacent branchpoint pixels into single representative nodes.
    Returns (updated_graph, node_map) where node_map[old_node] = representative.
    """
    bp_set = {n for n, nbrs in graph.items() if len(nbrs) >= 3}

    visited = set()
    clusters = []
    for bp in bp_set:
        if bp in visited:
            continue
        cluster = []
        q = deque([bp])
        visited.add(bp)
        while q:
            node = q.popleft()
            cluster.append(node)
            for nbr in graph[node]:
                if nbr in bp_set and nbr not in visited:
                    visited.add(nbr)
                    q.append(nbr)
        clusters.append(cluster)

    node_map = {n: n for n in graph}

    for cluster in clusters:
        if len(cluster) == 1:
            continue
        cy = sum(n[0] for n in cluster) / len(cluster)
        cx = sum(n[1] for n in cluster) / len(cluster)
        rep = min(cluster, key=lambda n: (n[0] - cy) ** 2 + (n[1] - cx) ** 2)
        cluster_set = set(cluster)

        # Collect external connections, keeping shortest distance per neighbor
        ext = {}
        for node in cluster:
            for nbr, _ in graph[node].items():
                if nbr not in cluster_set:
                    d = math.hypot(rep[0] - nbr[0], rep[1] - nbr[1])
                    if nbr not in ext or d < ext[nbr]:
                        ext[nbr] = d

        for node in cluster:
            del graph[node]
            node_map[node] = rep
        for nbr in ext:
            if nbr in graph:
                for node in cluster:
                    graph[nbr].pop(node, None)

        graph[rep] = ext
        for nbr, d in ext.items():
            if nbr in graph:
                graph[nbr][rep] = d

    return graph, node_map


def compress_edges(graph, protected=None):
    """
    Collapse degree-2 chains into single edges.
    protected: set of nodes that must not be removed (e.g. root).
    """
    if protected is None:
        protected = set()

    changed = True
    while changed:
        changed = False
        for node in list(graph.keys()):
            if node not in graph or node in protected:
                continue
            nbrs = list(graph[node].keys())
            if len(nbrs) != 2:
                continue
            a, b = nbrs[0], nbrs[1]
            if a not in graph or b not in graph:
                continue
            if a == b:
                # Isolated self-loop remnant — remove node
                del graph[node]
                graph[a].pop(node, None)
                changed = True
                continue
            new_dist = graph[node][a] + graph[node][b]
            graph[a].pop(node, None)
            graph[b].pop(node, None)
            graph[a][b] = new_dist
            graph[b][a] = new_dist
            del graph[node]
            changed = True

    return graph


def prune_spurs(graph, min_spur_length=5.0, protected=None):
    """
    Iteratively remove terminal branches shorter than min_spur_length pixels.
    protected: set of nodes that must not be removed (e.g. root).
    """
    if protected is None:
        protected = set()

    changed = True
    while changed:
        changed = False
        for node in list(graph.keys()):
            if node not in graph or node in protected:
                continue
            if len(graph[node]) != 1:
                continue
            nbr = next(iter(graph[node]))
            if graph[node][nbr] < min_spur_length:
                graph[nbr].pop(node, None)
                del graph[node]
                changed = True

    return graph


def export_clean_swc(graph, root, output_path):
    """
    Export clean SWC containing only junction nodes (root, branchpoints, endpoints).
    Tree structure is recovered via BFS from root.
    """
    if root not in graph:
        raise RuntimeError(f"Root {root} not found in cleaned graph.")

    parent = {root: None}
    order = [root]
    q = deque([root])
    while q:
        node = q.popleft()
        for nbr in graph[node]:
            if nbr not in parent:
                parent[nbr] = node
                order.append(nbr)
                q.append(nbr)

    node_id = {node: i + 1 for i, node in enumerate(order)}
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# Clean SWC: branchpoint clusters merged, degree-2 chains compressed, spurs pruned\n")
        f.write("# id type x y z radius parent\n")
        for node in order:
            y, x = node
            pid = -1 if parent[node] is None else node_id[parent[node]]
            ntype = 1 if pid == -1 else 3
            f.write(f"{node_id[node]} {ntype} {x:.3f} {-y:.3f} 0.000 1.000 {pid}\n")


def path_length_to_root(node, parent):
    total = 0.0
    cur = node
    while parent[cur] is not None:
        p = parent[cur]
        total += math.hypot(cur[0] - p[0], cur[1] - p[1])
        cur = p
    return total


def export_pixel_swc(order, parent, output_path):
    node_id = {node: i + 1 for i, node in enumerate(order)}
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# Raw pixel-wise SWC extracted from binary skeleton\n")
        f.write("# id type x y z radius parent\n")
        for node in order:
            y, x = node
            pid = -1 if parent[node] is None else node_id[parent[node]]
            ntype = 1 if pid == -1 else 3
            f.write(f"{node_id[node]} {ntype} {x:.3f} {-y:.3f} 0.000 1.000 {pid}\n")


def compute_basic_metrics(mask, skel, endpoints, branchpoints, parent):
    ys, xs = np.nonzero(mask)
    bbox_h = int(ys.max() - ys.min() + 1)
    bbox_w = int(xs.max() - xs.min() + 1)
    foreground_area = int(mask.sum())

    total_length = 0.0
    for node, p in parent.items():
        if p is not None:
            total_length += math.hypot(node[0] - p[0], node[1] - p[1])

    max_path = 0.0
    for tip in endpoints:
        if tip in parent:
            max_path = max(max_path, path_length_to_root(tip, parent))

    return {
        "foreground_area_px": foreground_area,
        "bounding_box_h_px": bbox_h,
        "bounding_box_w_px": bbox_w,
        "skeleton_pixels": int(skel.sum()),
        "raw_endpoints": int(len(endpoints)),
        "raw_branch_pixels": int(len(branchpoints)),
        "total_skeleton_length_px": float(total_length),
        "max_root_to_tip_path_px": float(max_path),
    }


def run_pipeline_one_png(image_path, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    img = load_rgba_or_rgb(image_path)
    mask = segment_arbor(img)
    skel = compute_skeleton(mask)
    endpoints, branchpoints, _ = classify_pixels(skel)
    root = choose_root(mask, skel)
    order, parent = raw_pixel_tree_from_skeleton(skel, root)

    # Pixel-level SWC kept as debug reference
    export_pixel_swc(order, parent, out / "pixel_skeleton.swc")

    # Graph-cleaning pipeline → clean_skeleton.swc
    graph = build_graph_from_skeleton(skel)
    graph, node_map = merge_branchpoints(graph)
    clean_root = node_map.get(root, root)
    graph = compress_edges(graph, protected={clean_root})
    graph = prune_spurs(graph, min_spur_length=5.0, protected={clean_root})
    export_clean_swc(graph, clean_root, out / "clean_skeleton.swc")

    clean_tips = sum(1 for n in graph if len(graph[n]) == 1 and n != clean_root)
    clean_bp   = sum(1 for n in graph if len(graph[n]) >= 3)
    print(f"    skeleton: {int(skel.sum())} px  →  clean graph: {len(graph)} nodes "
          f"({clean_bp} BP, {clean_tips} tips)")

    metrics = compute_basic_metrics(mask, skel, endpoints, branchpoints, parent)
    with open(out / "metrics.txt", "w", encoding="utf-8") as f:
        for k, v in metrics.items():
            f.write(f"{k}: {v}\n")
        f.write(f"root_pixel: {root}\n")


def export_one(path, outdir):
    print(f"Exporting SWC from: {path}")
    # path is already the folder containing the PNG files
    png_dir = path
    os.makedirs(outdir, exist_ok=True)
    print(f"  png_dir = {png_dir}")
    print(f"  outdir = {outdir}")
    
    png_files = sorted(glob.glob(os.path.join(png_dir, "*.png")))
    if not png_files:
        print(f"  WARNING: no PNG files found in {png_dir}, skipping.")
        return

    for png_path in png_files:
        idx = os.path.splitext(os.path.basename(png_path))[0]
        out_dir = os.path.join(outdir, idx)
        print(f"  Processing {os.path.basename(png_path)} -> {out_dir}")
        run_pipeline_one_png(png_path, out_dir)

    print(f"Done: {path}")
