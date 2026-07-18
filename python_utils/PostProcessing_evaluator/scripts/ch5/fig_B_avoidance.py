#!/usr/bin/env python3
"""
fig_B_avoidance.py — CH5 Part B (avoidance regime).

§5.6 κ-vs-log10(I) scatter (curvature approximated from path topology since
    `curvature` is not in the enriched CSV).
§5.7 GC-region inhibitor distribution over time (heatmap), as a proxy for the
    cos-θ heatmap described in the methods doc (cos θ requires ∇I vector data
    that is not propagated through the SWC sampling).

For the current dataset, inhibitor is essentially zero everywhere
(max ≲ 5e-7). The figures are produced for completeness and as scaffolding
for future avoidance-active datasets.
"""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from utils_ch5 import (  # noqa: E402
    DT_SECONDS, add_exploratory_footer, load_enriched, per_step_value,
)


def approx_curvature_from_swc(df: pd.DataFrame) -> pd.Series:
    """
    Approximate per-node curvature using consecutive SWC arc-length triples:
    for nodes with both a parent and >=1 child within the same frame, compute
    discrete curvature from (parent, node, child) using the Menger formula.
    Returns Series indexed like the original df rows, NaN where undefined.
    """
    kappa = pd.Series(np.nan, index=df.index)
    if not {'x_phys', 'y_phys', 'parent_id', 'swc_node_id', 'time_step'}.issubset(df.columns):
        return kappa
    for t, sub in df.groupby('time_step'):
        sub = sub.reset_index().rename(columns={'index': 'orig_idx'})
        # parent map
        parent_of = dict(zip(sub['swc_node_id'], sub['parent_id']))
        coord = dict(zip(sub['swc_node_id'],
                          zip(sub['x_phys'].values, sub['y_phys'].values)))
        # children list
        children = {}
        for nid, pid in parent_of.items():
            children.setdefault(pid, []).append(nid)
        for _, row in sub.iterrows():
            nid = int(row['swc_node_id']); pid = int(row['parent_id'])
            if pid <= 0:
                continue
            if pid not in coord:
                continue
            childs = children.get(nid, [])
            if not childs:
                continue
            cid = childs[0]  # take first child for kappa estimate
            if cid not in coord:
                continue
            p = np.array(coord[pid]); q = np.array((row['x_phys'], row['y_phys']))
            c = np.array(coord[cid])
            a = np.linalg.norm(q - p); b = np.linalg.norm(c - q); d = np.linalg.norm(c - p)
            if a < 1e-9 or b < 1e-9 or d < 1e-9:
                continue
            # Menger curvature
            # area = abs((q-p) x (c-p)) / 2
            cross = (q[0] - p[0]) * (c[1] - p[1]) - (q[1] - p[1]) * (c[0] - p[0])
            area = abs(cross) / 2.0
            if area < 1e-12:
                continue
            k = 4 * area / (a * b * d)
            kappa.iloc[int(row['orig_idx'])] = k
    return kappa


def figure_kappa_vs_inhibitor(df, out_dir, csv_path):
    """§5.6 — κ vs log10(I) scatter at SWC nodes (regime stratification limited
    since I≈0 everywhere; we still color by GC vs non-GC for visual reference).
    """
    kappa = approx_curvature_from_swc(df)
    df = df.assign(kappa=kappa)
    I = pd.to_numeric(df['inhibitor'], errors='coerce')
    mask = np.isfinite(I) & np.isfinite(kappa)
    sub = df[mask]
    if len(sub) == 0:
        return None

    log10I = np.log10(np.maximum(sub['inhibitor'].values, 1e-15))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4),
                              gridspec_kw=dict(width_ratios=[1.0, 1.0]))

    # (a) κ vs log10(I)
    gc_mask = sub['is_growth_cone_region'] == 1
    axes[0].scatter(log10I[~gc_mask], sub.loc[~gc_mask, 'kappa'],
                     s=8, alpha=0.35, color='#888888', label='shaft / non-GC nodes')
    axes[0].scatter(log10I[gc_mask], sub.loc[gc_mask, 'kappa'],
                     s=24, alpha=0.85, color='#cc0066', edgecolor='black',
                     linewidth=0.3, label='growth-cone region')
    axes[0].set_xlabel(r'$\log_{10}(I)$  (inhibitor)')
    axes[0].set_ylabel(r'$\hat\kappa$  (Menger curvature, µm$^{-1}$)')
    axes[0].set_title('§5.6 κ vs log10(I)  — scaffolding for active-avoidance datasets')
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(loc='best', fontsize=8)
    txt = (f'inhibitor range: {sub["inhibitor"].min():.3g} .. '
           f'{sub["inhibitor"].max():.3g}\n'
           'I is essentially zero in this dataset → \n'
           'I_avoid / I_prune thresholds not active')
    axes[0].text(0.02, 0.98, txt, transform=axes[0].transAxes,
                  ha='left', va='top', fontsize=8,
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                            edgecolor='#666', alpha=0.85))

    # (b) κ vs path_length, color by time
    sc = axes[1].scatter(sub['path_length_from_soma'], sub['kappa'],
                          c=sub['time_step'], cmap='viridis',
                          s=14, alpha=0.8, edgecolor='black', linewidth=0.2)
    axes[1].set_xlabel('path_length_from_soma  (µm)')
    axes[1].set_ylabel(r'$\hat\kappa$  (µm$^{-1}$)')
    axes[1].set_title('curvature distribution along arclength')
    cbar = fig.colorbar(sc, ax=axes[1], pad=0.02, fraction=0.04)
    cbar.set_label('time step')
    axes[1].grid(True, alpha=0.25)

    fig.suptitle('CH5 §5.6 — curvature × inhibitor (proxy κ from SWC; '
                  'inhibitor flat in D7.5_V3)', fontsize=11)
    add_exploratory_footer(fig, csv_path,
                            extra='κ approximated via Menger curvature on SWC triplets; not from level-set')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_B_kappa_vs_inhibitor.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)

    return dict(n_nodes_used=int(len(sub)),
                inhibitor_min=float(sub['inhibitor'].min()),
                inhibitor_max=float(sub['inhibitor'].max()),
                kappa_min=float(sub['kappa'].min()),
                kappa_max=float(sub['kappa'].max()))


def figure_inhibitor_distribution(df, out_dir, csv_path, n_bins=24):
    """
    §5.7-proxy — distribution of inhibitor at GC nodes over time as a 2-D
    density heatmap (frame × log10(I) bin).
    """
    gc = df[df['is_growth_cone_region'] == 1]
    if len(gc) == 0:
        return None
    I = pd.to_numeric(gc['inhibitor'], errors='coerce').values
    I = np.maximum(I, 1e-15)
    log10I = np.log10(I)
    edges = np.linspace(log10I.min(), log10I.max(), n_bins + 1)
    steps = sorted(gc['time_step'].unique())
    grid = np.zeros((len(steps), n_bins))
    for i, t in enumerate(steps):
        sub = gc[gc['time_step'] == t]
        vals = np.log10(np.maximum(
            pd.to_numeric(sub['inhibitor'], errors='coerce').values, 1e-15))
        h, _ = np.histogram(vals, bins=edges)
        grid[i] = h / max(1, h.sum())

    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    extent = [edges[0], edges[-1], steps[-1], steps[0]]
    im = ax.imshow(grid, aspect='auto', origin='upper', extent=extent,
                    cmap='viridis')
    ax.invert_yaxis()
    ax.set_xlabel(r'$\log_{10}(I)$  at GC nodes')
    ax.set_ylabel('time step')
    ax.set_title('§5.7 proxy: inhibitor distribution at GC over time')
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label('per-step density')
    add_exploratory_footer(fig, csv_path,
                            extra=('proxy for cos-θ heatmap; '
                                   'true cos-θ needs ∇I vectors from VTU'))
    fig.tight_layout()
    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_B_inhibitor_distribution_over_time.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)
    return dict(I_range=[float(I.min()), float(I.max())],
                n_steps=len(steps))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    df = load_enriched(args.csv)
    s_kappa = figure_kappa_vs_inhibitor(df, out_dir, args.csv)
    s_dist  = figure_inhibitor_distribution(df, out_dir, args.csv)
    with open(out_dir / 'fig_B_summary.json', 'w') as fh:
        json.dump(dict(csv=args.csv, kappa_vs_I=s_kappa,
                        inhibitor_dist=s_dist), fh, indent=2, default=float)
    print(f'fig_B done. summary -> {out_dir / "fig_B_summary.json"}')


if __name__ == '__main__':
    main()
