#!/usr/bin/env python3
"""
figdesign_D_new_region.py — new-region profiles, simplified.

For each meaningful field, produce one figure with two panels:
    LEFT:  η-profile (normalized_new_region_pos), 5 selected timepoints —
           early, early-mid, mid, late-mid, late — chosen from the available
           range of timesteps in the new region.
    RIGHT: heatmap (η × time_step) with the same field; cleaner replacement
           for the original figE which overlaid every timestep.

η = 0 is the initial front (just past L0). η = 1 is the current moving tip.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from figdesign_common import (
    FIELD_COLORS, FIELD_LABELS, MAIN_FIELDS, add_exploratory_footer,
    detect_meaningful_fields, load_enriched, percentile_clip, save_figure,
)


def select_timepoints(steps, k=5):
    if len(steps) <= k:
        return list(steps)
    qs = np.linspace(0.1, 0.95, k)
    idx = (qs * (len(steps) - 1)).astype(int)
    return [int(steps[i]) for i in idx]


def plot_field(df, field, out_dir, n_bins=24):
    new = df[df['is_new_region'] == 1].dropna(
        subset=['normalized_new_region_pos', field])
    if len(new) == 0:
        return []
    steps_with_new = sorted(new['time_step'].unique())
    selected = select_timepoints(steps_with_new, k=5)

    edges = np.linspace(0.0, 1.0, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8),
                             gridspec_kw=dict(width_ratios=[1.0, 1.15]))

    # LEFT: selected timepoint profiles
    cmap_t = plt.colormaps.get_cmap('plasma')
    for t in selected:
        sub = new[new['time_step'] == t]
        if len(sub) == 0:
            continue
        idx = np.digitize(sub['normalized_new_region_pos'].values, edges) - 1
        idx = np.clip(idx, 0, n_bins - 1)
        vals = pd.to_numeric(sub[field], errors='coerce').values
        prof = np.full(n_bins, np.nan)
        for j in range(n_bins):
            m = (idx == j) & np.isfinite(vals)
            if m.any():
                prof[j] = float(vals[m].mean())
        color = cmap_t(steps_with_new.index(t) / max(1, len(steps_with_new) - 1))
        axes[0].plot(centers, prof, marker='o', color=color, linewidth=1.6,
                     label=f't={t}')
    axes[0].set_xlabel('η = (s − L0)/(Lmax − L0)   (0 = initial front, 1 = tip)')
    axes[0].set_ylabel(f'{FIELD_LABELS.get(field, field)}  (bin mean)')
    axes[0].set_title('selected timepoints')
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(fontsize=7, loc='best')

    # RIGHT: full heatmap (η, time)
    grid = np.full((len(steps_with_new), n_bins), np.nan)
    for i, t in enumerate(steps_with_new):
        sub = new[new['time_step'] == t]
        idx = np.digitize(sub['normalized_new_region_pos'].values, edges) - 1
        idx = np.clip(idx, 0, n_bins - 1)
        vals = pd.to_numeric(sub[field], errors='coerce').values
        for j in range(n_bins):
            m = (idx == j) & np.isfinite(vals)
            if m.any():
                grid[i, j] = float(vals[m].mean())
    vlo, vhi = percentile_clip(grid, 2.0, 98.0)
    if vhi <= vlo:
        vhi = vlo + 1e-12
    im = axes[1].imshow(grid, origin='upper', aspect='auto',
                        extent=[0, 1, steps_with_new[-1], steps_with_new[0]],
                        cmap='viridis', norm=Normalize(vmin=vlo, vmax=vhi))
    axes[1].invert_yaxis()
    axes[1].set_xlabel('η')
    axes[1].set_ylabel('time step')
    axes[1].set_title('all timepoints (heatmap)')
    cbar = fig.colorbar(im, ax=axes[1], fraction=0.05, pad=0.02)
    cbar.set_label('bin mean')

    fig.suptitle(f'New-region profile (candidate redesign) — {FIELD_LABELS.get(field, field)}',
                 fontsize=11)
    add_exploratory_footer(fig, '<enriched CSV>',
                           extra=f'field={field}, n_bins={n_bins}')
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    paths = save_figure(fig, out_dir, f'D_new_region_profile_{field}')
    return [str(p) for p in paths]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--n-bins', type=int, default=24)
    args = ap.parse_args()
    df = load_enriched(args.csv)
    kept, dropped, _ = detect_meaningful_fields(df, MAIN_FIELDS)
    for f in kept:
        for p in plot_field(df, f, args.out_dir, n_bins=args.n_bins):
            print(f'  wrote {p}')


if __name__ == '__main__':
    main()
