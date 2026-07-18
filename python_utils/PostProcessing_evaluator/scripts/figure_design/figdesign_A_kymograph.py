#!/usr/bin/env python3
"""
figdesign_A_kymograph.py — clearer per-field kymographs (candidate redesign).

For each meaningful field, produce two stacked panels:
    top:    raw concentration heatmap with colorbar clipped to 2-98th percentile
    bottom: same data normalized to the field's value at t=0 in each path bin
            (relative change vs t=0; diverging colormap centered on 1.0)

Differences from the original figA kymograph:
    - colorbar clipped to per-field 2-98 percentiles (prevents one outlier
      from washing out the rest);
    - grey "no skeleton here" cells are hatched + labeled in the legend
      instead of dominating the visual field;
    - explicit L0 line + Lmax(t) line, both labeled in the legend;
    - x-axis caption explains that path_length is the cumulative SWC
      arc-length from soma, summed over all branches (NOT a single ray);
    - exploratory footer; cannot be mistaken for a publication figure.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm, Normalize
from matplotlib.patches import Patch

from figdesign_common import (
    FIELD_LABELS, MAIN_FIELDS, add_exploratory_footer, detect_meaningful_fields,
    l0_value, lmax_per_step, load_enriched, percentile_clip, save_figure,
    time_values,
)


def bin_field_to_grid(df, field, n_bins=80):
    """Return (time_index, edges, mean_grid, count_grid) for the field."""
    steps = sorted(df['time_step'].unique())
    pmax = float(df['path_length_from_soma'].max())
    edges = np.linspace(0.0, pmax, n_bins + 1)
    mean_grid = np.full((len(steps), n_bins), np.nan, dtype=float)
    count_grid = np.zeros((len(steps), n_bins), dtype=int)
    for i, s in enumerate(steps):
        sub = df[df['time_step'] == s]
        if len(sub) == 0:
            continue
        idx = np.digitize(sub['path_length_from_soma'].values, edges) - 1
        idx = np.clip(idx, 0, n_bins - 1)
        vals = pd.to_numeric(sub[field], errors='coerce').values
        for j in range(n_bins):
            mask = idx == j
            count_grid[i, j] = int(mask.sum())
            if mask.any():
                v = vals[mask]
                v = v[np.isfinite(v)]
                if v.size:
                    mean_grid[i, j] = float(v.mean())
    return steps, edges, mean_grid, count_grid


def plot_kymograph_pair(df, field, out_dir, n_bins=80):
    steps, edges, grid, counts = bin_field_to_grid(df, field, n_bins)
    if np.all(np.isnan(grid)):
        return []
    centers = 0.5 * (edges[:-1] + edges[1:])
    extent = [edges[0], edges[-1], steps[-1], steps[0]]

    # Raw
    vlo, vhi = percentile_clip(grid[np.isfinite(grid)], 2.0, 98.0)
    if vhi <= vlo:
        vhi = vlo + 1e-12
    norm = Normalize(vmin=vlo, vmax=vhi, clip=True)

    # Relative-to-t0
    base = grid[0]                                    # bin profile at first available step
    rel = np.full_like(grid, np.nan, dtype=float)
    for i in range(grid.shape[0]):
        with np.errstate(divide='ignore', invalid='ignore'):
            rel[i] = np.where(np.isfinite(base) & (np.abs(base) > 1e-12),
                              grid[i] / base, np.nan)
    rng = np.nanpercentile(np.abs(rel[np.isfinite(rel)] - 1.0), 95) \
        if np.isfinite(rel).any() else 1.0
    half = max(rng, 1e-3)
    div_norm = TwoSlopeNorm(vcenter=1.0,
                            vmin=max(1.0 - half, 0.0),
                            vmax=1.0 + half)

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 6.6), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.0, 1.0]))

    im0 = axes[0].imshow(grid, aspect='auto', origin='upper', extent=extent,
                         cmap='viridis', norm=norm,
                         interpolation='nearest')
    axes[0].set_title(f'{FIELD_LABELS.get(field, field)} — raw '
                      f'(colorbar clipped to 2-98th percentile)')
    axes[0].set_ylabel('time step (later → top)')

    im1 = axes[1].imshow(rel, aspect='auto', origin='upper', extent=extent,
                         cmap='RdBu_r', norm=div_norm,
                         interpolation='nearest')
    axes[1].set_title(f'{FIELD_LABELS.get(field, field)} — '
                      f'relative to t=0 in same bin (1.0 = no change)')
    axes[1].set_xlabel('path_length_from_soma — cumulative SWC arc-length '
                       '(branches concatenated, NOT a single ray)')
    axes[1].set_ylabel('time step')

    # Overlay Lmax(t) and L0 on both panels.
    L0 = l0_value(df)
    lmax = lmax_per_step(df).reindex(steps).values
    for ax in axes:
        ax.axvline(L0, color='#222222', linewidth=1.2, linestyle='--',
                   label='L0 (initial front)')
        ax.plot(lmax, np.array(steps), color='#ffffff', linewidth=2.4,
                solid_capstyle='round')
        ax.plot(lmax, np.array(steps), color='#000000', linewidth=1.4,
                solid_capstyle='round', label='Lmax(t) — moving tip')
        ax.invert_yaxis()  # time grows downward like the original

    # Mask hatch legend handle.
    masked_handle = Patch(facecolor='#dddddd', edgecolor='#888888',
                          hatch='///', label='no SWC nodes in bin')

    handles_top, labels_top = axes[0].get_legend_handles_labels()
    handles_top.append(masked_handle); labels_top.append(masked_handle.get_label())
    axes[0].legend(handles_top, labels_top, loc='upper right',
                   fontsize=7, framealpha=0.85)

    cbar0 = fig.colorbar(im0, ax=axes[0], fraction=0.04, pad=0.02)
    cbar0.set_label('concentration')
    cbar1 = fig.colorbar(im1, ax=axes[1], fraction=0.04, pad=0.02)
    cbar1.set_label('c(t)/c(t=0)  (1.0 = no change)')

    # Render the missing-bin hatch beneath the heatmap.
    for ax, data in [(axes[0], grid), (axes[1], rel)]:
        miss = ~np.isfinite(data)
        if miss.any():
            ax.imshow(np.where(miss, 1, np.nan),
                      aspect='auto', origin='upper', extent=extent,
                      cmap='gray_r', alpha=0.18, interpolation='nearest')

    add_exploratory_footer(fig, '<enriched CSV>',
                           extra=f'field={field}, n_bins={n_bins}')
    fig.suptitle(f'kymograph (candidate redesign) — {FIELD_LABELS.get(field, field)}',
                 fontsize=11)
    fig.tight_layout(rect=[0, 0.025, 1, 0.96])
    paths = save_figure(fig, out_dir,
                        f'A_kymograph_{field}_raw_and_relative')
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--n-bins', type=int, default=80)
    args = ap.parse_args()

    df = load_enriched(args.csv)
    kept, dropped, reasons = detect_meaningful_fields(df, MAIN_FIELDS)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f'[A_kymograph] kept fields: {kept}')
    print(f'[A_kymograph] dropped (flat/missing): '
          f'{[(f, reasons[f]) for f in dropped]}')

    written = []
    for field in kept:
        paths = plot_kymograph_pair(df, field, out_dir, n_bins=args.n_bins)
        for p in paths:
            print(f'  wrote {p}')
            written.append(str(p))
    return written


if __name__ == '__main__':
    main()
