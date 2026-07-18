#!/usr/bin/env python3
"""
figdesign_E_branch_summary.py — per-branch heatmaps, sorted by final length.

Reorders branches by their final path_length_from_soma (descending). For each
of the main fields with nontrivial dynamic range, produces a heatmap with
branches on the y-axis and time on the x-axis.

Excludes inhibitor (the flat-signal case Nicole flagged) and any other field
that fails the dynamic-range check.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

from figdesign_common import (
    FIELD_LABELS, MAIN_FIELDS, add_exploratory_footer,
    detect_meaningful_fields, load_enriched, percentile_clip, save_figure,
)


def branch_order(df):
    """Sort branches by their max path_length at the LAST time_step."""
    last_t = df['time_step'].max()
    g = (df[df['time_step'] == last_t]
         .groupby('branch_id')['path_length_from_soma'].max())
    return g.sort_values(ascending=False).index.tolist()


def plot_field(df, field, out_dir):
    branches = branch_order(df)
    if not branches:
        return []
    steps = sorted(df['time_step'].unique())
    grid = np.full((len(branches), len(steps)), np.nan)
    for r, b in enumerate(branches):
        for c, t in enumerate(steps):
            sub = df[(df['branch_id'] == b) & (df['time_step'] == t)]
            v = pd.to_numeric(sub[field], errors='coerce').dropna()
            if len(v):
                grid[r, c] = float(v.mean())
    if np.all(np.isnan(grid)):
        return []
    vlo, vhi = percentile_clip(grid, 2.0, 98.0)
    if vhi <= vlo:
        vhi = vlo + 1e-12

    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    im = ax.imshow(grid, origin='upper', aspect='auto',
                   extent=[steps[0], steps[-1], len(branches), 0],
                   cmap='viridis', norm=Normalize(vmin=vlo, vmax=vhi),
                   interpolation='nearest')
    ax.set_yticks(np.arange(len(branches)) + 0.5)
    ax.set_yticklabels([f'branch {b}' for b in branches], fontsize=7)
    ax.set_xlabel('time step')
    ax.set_ylabel('branch (sorted by final length, longest at top)')
    ax.set_title(f'Per-branch mean — {FIELD_LABELS.get(field, field)} '
                 f'(candidate redesign)')
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label('branch mean')
    miss = ~np.isfinite(grid)
    if miss.any():
        ax.imshow(np.where(miss, 1, np.nan),
                  origin='upper', aspect='auto',
                  extent=[steps[0], steps[-1], len(branches), 0],
                  cmap='gray_r', alpha=0.25)
    add_exploratory_footer(fig, '<enriched CSV>',
                           extra=f'field={field} · branch order = final length descending')
    fig.tight_layout()
    paths = save_figure(fig, out_dir, f'E_branch_heatmap_{field}')
    return [str(p) for p in paths]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    df = load_enriched(args.csv)
    # Per Nicole's guidance: skip inhibitor unless meaningful.
    fields_to_try = ['u_t', 'u_u', 'u_b', 'u_p', 'u_ca_cyt']
    kept, dropped, _ = detect_meaningful_fields(df, fields_to_try)
    for f in kept:
        for p in plot_field(df, f, args.out_dir):
            print(f'  wrote {p}')


if __name__ == '__main__':
    main()
