#!/usr/bin/env python3
"""
fig_D_global.py — CH5 Part D (global coupling).

§5.12 — spatial correlation map ρ(s) along arclength for each species pair.
§5.14-proxy — total-tubulin "mass" proxy (sum of u_t per step), and its
              finite-difference time-derivative compared with sources/sinks
              we can compute from the enriched CSV.

Neither of these uses VTU data; both are scaffolding for later VTU-side
Reynolds-budget work.
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
from matplotlib.colors import TwoSlopeNorm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from utils_ch5 import (   # noqa: E402
    DT_SECONDS, add_exploratory_footer, available_fields, load_enriched,
    make_correlation_norm, per_step_value, spatial_corr_per_bin,
)


def figure_spatial_corr(df, out_dir, csv_path, n_bins=None):
    """
    §5.12 — spatial correlation per arclength bin for species pairs.

    Case-adaptive design:
      - species pairs are derived from `available_fields()` (no hard-coded list).
      - `n_bins` defaults to sqrt(n_unique_steps) clipped to [6, 20].
      - colormap is built via `make_correlation_norm()` — always anchored at 0,
        symmetric vmin = -vmax, range fitted to data with a 0.05 pad and
        clipped to [-1, +1].
    """
    from itertools import combinations
    # Dynamic species list — Rule D + Rule E
    candidates = ['u_t', 'u_u', 'u_b', 'u_p', 'u_ca_cyt', 'inhibitor']
    fields, skipped = available_fields(df, candidates,
                                        min_finite=10, min_std=1e-12)
    pairs = list(combinations(fields, 2))
    if not pairs:
        print('[fig_D] no available species pairs; skipping spatial-corr.')
        return None
    if n_bins is None:
        steps_n = int(df['time_step'].nunique())
        n_bins = int(max(6, min(20, round(steps_n ** 0.5))))

    rho_matrix = np.full((len(pairs), n_bins), np.nan)
    n_matrix = np.zeros((len(pairs), n_bins), dtype=int)
    centers = None
    for i, (fx, fy) in enumerate(pairs):
        c, rho, nused = spatial_corr_per_bin(df, fx, fy, n_bins=n_bins)
        rho_matrix[i] = rho; n_matrix[i] = nused; centers = c

    fig, ax = plt.subplots(figsize=(8.5, 0.45 * len(pairs) + 2.5))
    extent = [0, n_bins, len(pairs), 0]
    norm = make_correlation_norm(rho_matrix, pad=0.05)
    im = ax.imshow(rho_matrix, aspect='auto', cmap='RdBu_r', norm=norm,
                    extent=extent, origin='upper')
    ax.set_yticks(np.arange(len(pairs)) + 0.5)
    ax.set_yticklabels([f'{a}·{b}' for a, b in pairs])
    ax.set_xticks(np.arange(n_bins) + 0.5)
    ax.set_xticklabels([f'{centers[i]:.2f}' for i in range(n_bins)],
                        fontsize=7, rotation=45)
    ax.set_xlabel('arclength bin center  (path_length_from_soma; µm)')
    ax.set_title('§5.12 spatial correlation per arclength bin  (Pearson ρ over time)')
    # Hatch bins with too few data points
    for i in range(rho_matrix.shape[0]):
        for j in range(rho_matrix.shape[1]):
            if n_matrix[i, j] < 4:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                              hatch='///',
                                              edgecolor='#888888'))
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label(r'Pearson $\rho$ (over time)')
    # Caption noting auto-scaled colormap range
    vmax_used = float(norm.vmax)
    add_exploratory_footer(fig, csv_path,
                            extra=(f'colormap symmetric ±{vmax_used:.2f} '
                                   f'anchored at ρ=0; hatch = <4 timepoints; '
                                   f'n_pairs={len(pairs)}, n_bins={n_bins}'))
    fig.tight_layout()
    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_D_spatial_correlation_map.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)
    return dict(pairs=pairs, centers=centers.tolist(),
                rho=rho_matrix.tolist(),
                n=n_matrix.tolist())


def figure_mass_proxy(df, out_dir, csv_path):
    """
    §5.14 proxy — sum of u_t over all sampled nodes per timestep, and its
    finite difference. Not a true Reynolds budget (needs VTU).
    """
    steps = sorted(df['time_step'].unique())
    M = np.array([float(df.loc[df['time_step'] == t, 'u_t'].sum())
                   for t in steps])
    Mn = np.array([float((df['time_step'] == t).sum()) for t in steps])
    dM_dt = np.gradient(M) / 1.0  # per frame

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(steps, M, color='#1f77b4', linewidth=2, marker='o',
                  markersize=4, label=r'$M(t)=\sum_{i}u_t$ across sampled nodes')
    axes[0].set_xlabel('time step')
    axes[0].set_ylabel(r'$M_{T}(t)$  (sum of $u_t$, sampled nodes)')
    axes[0].set_title('§5.14 proxy — tubulin "mass" sum')
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(loc='best', fontsize=8)

    axes[1].plot(steps, dM_dt, color='#cc0066', linewidth=1.8, marker='s',
                  markersize=4)
    axes[1].axhline(0, color='#888', linewidth=0.8)
    axes[1].set_xlabel('time step')
    axes[1].set_ylabel(r'$dM_{T}/dt$  (per frame)')
    axes[1].set_title('finite-difference derivative')
    axes[1].grid(True, alpha=0.25)
    txt = (f'CAVEAT: this is the sum across SWC-sampled nodes,\n'
           f'NOT the true ∫u_t dV. Node count varies with time\n'
           f'as the skeleton lengthens (n step 0 = {int(Mn[0])}, '
           f'n last = {int(Mn[-1])}).')
    axes[1].text(0.02, 0.02, txt, transform=axes[1].transAxes,
                  ha='left', va='bottom', fontsize=7.5,
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                            edgecolor='#666', alpha=0.85))

    fig.suptitle('CH5 §5.14 proxy — tubulin mass proxy '
                  '(true Reynolds budget needs VTU)', fontsize=11)
    add_exploratory_footer(fig, csv_path,
                            extra='node-count-weighted sum; not a volume integral')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_D_tubulin_mass_proxy.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)
    return dict(M_per_step=M.tolist(), n_per_step=Mn.tolist(),
                dM_dt=dM_dt.tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    df = load_enriched(args.csv)
    s12 = figure_spatial_corr(df, out_dir, args.csv)
    s14 = figure_mass_proxy(df, out_dir, args.csv)
    with open(out_dir / 'fig_D_summary.json', 'w') as fh:
        json.dump(dict(csv=args.csv, spatial_correlation=s12,
                        mass_proxy=s14), fh, indent=2, default=float)
    print(f'fig_D done. summary -> {out_dir / "fig_D_summary.json"}')


if __name__ == '__main__':
    main()
