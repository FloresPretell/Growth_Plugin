#!/usr/bin/env python3
"""
figdesign_C_lag_corr.py — concentration-vs-growth lag analysis (exploratory).

For each of {u_t, u_u, u_b, u_p, u_ca_cyt}, plot growth-cone mean concentration
at time t against elongation rate dL/dt at time t+lag, for lag in {0, 1, 2}.
Annotate Pearson r, Spearman rho, and N. Add an explicit autocorrelation
caveat box. Overlay a rolling-mean trend.

This is exploratory: each timepoint is not an independent sample, so the
p-values that would normally accompany r are misleading. We deliberately
do NOT print p-values.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from figdesign_common import (
    FIELD_COLORS, FIELD_LABELS, add_exploratory_footer,
    detect_meaningful_fields, load_enriched, save_figure,
)


def gc_mean(df, field):
    gc = df[df['is_growth_cone_region'] == 1]
    g = gc.groupby('time_step')[field].mean().sort_index()
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--lags', type=int, nargs='+', default=[0, 1, 2])
    args = ap.parse_args()

    df = load_enriched(args.csv)
    elong = df.groupby('time_step')['elongation_rate'].first().sort_index()
    fields = [f for f in ['u_t', 'u_u', 'u_b', 'u_p', 'u_ca_cyt']
              if f in df.columns]
    kept, dropped, reasons = detect_meaningful_fields(df, fields)
    if not kept:
        print('[C] no meaningful fields; skipping.')
        return

    n_rows = len(kept)
    n_cols = len(args.lags)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.5 * n_cols, 2.7 * n_rows),
                             sharex=False)
    if n_rows == 1:
        axes = np.array([axes])
    if n_cols == 1:
        axes = axes.reshape(-1, 1)

    for i, f in enumerate(kept):
        conc = gc_mean(df, f)
        steps = sorted(set(conc.index) & set(elong.index))
        c_all = conc.reindex(steps).values
        for j, lag in enumerate(args.lags):
            ax = axes[i][j]
            # Pair c(t) with elong(t + lag)
            t_arr = np.array(steps)
            valid = t_arr + lag <= t_arr.max()
            t_pair = t_arr[valid]
            c_pair = c_all[valid]
            e_pair = elong.reindex(t_pair + lag).values
            m = np.isfinite(c_pair) & np.isfinite(e_pair)
            cp = c_pair[m]; ep = e_pair[m]
            n_pts = len(cp)
            r_p = rho_s = np.nan
            if n_pts >= 4:
                r_p = stats.pearsonr(cp, ep).statistic
                rho_s = stats.spearmanr(cp, ep).statistic

            sc = ax.scatter(cp, ep, c=t_pair[m], cmap='viridis', s=24,
                            edgecolor='black', linewidth=0.4, alpha=0.85)
            # rolling-window trend in time order (NOT a fit)
            order = np.argsort(t_pair[m])
            tw = 4
            if len(cp) >= tw:
                cp_ord = cp[order]; ep_ord = ep[order]
                roll_c = pd.Series(cp_ord).rolling(tw, min_periods=1).mean().values
                roll_e = pd.Series(ep_ord).rolling(tw, min_periods=1).mean().values
                ax.plot(roll_c, roll_e, color='#cc0066', linewidth=1.4,
                        alpha=0.75)
            ax.set_title(f'{FIELD_LABELS.get(f, f)}  ·  lag={lag}',
                         fontsize=9)
            if i == n_rows - 1:
                ax.set_xlabel(f'GC mean {f}(t)')
            if j == 0:
                ax.set_ylabel(f'dL/dt(t+lag)')
            ax.grid(True, alpha=0.25)
            txt = (f'N={n_pts}\n'
                   f"Pearson r={r_p:.2f}\n"
                   f"Spearman ρ={rho_s:.2f}\n"
                   '(autocorrelated; no p-values reported)')
            ax.text(0.02, 0.98, txt, transform=ax.transAxes,
                    ha='left', va='top', fontsize=7,
                    bbox=dict(boxstyle='round,pad=0.25',
                              facecolor='white', edgecolor='#888888',
                              alpha=0.85))

    fig.suptitle('GC concentration(t) vs elongation rate(t+lag) — exploratory '
                 '(no p-values; points temporally autocorrelated)',
                 fontsize=10)
    add_exploratory_footer(fig, '<enriched CSV>',
                           extra='colormap = time step; rolling-mean line is time-ordered, not a fit')
    fig.tight_layout(rect=[0, 0.025, 1, 0.96])
    paths = save_figure(fig, args.out_dir, 'C_lag_corr_concentration_vs_dLdt')
    for p in paths:
        print(f'  wrote {p}')


if __name__ == '__main__':
    main()
