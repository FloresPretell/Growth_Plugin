#!/usr/bin/env python3
"""
figdesign_B_gc_timecourse.py — growth-cone-only timecourses, grouped panels.

Single figure with four panels:
    1) Tubulin (u_t)              — its own y-axis
    2) MAP2 free / bound / phospho (u_u, u_b, u_p) — shared y-axis
    3) Calcium and Inhibitor (u_ca_cyt, inhibitor) — own panel; each only
       included if it has nontrivial dynamic range
    4) Lmax(t) and elongation_rate (dL/dt) on twin axes

Each line: GC-region mean, with ±SD shaded band, plotted against time_step.
Plateau onset (last step where dL/dt > epsilon · max(dL/dt)) marked.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figdesign_common import (
    FIELD_COLORS, FIELD_LABELS, MAP2_GROUP, add_exploratory_footer,
    detect_meaningful_fields, l0_value, lmax_per_step, load_enriched,
    save_figure,
)


def gc_stats(df, field):
    """Return (steps, mean, std, n) over rows where is_growth_cone_region==1."""
    gc = df[df['is_growth_cone_region'] == 1]
    g = gc.groupby('time_step')[field].agg(['mean', 'std', 'count']).sort_index()
    return g.index.values, g['mean'].values, g['std'].fillna(0).values, \
           g['count'].values


def detect_plateau(time_steps, dLdt, eps=0.1):
    if len(dLdt) == 0:
        return None
    dLdt = np.asarray(dLdt, dtype=float)
    a = np.where(np.isfinite(dLdt), np.abs(dLdt), 0.0)
    peak = a.max()
    if peak <= 0:
        return None
    last_active = np.where(a > eps * peak)[0]
    if len(last_active) == 0:
        return None
    return int(time_steps[last_active[-1]])


def plot_one(df, out_dir):
    fig, axes = plt.subplots(4, 1, figsize=(7.5, 9.0), sharex=True)

    # ---- panel 1: tubulin
    if 'u_t' in df.columns:
        s, m, sd, n = gc_stats(df, 'u_t')
        axes[0].plot(s, m, color=FIELD_COLORS['u_t'], linewidth=2,
                     label=FIELD_LABELS['u_t'])
        axes[0].fill_between(s, m - sd, m + sd, color=FIELD_COLORS['u_t'],
                             alpha=0.18, label='±1 SD')
        axes[0].set_ylabel(FIELD_LABELS['u_t'])
        axes[0].legend(loc='best', fontsize=8)

    # ---- panel 2: MAP2 group (free / bound / phospho)
    map2_present = [f for f in MAP2_GROUP if f in df.columns]
    for f in map2_present:
        s, m, sd, n = gc_stats(df, f)
        axes[1].plot(s, m, color=FIELD_COLORS[f], linewidth=2,
                     label=FIELD_LABELS[f])
        axes[1].fill_between(s, m - sd, m + sd, color=FIELD_COLORS[f],
                             alpha=0.15)
    axes[1].set_ylabel('MAP2 variants (GC mean)')
    axes[1].legend(loc='best', fontsize=8)

    # ---- panel 3: calcium + inhibitor (only if meaningful)
    keep_p3, drop_p3, reasons = detect_meaningful_fields(
        df, ['u_ca_cyt', 'inhibitor'])
    note3 = ''
    if keep_p3:
        for f in keep_p3:
            s, m, sd, n = gc_stats(df, f)
            axes[2].plot(s, m, color=FIELD_COLORS[f], linewidth=2,
                         label=FIELD_LABELS[f])
            axes[2].fill_between(s, m - sd, m + sd, color=FIELD_COLORS[f],
                                 alpha=0.18)
    if drop_p3:
        note3 = ' (omitted as flat: ' + ', '.join(
            f"{f} [{reasons.get(f, '')}]" for f in drop_p3) + ')'
    axes[2].set_ylabel('Calcium / Inhibitor (GC mean)' + note3)
    if keep_p3:
        axes[2].legend(loc='best', fontsize=8)
    else:
        axes[2].text(0.5, 0.5,
                     'all signals in this group flat across timesteps\n'
                     '(see panel title)', transform=axes[2].transAxes,
                     ha='center', va='center', color='#666666', fontsize=10)

    # ---- panel 4: Lmax and dL/dt on twin axes
    lmax = lmax_per_step(df)
    elong = df.groupby('time_step')['elongation_rate'].first().reindex(
        lmax.index)
    ax4 = axes[3]
    ax4.plot(lmax.index.values, lmax.values, color='#333333', linewidth=2,
             label='Lmax(t)')
    ax4.set_ylabel('Lmax(t) (path length)')
    ax4.axhline(l0_value(df), color='#333333', linewidth=0.8, linestyle='--',
                alpha=0.7, label='L0 (initial)')

    ax4b = ax4.twinx()
    ax4b.plot(elong.index.values, elong.values, color='#cc0066', linewidth=1.8,
              alpha=0.85, label='dL/dt (elongation rate)')
    ax4b.set_ylabel('dL/dt', color='#cc0066')
    ax4b.tick_params(axis='y', labelcolor='#cc0066')
    plateau_t = detect_plateau(elong.index.values, elong.values, eps=0.1)
    if plateau_t is not None:
        for a in axes:
            a.axvline(plateau_t, color='#cc0066', linewidth=0.8,
                      linestyle=':', alpha=0.55)
        ax4b.text(plateau_t, 0.05, f' plateau ≈ t={plateau_t} ',
                  transform=ax4b.get_xaxis_transform(),
                  ha='left', va='bottom', fontsize=8, color='#cc0066',
                  bbox=dict(facecolor='white', edgecolor='#cc0066',
                             alpha=0.7, pad=1.5))
    h4, l4 = ax4.get_legend_handles_labels()
    h4b, l4b = ax4b.get_legend_handles_labels()
    ax4.legend(h4 + h4b, l4 + l4b, loc='best', fontsize=8)
    ax4.set_xlabel('time step')

    fig.suptitle('Growth-cone region — multi-panel timecourse (candidate redesign)',
                 fontsize=11)
    add_exploratory_footer(fig, '<enriched CSV>',
                           extra='GC region = nodes within r_gc of any descendant tip')
    fig.tight_layout(rect=[0, 0.025, 1, 0.96])
    paths = save_figure(fig, out_dir, 'B_growthcone_multipanel')
    return [str(p) for p in paths]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    df = load_enriched(args.csv)
    written = plot_one(df, args.out_dir)
    for p in written:
        print(f'  wrote {p}')


if __name__ == '__main__':
    main()
