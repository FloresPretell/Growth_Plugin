#!/usr/bin/env python3
"""
plot_front_aware_analysis.py
============================
Front-aware kymograph analysis of SWC-sampled simulation fields.

Requires the enriched CSV produced by enrich_swc_samples.py.
Runs in swc_env (numpy, pandas, matplotlib).

Usage:
    conda activate swc_env
    python plot_front_aware_analysis.py \\
        --csv  /scratch/flore0a/AnalysisResults/D7p5_V3_fields_enriched.csv \\
        --out-dir /scratch/flore0a/AnalysisResults/figures/D7p5_V3/front_aware_analysis

Figures produced (in <out-dir>/):
    figA_<field>_kymograph.png/svg     — Front-aware kymograph (one per field)
    figB_region_timecourse.png/svg     — Old vs new vs GC region mean ± SD
    figC_growthcone_dynamics.png/svg   — GC-region concentrations over time
    figC_growthcone_zscore.png/svg     — Same, z-score normalised
    figD_growthrate_vs_field.png/svg   — Elongation rate vs GC concentration (lag-1)
    figE_new_region_profile.png/svg    — Normalised new-region concentration profiles
    summary_report.txt                 — Node counts and region statistics
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats


# ─────────────────────────────────────────────────────────────────────────────
# Styling constants
# ─────────────────────────────────────────────────────────────────────────────

FIELD_LABELS = {
    'u_t':      'Tubulin',
    'u_u':      'MAP2 free',
    'u_b':      'MAP2 bound',
    'u_p':      'MAP2 phospho.',
    'u_ca_cyt': 'Calcium (cyt.)',
    'inhibitor': 'Inhibitor',
    'lsf':      'Level-set φ',
}
FIELD_COLORS = {
    'u_t':      '#e6194b',
    'u_u':      '#3cb44b',
    'u_b':      '#4363d8',
    'u_p':      '#f58231',
    'u_ca_cyt': '#911eb4',
    'inhibitor': '#42d4f4',
    'lsf':      '#999999',
}
REGION_STYLES = {
    'initial': dict(color='steelblue',  linestyle='-',  label='Initial region (s ≤ L0)'),
    'new':     dict(color='firebrick',  linestyle='-',  label='New region (s > L0)'),
    'gc':      dict(color='darkorange', linestyle='--', label='Growth-cone region'),
}

def flabel(f):  return FIELD_LABELS.get(f, f)
def fcolor(f):  return FIELD_COLORS.get(f, '#333333')


def save(fig, out_dir, stem, dpi=150):
    for ext in ('png', 'svg'):
        p = Path(out_dir) / f"{stem}.{ext}"
        fig.savefig(p, dpi=dpi, bbox_inches='tight')
    print(f"  → {Path(out_dir) / stem}.png")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

def load(csv_path):
    df = pd.read_csv(csv_path)
    numeric = [
        'time_step', 'time', 'path_length_from_soma', 'distance_from_soma',
        'branch_order', 'is_tip', 'valid_point',
        'L0', 'Lmax_t', 'elongation_rate',
        'is_initial_region', 'is_new_region', 'distance_beyond_initial_front',
        'normalized_global_position', 'normalized_new_region_pos',
        'is_growth_cone_region', 'branch_id',
    ]
    for c in numeric:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')

    # Auto-detect field columns (numeric, not metadata)
    meta = set(numeric + [
        'case_id', 'swc_frame', 'swc_node_id', 'parent_id', 'node_type',
        'x_px', 'y_px', 'x_phys', 'y_phys', 'z_phys',
    ])
    auto_fields = [c for c in df.columns
                   if c not in meta and pd.api.types.is_numeric_dtype(df[c])]

    if 'valid_point' in df.columns:
        n_inv = (df['valid_point'] == 0).sum()
        if n_inv:
            print(f"[info] Dropping {n_inv} invalid-probe rows")
        df = df[df['valid_point'] != 0].copy()

    return df, auto_fields


# ─────────────────────────────────────────────────────────────────────────────
# Figure A  —  front-aware kymograph (one per field)
# ─────────────────────────────────────────────────────────────────────────────

def fig_kymograph(df, fields, out_dir, case_id, n_bins=40):
    """
    Heatmap: x = path_length_from_soma, y = time_step, colour = concentration.
    Overlays:
      - vertical dashed white line at L0
      - cyan curve tracing Lmax(t)
    Bins with no SWC node are NaN (shown as white/grey).
    """
    L0 = float(df['L0'].iloc[0])
    path_max = df['path_length_from_soma'].max()
    bins = np.linspace(0, path_max * 1.02, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    df = df.copy()
    df['path_bin'] = pd.cut(df['path_length_from_soma'], bins=bins,
                             labels=bin_centers, include_lowest=True).astype(float)

    lmax_by_ts = df.groupby('time_step')['path_length_from_soma'].max()
    ts_values = np.array(sorted(df['time_step'].unique()))

    cmap = plt.cm.inferno.copy()
    cmap.set_bad(color='#cccccc')  # NaN → light grey

    for fname in fields:
        if fname not in df.columns:
            continue
        pivot = df.pivot_table(
            values=fname, index='time_step', columns='path_bin', aggfunc='mean'
        ).reindex(columns=bin_centers)  # keep all bins, missing → NaN

        fig, ax = plt.subplots(figsize=(10, 6))
        extent = [
            bin_centers[0], bin_centers[-1],
            ts_values[0] - 0.5, ts_values[-1] + 0.5,
        ]
        im = ax.imshow(
            pivot.values,
            aspect='auto', origin='lower', extent=extent, cmap=cmap,
        )
        plt.colorbar(im, ax=ax, label=flabel(fname), fraction=0.03, pad=0.02)

        # L0 — initial morphology boundary
        ax.axvline(L0, color='white', linestyle='--', linewidth=2.0,
                   label=f'L₀ = {L0:.2f} (initial front)')

        # Lmax(t) — moving front
        lmax_x = lmax_by_ts.reindex(ts_values).values
        ax.plot(lmax_x, ts_values, color='cyan', linewidth=1.8,
                marker='o', markersize=2.5, label='Lₘₐₓ(t)')

        ax.set_xlabel('Path length from soma (phys. units)', fontsize=11)
        ax.set_ylabel('Time step', fontsize=11)
        ax.set_title(f'{case_id}  —  {flabel(fname)} kymograph', fontsize=12)
        ax.legend(loc='upper left', fontsize=9, framealpha=0.7)

        save(fig, out_dir, f'figA_{fname}_kymograph')


# ─────────────────────────────────────────────────────────────────────────────
# Figure B  —  old vs new vs GC region time course
# ─────────────────────────────────────────────────────────────────────────────

def fig_region_timecourse(df, fields, out_dir, case_id):
    """
    For each field: mean ± 1 SD over time for three regions:
      initial (path ≤ L0), new (path > L0), growth-cone.
    """
    available = [f for f in fields if f in df.columns]
    if not available:
        return

    L0 = float(df['L0'].iloc[0])
    has_gc = 'is_growth_cone_region' in df.columns
    has_new = df['is_new_region'].sum() > 0

    n_fields = len(available)
    ncols = min(3, n_fields)
    nrows = (n_fields + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    fig.suptitle(f'{case_id} — Concentration by morphological region (mean ± 1 SD)',
                 fontsize=13)

    for fi, fname in enumerate(available):
        ax = axes[fi // ncols][fi % ncols]

        def plot_region(mask, style_key):
            sub = df[mask].groupby('time_step')[fname]
            mn = sub.mean()
            sd = sub.std().fillna(0)
            ts = mn.index.values
            kw = REGION_STYLES[style_key]
            ax.plot(ts, mn.values, **{k: v for k, v in kw.items() if k != 'label'},
                    linewidth=1.8, label=kw['label'])
            ax.fill_between(ts, (mn - sd).values, (mn + sd).values,
                            color=kw['color'], alpha=0.18)

        plot_region(df['is_initial_region'] == 1, 'initial')
        if has_new:
            plot_region(df['is_new_region'] == 1, 'new')
        if has_gc:
            plot_region(df['is_growth_cone_region'] == 1, 'gc')

        ax.set_xlabel('Time step')
        ax.set_ylabel(flabel(fname))
        ax.set_title(flabel(fname), fontsize=10)
        ax.grid(True, alpha=0.3)
        if fi == 0:
            ax.legend(fontsize=8)

    for fi in range(n_fields, nrows * ncols):
        axes[fi // ncols][fi % ncols].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, out_dir, 'figB_region_timecourse')


# ─────────────────────────────────────────────────────────────────────────────
# Figure C  —  growth-cone field dynamics
# ─────────────────────────────────────────────────────────────────────────────

def fig_growthcone_dynamics(df, fields, out_dir, case_id):
    """
    Mean concentration in the growth-cone region per time step, for all fields.
    Panel 1: raw concentrations.
    Panel 2: z-score normalised (zero-mean, unit-SD across time).
    """
    if 'is_growth_cone_region' not in df.columns:
        print("[warn] is_growth_cone_region not found; skipping fig C")
        return

    gc = df[df['is_growth_cone_region'] == 1]
    available = [f for f in fields if f in gc.columns]
    if not available:
        return

    gc_mean = gc.groupby('time_step')[available].mean()
    if gc_mean.empty:
        print("[warn] No growth-cone nodes found; skipping fig C")
        return

    ts = gc_mean.index.values

    # ── Raw ───────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    fig.suptitle(f'{case_id} — Growth-cone region: mean concentration over time', fontsize=12)

    ax = axes[0]
    for fname in available:
        ax.plot(ts, gc_mean[fname].values, color=fcolor(fname),
                linewidth=1.8, label=flabel(fname), marker='o', markersize=2)
    ax.set_xlabel('Time step')
    ax.set_ylabel('Concentration')
    ax.set_title('Raw concentrations')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)

    # ── Z-score ──────────────────────────────────────────────────────────────
    ax2 = axes[1]
    for fname in available:
        vals = gc_mean[fname].values.astype(float)
        mu, sigma = np.nanmean(vals), np.nanstd(vals)
        z = (vals - mu) / sigma if sigma > 1e-15 else np.zeros_like(vals)
        ax2.plot(ts, z, color=fcolor(fname), linewidth=1.8,
                 label=flabel(fname), marker='o', markersize=2)
    ax2.axhline(0, color='gray', linewidth=0.8, linestyle='--')
    ax2.set_xlabel('Time step')
    ax2.set_ylabel('Z-score')
    ax2.set_title('Z-score normalised (zero-mean, unit-SD across time)')
    ax2.legend(fontsize=8, ncol=2)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, out_dir, 'figC_growthcone_dynamics')


# ─────────────────────────────────────────────────────────────────────────────
# Figure D  —  elongation rate vs growth-cone concentration (lag-1)
# ─────────────────────────────────────────────────────────────────────────────

def fig_growthrate_vs_field(df, fields, out_dir, case_id):
    """
    Scatter: growth-cone mean concentration at time t  vs  elongation_rate at t+1.
    Tests whether local biochemical state predicts subsequent elongation.
    Pearson r and Spearman ρ reported.
    """
    if 'elongation_rate' not in df.columns:
        print("[warn] elongation_rate not found; skipping fig D")
        return
    if 'is_growth_cone_region' not in df.columns:
        print("[warn] is_growth_cone_region not found; skipping fig D")
        return

    gc = df[df['is_growth_cone_region'] == 1]
    available = [f for f in fields if f in gc.columns]
    if not available:
        return

    gc_mean = gc.groupby('time_step')[available].mean()
    elong_by_ts = df.groupby('time_step')['elongation_rate'].mean()

    # Lag-1: field at t vs rate at t+1
    ts_all = sorted(gc_mean.index)
    lag_data = []
    for ts in ts_all[:-1]:
        ts_next = ts + 1
        if ts_next not in elong_by_ts.index:
            continue
        row = {'ts': ts, 'rate_next': float(elong_by_ts[ts_next])}
        for fname in available:
            row[fname] = float(gc_mean.loc[ts, fname]) if ts in gc_mean.index else np.nan
        lag_data.append(row)

    if not lag_data:
        print("[warn] Not enough time steps for lag-1 analysis; skipping fig D")
        return

    lag_df = pd.DataFrame(lag_data).dropna(subset=['rate_next'])

    n_fields = len(available)
    ncols = min(3, n_fields)
    nrows = (n_fields + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)
    fig.suptitle(f'{case_id} — Elongation rate (t+1) vs GC concentration (t)\n'
                 f'[Tests whether GC state predicts subsequent growth]', fontsize=11)

    for fi, fname in enumerate(available):
        ax = axes[fi // ncols][fi % ncols]
        if fname not in lag_df.columns:
            ax.set_visible(False)
            continue

        sub = lag_df[['rate_next', fname]].dropna()
        if len(sub) < 4:
            ax.set_visible(False)
            continue

        x, y = sub[fname].values, sub['rate_next'].values
        ax.scatter(x, y, color=fcolor(fname), alpha=0.8, s=40, edgecolors='none')

        # Regression line
        slope, intercept, r, p, _ = stats.linregress(x, y)
        xs = np.linspace(x.min(), x.max(), 100)
        ax.plot(xs, slope * xs + intercept, color='black', linewidth=1.2, linestyle='--')

        rho, p_sp = stats.spearmanr(x, y)
        ax.set_title(f'{flabel(fname)}\nPearson r={r:.2f}  Spearman ρ={rho:.2f}', fontsize=9)
        ax.set_xlabel(f'GC {flabel(fname)} at t')
        ax.set_ylabel('Elongation rate at t+1')
        ax.grid(True, alpha=0.3)
        ax.text(0.97, 0.05, f'p={p:.3f}', transform=ax.transAxes,
                ha='right', va='bottom', fontsize=8, color='grey')

    for fi in range(n_fields, nrows * ncols):
        axes[fi // ncols][fi % ncols].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save(fig, out_dir, 'figD_growthrate_vs_field')


# ─────────────────────────────────────────────────────────────────────────────
# Figure E  —  normalised new-region concentration profile
# ─────────────────────────────────────────────────────────────────────────────

def fig_new_region_profile(df, fields, out_dir, case_id, n_eta_bins=20):
    """
    For the newly generated region only (path > L0), plot concentration vs
        η = (path_length - L0) / (Lmax(t) - L0)
    where η=0 is the initial-front boundary and η=1 is the current distal tip.
    Colour indicates time step.  Allows comparison across frames despite growing length.
    """
    new = df[df['is_new_region'] == 1].copy()
    if new.empty:
        print("[warn] No new-region nodes found; skipping fig E")
        return

    available = [f for f in fields if f in new.columns]
    if not available:
        return

    # Bin η
    eta_bins = np.linspace(0, 1, n_eta_bins + 1)
    eta_centers = 0.5 * (eta_bins[:-1] + eta_bins[1:])
    new['eta_bin'] = pd.cut(
        new['normalized_new_region_pos'], bins=eta_bins,
        labels=eta_centers, include_lowest=True
    ).astype(float)

    ts_all = sorted(new['time_step'].unique())
    try:
        cmap = matplotlib.colormaps.get_cmap('plasma')
    except AttributeError:
        cmap = cm.get_cmap('plasma')
    ts_norm = {ts: cmap(i / max(len(ts_all) - 1, 1)) for i, ts in enumerate(ts_all)}

    n_fields = len(available)
    ncols = min(3, n_fields)
    nrows = (n_fields + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    fig.suptitle(f'{case_id} — New-region concentration profile\n'
                 f'η = 0: initial-front boundary  │  η = 1: current distal tip',
                 fontsize=12)

    for fi, fname in enumerate(available):
        ax = axes[fi // ncols][fi % ncols]
        if fname not in new.columns:
            ax.set_visible(False)
            continue

        for ts in ts_all:
            sub_ts = new[new['time_step'] == ts].groupby('eta_bin')[fname].mean()
            if sub_ts.empty:
                continue
            ax.plot(sub_ts.index, sub_ts.values,
                    color=ts_norm[ts], linewidth=1.5, alpha=0.85, marker='o', markersize=3)

        ax.axvline(0.0, color='gray', linestyle=':', linewidth=1.0)
        ax.axvline(1.0, color='gray', linestyle=':', linewidth=1.0)
        ax.set_xlabel('η  (normalised position in new region)')
        ax.set_ylabel(flabel(fname))
        ax.set_title(flabel(fname), fontsize=10)
        ax.set_xlim(-0.05, 1.05)
        ax.grid(True, alpha=0.3)

    # Shared colorbar for time
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=plt.Normalize(vmin=ts_all[0], vmax=ts_all[-1]))
    sm.set_array([])
    fig.colorbar(sm, ax=axes.ravel().tolist(), label='Time step',
                 shrink=0.5, pad=0.02, fraction=0.02)

    for fi in range(n_fields, nrows * ncols):
        axes[fi // ncols][fi % ncols].set_visible(False)

    fig.tight_layout(rect=[0, 0, 0.93, 0.93])
    save(fig, out_dir, 'figE_new_region_profile')


# ─────────────────────────────────────────────────────────────────────────────
# Figure F  —  branch-level summary (if branch_id available)
# ─────────────────────────────────────────────────────────────────────────────

def fig_branch_summary(df, fields, out_dir, case_id):
    """
    Branch-by-time heatmap for each field: rows = branch_id, cols = time_step,
    colour = tip concentration on that branch at that time.

    Branch IDs are per-frame (inferred from SWC tree).  Cross-frame branch
    tracking is not implemented; branches are matched by their rank in
    path_length_from_soma at each time step (longest tip = branch 0, etc.).

    Skipped if branch_id == -1 for all rows (branch_id not computed).
    """
    if 'branch_id' not in df.columns or (df['branch_id'] == -1).all():
        print("[info] branch_id not available; skipping fig F "
              "(re-run enrich_swc_samples.py with --swc-dir to enable)")
        return

    tips = df[(df['is_tip'] == 1) & (df['branch_id'] >= 0)].copy()
    if tips.empty:
        return

    available = [f for f in fields if f in tips.columns]

    # Since branch IDs change per frame, we rank-sort tips by path_length
    # within each time step to create a stable "branch rank" for the heatmap.
    def rank_tips(group):
        group = group.sort_values('path_length_from_soma', ascending=False)
        group['branch_rank'] = range(len(group))
        return group

    tips = tips.groupby('time_step', group_keys=False).apply(rank_tips)
    max_rank = int(tips['branch_rank'].max())
    ts_all = sorted(tips['time_step'].unique())

    for fname in available:
        pivot = tips.pivot_table(
            values=fname, index='branch_rank', columns='time_step', aggfunc='mean'
        ).reindex(index=range(max_rank + 1), columns=ts_all)

        fig, ax = plt.subplots(figsize=(max(8, len(ts_all) * 0.4), max(4, max_rank + 2)))
        cmap = plt.cm.inferno.copy()
        cmap.set_bad('#cccccc')
        im = ax.imshow(pivot.values, aspect='auto', origin='upper', cmap=cmap)
        plt.colorbar(im, ax=ax, label=flabel(fname), fraction=0.03, pad=0.02)
        ax.set_xticks(range(len(ts_all)))
        ax.set_xticklabels(ts_all, fontsize=7, rotation=90)
        ax.set_yticks(range(max_rank + 1))
        ax.set_yticklabels([f'B{r}' for r in range(max_rank + 1)], fontsize=8)
        ax.set_xlabel('Time step')
        ax.set_ylabel('Branch rank (0 = longest)')
        ax.set_title(f'{case_id} — Tip {flabel(fname)} per branch over time', fontsize=11)
        fig.tight_layout()
        save(fig, out_dir, f'figF_branch_{fname}')


# ─────────────────────────────────────────────────────────────────────────────
# Summary report
# ─────────────────────────────────────────────────────────────────────────────

def write_summary(df, fields, out_dir, case_id, r_gc):
    lines = []
    lines.append(f"Front-aware analysis summary — {case_id}")
    lines.append("=" * 60)

    L0 = float(df['L0'].iloc[0])
    lines.append(f"\nL0 (initial max path length): {L0:.4f}")
    lines.append("  Computed as: max(path_length_from_soma) at the earliest time step")
    lines.append(f"\nGrowth-cone radius (r_gc): {r_gc:.4f} physical units")

    gc_method = ("SWC subtree" if (df['branch_id'] != -1).any()
                 else "CSV approximation (min downstream tip distance)")
    lines.append(f"Growth-cone method: {gc_method}")

    lines.append(f"\n{'ts':>4}  {'n':>4}  {'init':>6}  {'new':>5}  "
                 f"{'tip':>4}  {'gc':>5}  {'Lmax':>7}  {'dL/dt':>8}")
    for ts in sorted(df['time_step'].unique()):
        sub = df[df['time_step'] == ts]
        elong = sub['elongation_rate'].mean() if 'elongation_rate' in sub else np.nan
        lines.append(
            f"  {int(ts):4d}  {len(sub):4d}  "
            f"{int(sub['is_initial_region'].sum()):6d}  "
            f"{int(sub['is_new_region'].sum()):5d}  "
            f"{int(sub['is_tip'].sum()) if 'is_tip' in sub else -1:4d}  "
            f"{int(sub['is_growth_cone_region'].sum()):5d}  "
            f"{sub['Lmax_t'].iloc[0]:7.4f}  "
            f"{elong:8.5f}"
        )

    lines.append("\nField stats in growth-cone region (mean over all time steps):")
    gc = df[df['is_growth_cone_region'] == 1]
    for f in fields:
        if f in gc.columns:
            lines.append(f"  {f:12s}  mean={gc[f].mean():.4g}  std={gc[f].std():.4g}")

    report_path = Path(out_dir) / 'summary_report.txt'
    with open(report_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"  → {report_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--csv', required=True,
                        help='Enriched CSV from enrich_swc_samples.py')
    parser.add_argument('--out-dir', default='./front_aware_analysis',
                        help='Output directory for figures (default: ./front_aware_analysis)')
    parser.add_argument('--case-id', default=None,
                        help='Case label for titles (auto-inferred from CSV if omitted)')
    parser.add_argument('--fields',
                        default='u_t,u_u,u_b,u_p,u_ca_cyt,inhibitor',
                        help='Comma-separated fields to plot '
                             '(default: u_t,u_u,u_b,u_p,u_ca_cyt,inhibitor)')
    parser.add_argument('--kymograph-fields', default='u_t,u_ca_cyt,u_b,u_p',
                        help='Fields for kymograph figures A (default: u_t,u_ca_cyt,u_b,u_p)')
    parser.add_argument('--r-gc', type=float, default=0.1,
                        help='Growth-cone radius used in enrichment (for report only, '
                             'actual GC flags come from the CSV)')
    parser.add_argument('--n-kbins', type=int, default=40,
                        help='Path-length bins for kymograph (default: 40)')
    parser.add_argument('--n-eta-bins', type=int, default=20,
                        help='η bins for normalised new-region profile (default: 20)')
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[info] Loading {args.csv} ...")
    df, auto_fields = load(args.csv)
    print(f"[info] {len(df)} rows, auto-detected fields: {auto_fields}")

    case_id = args.case_id or (
        str(df['case_id'].iloc[0]) if 'case_id' in df.columns else 'case')

    # Required enriched columns
    required = {'L0', 'Lmax_t', 'is_initial_region', 'is_new_region',
                 'is_growth_cone_region', 'normalized_new_region_pos'}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(
            f"[ERROR] Missing enriched columns: {missing}\n"
            "        Run enrich_swc_samples.py first."
        )

    fields = [f.strip() for f in args.fields.split(',') if f.strip() and f in df.columns]
    kymo_fields = [f.strip() for f in args.kymograph_fields.split(',')
                   if f.strip() and f in df.columns]

    print(f"\n[fig A] Front-aware kymographs ({len(kymo_fields)} fields) ...")
    fig_kymograph(df, kymo_fields, out_dir, case_id, n_bins=args.n_kbins)

    print("[fig B] Old vs new vs GC region time course ...")
    fig_region_timecourse(df, fields, out_dir, case_id)

    print("[fig C] Growth-cone dynamics ...")
    fig_growthcone_dynamics(df, fields, out_dir, case_id)

    print("[fig D] Elongation rate vs GC concentration (lag-1) ...")
    fig_growthrate_vs_field(df, fields, out_dir, case_id)

    print("[fig E] Normalised new-region profile ...")
    fig_new_region_profile(df, fields, out_dir, case_id, n_eta_bins=args.n_eta_bins)

    print("[fig F] Branch-level summary ...")
    fig_branch_summary(df, fields, out_dir, case_id)

    print("[report] Writing summary ...")
    write_summary(df, fields, out_dir, case_id, r_gc=args.r_gc)

    print(f"\n[done] All figures in {out_dir}/")


if __name__ == '__main__':
    main()
