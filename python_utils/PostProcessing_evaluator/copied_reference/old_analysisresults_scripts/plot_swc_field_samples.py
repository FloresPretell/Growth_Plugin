#!/usr/bin/env python3
"""
plot_swc_field_samples.py
=========================
Generate figures from the long-format CSV produced by sample_fields_on_swc.py.

Requires: numpy, pandas, matplotlib  (available in swc_env)
    conda activate swc_env
    python plot_swc_field_samples.py --csv output.csv --out-dir ./figures

────────────────────────────────────────────────────────────────────────────
FIGURES PRODUCED
────────────────────────────────────────────────────────────────────────────
1. concentration_vs_distance.png
   — Concentration of each field vs. distance_from_soma,
     shown for a small set of selected time steps.
     Branchpoints and tips are highlighted.

2. tip_concentrations_over_time.png
   — Concentration of each field at tip nodes (growth cones)
     plotted as a function of time (one line per tip's path_length_from_soma).

3. heatmap_path_vs_time.png
   — 2D heatmap: x-axis = path_length_from_soma (binned),
     y-axis = time_step, colour = median concentration.
     Separate panels for tubulin (u_t) and calcium (u_ca_cyt).

4. field_distributions_over_time.png
   — Violin / box plot showing the distribution of each field value
     across all nodes at each time step.

────────────────────────────────────────────────────────────────────────────
USAGE
────────────────────────────────────────────────────────────────────────────

    python plot_swc_field_samples.py \\
        --csv /scratch/flore0a/AnalysisResults/D7p5_V3_fields.csv \\
        --out-dir ./figures/D7p5_V3 \\
        --case-id D7p5_V3

    # Choose specific fields to plot (comma-separated):
    python plot_swc_field_samples.py --csv ... --fields u_t,u_ca_cyt

    # Choose time steps for figure 1 (0-indexed, comma-separated):
    python plot_swc_field_samples.py --csv ... --time-steps 0,10,20,40
"""

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
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
    'lsf':      '#808080',
}


def label(fname):
    return FIELD_LABELS.get(fname, fname)


def color(fname):
    return FIELD_COLORS.get(fname, '#333333')


def load_csv(csv_path):
    df = pd.read_csv(csv_path)
    numeric_cols = [
        'time_step', 'time', 'swc_node_id', 'parent_id', 'node_type',
        'x_px', 'y_px', 'x_phys', 'y_phys', 'z_phys',
        'distance_from_soma', 'path_length_from_soma',
        'branch_order', 'is_tip', 'valid_point',
    ]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')

    all_fields = [c for c in df.columns if c not in (
        numeric_cols + ['case_id', 'swc_frame', 'swc_node_id']
    )]
    for f in all_fields:
        df[f] = pd.to_numeric(df[f], errors='coerce')

    # Keep only valid points
    if 'valid_point' in df.columns:
        n_invalid = (df['valid_point'] == 0).sum()
        if n_invalid:
            print(f"[info] Dropping {n_invalid} invalid-probe rows (valid_point=0)")
        df = df[df['valid_point'] == 1].copy()

    return df


def save_fig(fig, out_dir, name, dpi=150):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'svg'):
        fpath = out_dir / f"{name}.{ext}"
        fig.savefig(fpath, dpi=dpi, bbox_inches='tight')
        print(f"  saved: {fpath}")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: Concentration vs distance_from_soma
# ─────────────────────────────────────────────────────────────────────────────

def fig_concentration_vs_distance(df, fields, time_steps, out_dir, case_id):
    """
    For each field: scatter of concentration vs. path_length_from_soma,
    coloured by time step. Tip nodes are shown with larger markers.
    """
    available_ts = sorted(df['time_step'].unique())
    if time_steps is None:
        # Auto-select ~5 evenly spaced time steps
        n = min(5, len(available_ts))
        idxs = np.linspace(0, len(available_ts) - 1, n, dtype=int)
        time_steps = [available_ts[i] for i in idxs]

    n_fields = len(fields)
    ncols = min(3, n_fields)
    nrows = (n_fields + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    fig.suptitle(f'{case_id} — Concentration vs. path length from soma', fontsize=13)

    try:
        cmap = matplotlib.colormaps.get_cmap('viridis').resampled(len(time_steps))
    except AttributeError:
        cmap = cm.get_cmap('viridis', len(time_steps))
    ts_colors = {ts: cmap(i) for i, ts in enumerate(time_steps)}

    for fi, fname in enumerate(fields):
        ax = axes[fi // ncols][fi % ncols]
        if fname not in df.columns:
            ax.set_visible(False)
            continue

        for ts in time_steps:
            sub = df[df['time_step'] == ts]
            if sub.empty:
                continue
            # Regular nodes
            reg = sub[sub['is_tip'] == 0]
            tip = sub[sub['is_tip'] == 1]
            c = ts_colors[ts]
            ax.scatter(reg['path_length_from_soma'], reg[fname],
                       s=18, alpha=0.7, color=c, label=f't={ts}', zorder=2)
            ax.scatter(tip['path_length_from_soma'], tip[fname],
                       s=60, alpha=0.9, color=c, marker='^', zorder=3)

        ax.set_xlabel('Path length from soma')
        ax.set_ylabel(label(fname))
        ax.set_title(label(fname), fontsize=10)
        ax.grid(True, alpha=0.3)

    # Single legend for time steps
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker='o', color=ts_colors[ts],
                       linestyle='', label=f't = {ts}')
               for ts in time_steps]
    handles.append(Line2D([0], [0], marker='^', color='gray',
                           linestyle='', label='tip node', markersize=8))
    axes[0][-1].legend(handles=handles, fontsize=8, loc='upper right')

    # Hide unused subplots
    for fi in range(n_fields, nrows * ncols):
        axes[fi // ncols][fi % ncols].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, out_dir, 'fig1_concentration_vs_distance', dpi=150)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Tip concentrations over time
# ─────────────────────────────────────────────────────────────────────────────

def fig_tip_concentrations_over_time(df, fields, out_dir, case_id):
    """
    For each field: concentration at tip nodes over time.
    Each tip is a separate line, coloured by path_length_from_soma.
    """
    tips = df[df['is_tip'] == 1].copy()
    if tips.empty:
        print("[warn] No tip nodes found in data — skipping fig 2")
        return

    # Identify individual tips by their (path_length range, branch_order) at last frame
    last_ts = tips['time_step'].max()
    tip_ids = (
        tips[tips['time_step'] == last_ts][['swc_node_id', 'path_length_from_soma']]
        .sort_values('path_length_from_soma')
        .reset_index(drop=True)
    )

    n_fields = len(fields)
    ncols = min(3, n_fields)
    nrows = (n_fields + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    fig.suptitle(f'{case_id} — Tip node concentrations over time', fontsize=13)

    # colour tips by their final path length
    max_path = tip_ids['path_length_from_soma'].max()
    try:
        cmap = matplotlib.colormaps.get_cmap('plasma')
    except AttributeError:
        cmap = cm.get_cmap('plasma')

    for fi, fname in enumerate(fields):
        ax = axes[fi // ncols][fi % ncols]
        if fname not in df.columns:
            ax.set_visible(False)
            continue

        for _, row in tip_ids.iterrows():
            # Find rows for this node across time (nodes may be renumbered each frame;
            # track by proximity of path_length)
            pl = row['path_length_from_soma']
            tol = max_path * 0.05  # 5 % path length tolerance for tip tracking
            tip_ts = tips[
                (tips['path_length_from_soma'] >= pl - tol) &
                (tips['path_length_from_soma'] <= pl + tol)
            ].groupby('time_step')[fname].median().reset_index()

            if tip_ts.empty or len(tip_ts) < 2:
                continue

            norm_pl = pl / max_path if max_path > 0 else 0
            c = cmap(norm_pl)
            ax.plot(tip_ts['time_step'], tip_ts[fname],
                    color=c, alpha=0.8, linewidth=1.2)

        ax.set_xlabel('Time step')
        ax.set_ylabel(label(fname))
        ax.set_title(label(fname), fontsize=10)
        ax.grid(True, alpha=0.3)

    # Colorbar for path length
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=plt.Normalize(vmin=0, vmax=max_path))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.6, pad=0.02)
    cbar.set_label('Path length from soma at final frame')

    for fi in range(n_fields, nrows * ncols):
        axes[fi // ncols][fi % ncols].set_visible(False)

    fig.tight_layout(rect=[0, 0, 0.90, 0.97])
    save_fig(fig, out_dir, 'fig2_tip_concentrations_over_time', dpi=150)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: Heatmap — path length × time
# ─────────────────────────────────────────────────────────────────────────────

def fig_heatmap_path_vs_time(df, fields, out_dir, case_id, n_bins=30):
    """
    2D heatmap: x = path_length_from_soma (binned), y = time_step,
    colour = median concentration.
    Shown for fields specified (default: u_t and u_ca_cyt if present).
    """
    heatmap_fields = [f for f in fields if f in df.columns]
    if not heatmap_fields:
        return

    max_path = df['path_length_from_soma'].max()
    bins = np.linspace(0, max_path, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    df = df.copy()
    df['path_bin'] = pd.cut(df['path_length_from_soma'], bins=bins,
                             labels=bin_centers, include_lowest=True)
    df['path_bin'] = df['path_bin'].astype(float)

    time_steps = sorted(df['time_step'].unique())

    n_fields = len(heatmap_fields)
    fig, axes = plt.subplots(1, n_fields, figsize=(7 * n_fields, 6), squeeze=False)
    fig.suptitle(f'{case_id} — Concentration heatmap (path length × time)', fontsize=13)

    for fi, fname in enumerate(heatmap_fields):
        ax = axes[0][fi]
        pivot = (
            df.groupby(['time_step', 'path_bin'])[fname]
            .median()
            .unstack('path_bin')
            .reindex(time_steps)
        )
        im = ax.imshow(
            pivot.values,
            aspect='auto', origin='lower',
            extent=[0, max_path, time_steps[0] - 0.5, time_steps[-1] + 0.5],
            cmap='inferno',
        )
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_xlabel('Path length from soma')
        ax.set_ylabel('Time step')
        ax.set_title(label(fname), fontsize=11)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, out_dir, 'fig3_heatmap_path_vs_time', dpi=150)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: Field distributions over time (violin)
# ─────────────────────────────────────────────────────────────────────────────

def fig_field_distributions(df, fields, out_dir, case_id):
    """
    Violin plot of each field's distribution across all nodes at each time step.
    """
    available = [f for f in fields if f in df.columns]
    if not available:
        return

    time_steps = sorted(df['time_step'].unique())
    # Subsample time steps for readability if too many
    if len(time_steps) > 15:
        idxs = np.linspace(0, len(time_steps) - 1, 15, dtype=int)
        time_steps = [time_steps[i] for i in idxs]

    n_fields = len(available)
    ncols = min(3, n_fields)
    nrows = (n_fields + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    fig.suptitle(f'{case_id} — Field distributions over time', fontsize=13)

    for fi, fname in enumerate(available):
        ax = axes[fi // ncols][fi % ncols]
        data_by_ts = []
        valid_ts = []
        for ts in time_steps:
            vals = df[df['time_step'] == ts][fname].dropna().values
            if len(vals) >= 2:
                data_by_ts.append(vals)
                valid_ts.append(ts)

        if not data_by_ts:
            ax.set_visible(False)
            continue

        parts = ax.violinplot(data_by_ts, positions=valid_ts,
                               showmedians=True, widths=0.8)
        for pc in parts['bodies']:
            pc.set_facecolor(color(fname))
            pc.set_alpha(0.6)
        for key in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            if key in parts:
                parts[key].set_color('black')
                parts[key].set_linewidth(0.8)

        ax.set_xlabel('Time step')
        ax.set_ylabel(label(fname))
        ax.set_title(label(fname), fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')

    for fi in range(n_fields, nrows * ncols):
        axes[fi // ncols][fi % ncols].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save_fig(fig, out_dir, 'fig4_field_distributions', dpi=150)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--csv', required=True,
                        help='Path to the CSV produced by sample_fields_on_swc.py')
    parser.add_argument('--out-dir', default='./swc_field_figures',
                        help='Output directory for figures (default: ./swc_field_figures)')
    parser.add_argument('--case-id', default='',
                        help='Case label for figure titles (auto-inferred from CSV if omitted)')
    parser.add_argument('--fields', default='u_t,u_u,u_b,u_p,u_ca_cyt,inhibitor',
                        help='Comma-separated list of fields to plot '
                             '(default: u_t,u_u,u_b,u_p,u_ca_cyt,inhibitor)')
    parser.add_argument('--time-steps', default=None,
                        help='Comma-separated list of time step indices for fig 1 '
                             '(default: auto-select 5 evenly spaced)')
    parser.add_argument('--heatmap-fields', default='u_t,u_ca_cyt',
                        help='Fields to show in heatmap (fig 3); '
                             'default: u_t,u_ca_cyt')
    parser.add_argument('--n-bins', type=int, default=30,
                        help='Number of path-length bins for heatmap (default: 30)')
    args = parser.parse_args()

    print(f"[info] Loading {args.csv} ...")
    df = load_csv(args.csv)
    print(f"[info] {len(df)} rows, {df['time_step'].nunique()} time steps, "
          f"{df['swc_frame'].nunique()} frames")

    case_id = args.case_id or (df['case_id'].iloc[0] if 'case_id' in df.columns else 'case')
    fields = [f.strip() for f in args.fields.split(',') if f.strip()]
    heatmap_fields = [f.strip() for f in args.heatmap_fields.split(',') if f.strip()]

    time_steps = None
    if args.time_steps:
        time_steps = [int(t.strip()) for t in args.time_steps.split(',')]

    print(f"\n[fig 1] Concentration vs. distance from soma ...")
    fig_concentration_vs_distance(df, fields, time_steps, args.out_dir, case_id)

    print(f"[fig 2] Tip concentrations over time ...")
    fig_tip_concentrations_over_time(df, fields, args.out_dir, case_id)

    print(f"[fig 3] Heatmap (path × time) ...")
    fig_heatmap_path_vs_time(df, heatmap_fields, args.out_dir, case_id, n_bins=args.n_bins)

    print(f"[fig 4] Field distributions over time ...")
    fig_field_distributions(df, fields, args.out_dir, case_id)

    print(f"\n[done] All figures written to {args.out_dir}/")


if __name__ == '__main__':
    main()
