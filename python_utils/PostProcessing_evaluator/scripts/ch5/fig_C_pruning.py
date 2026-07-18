#!/usr/bin/env python3
"""
fig_C_pruning.py — CH5 Part C (pruning regime).

§5.9 Kaplan–Meier  (branch survival).
§5.10 ROC of inhibitor as a pruning classifier.
§5.11 Event-triggered concentration average.

CAVEATS specific to D7.5_V3:
  - inhibitor ≈ 0 across the whole simulation; the ROC will likely give
    AUC ≈ 0.5 — the figure is scaffolding for active-pruning datasets.
  - branch_id is RE-ASSIGNED per frame by enrich_swc_samples.py — branches
    are not tracked across timesteps. The KM curve below uses
    elongation_rate < -ε for ≥ 3 consecutive frames as the pruning event,
    interpreted on the GLOBAL longest path (Lmax) rather than per-branch.
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
from sklearn.metrics import roc_curve, roc_auc_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from utils_ch5 import (  # noqa: E402
    DT_SECONDS, add_exploratory_footer, gc_timeseries, load_enriched,
    per_step_value, tip_velocity_proxy,
)


def detect_pruning_events(elong: pd.Series, eps: float = 0.005,
                            min_run: int = 3) -> list:
    """
    Return list of (start_step, end_step) tuples where elongation_rate
    stays below -eps for >= min_run consecutive frames.
    """
    e = elong.sort_index()
    steps = list(e.index)
    vals = e.values
    in_run = False
    run_start = None
    events = []
    for i, v in enumerate(vals):
        if v is None or pd.isna(v):
            in_run = False; run_start = None; continue
        if v < -eps:
            if not in_run:
                run_start = steps[i]; in_run = True
        else:
            if in_run:
                run_end = steps[i - 1]
                if (run_end - run_start + 1) >= min_run:
                    events.append((int(run_start), int(run_end)))
                in_run = False; run_start = None
    if in_run:
        run_end = steps[-1]
        if (run_end - run_start + 1) >= min_run:
            events.append((int(run_start), int(run_end)))
    return events


def figure_km_pruning(df, out_dir, csv_path):
    """
    §5.9 KM: treat the whole simulation as one "branch" (since branch_id
    isn't tracked across time). Define survival as "Lmax has not yet
    plateaued below 95% of its peak rate". Right-censor at final frame.

    This is a coarse proxy; the figure is primarily structural to show
    the KM layout for later datasets that have true pruning events.
    """
    elong = tip_velocity_proxy(df)
    events = detect_pruning_events(elong, eps=0.005, min_run=3)
    last_t = int(elong.index.max())
    if events:
        first_event_start = events[0][0]
    else:
        first_event_start = None
    # Single-branch KM (one "subject"; nothing useful — render as scaffolding)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))

    # (a) KM-style survival proxy
    times = np.array(sorted(elong.index))
    if first_event_start is None:
        S = np.ones_like(times, dtype=float)
    else:
        S = np.where(times < first_event_start, 1.0, 0.0)
    axes[0].step(times, S, where='post', color='#1f77b4', linewidth=2)
    axes[0].set_xlabel('time step')
    axes[0].set_ylabel('S(t)  (KM-style survival, single subject)')
    axes[0].set_title('§5.9 Kaplan–Meier (scaffolding; 1 subject in this dataset)')
    axes[0].set_ylim(-0.05, 1.05)
    axes[0].grid(True, alpha=0.25)
    txt = (f'pruning runs detected (rate < -0.005, ≥3 frames):\n  {events}\n'
           f'(branch_id NOT tracked across frames in enriched CSV;\n'
           ' true per-branch KM requires upstream tracking)')
    axes[0].text(0.02, 0.45, txt, transform=axes[0].transAxes,
                  ha='left', va='top', fontsize=8,
                  bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                            edgecolor='#666', alpha=0.85))

    # (b) elongation rate timeline with detected events shaded
    axes[1].plot(times, elong.values, color='#cc0066', linewidth=1.6,
                  marker='o', markersize=4)
    axes[1].axhline(0, color='#888', linewidth=0.8)
    axes[1].axhline(-0.005, color='#888', linewidth=0.6, linestyle=':')
    for s, e in events:
        axes[1].axvspan(s - 0.4, e + 0.4, color='#cc0066', alpha=0.18)
    axes[1].set_xlabel('time step')
    axes[1].set_ylabel('elongation_rate (µm/frame)')
    axes[1].set_title('detected pruning runs (shaded)')
    axes[1].grid(True, alpha=0.25)

    fig.suptitle('CH5 §5.9 — Kaplan–Meier scaffolding '
                  '(branch_id not tracked; structural figure only)',
                  fontsize=11)
    add_exploratory_footer(fig, csv_path,
                            extra='single-subject KM; events from rate timeline')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_C_kaplan_meier.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)
    return dict(events=events, n_events=len(events))


def figure_roc_inhibitor(df, out_dir, csv_path, pre_window=3):
    """
    §5.10 ROC — label = 1 if a step is within pre_window of a pruning event,
    else 0; score = mean inhibitor at GC nodes per step.

    With inhibitor ≈ 0 across the dataset, AUC is expected near 0.5.
    """
    elong = tip_velocity_proxy(df)
    events = detect_pruning_events(elong, eps=0.005, min_run=3)
    inhib_gc = gc_timeseries(df, 'inhibitor')['mean']
    steps = sorted(set(elong.index) & set(inhib_gc.index))
    labels = np.zeros(len(steps), dtype=int)
    for i, t in enumerate(steps):
        for s, e in events:
            if s - pre_window <= t < s:
                labels[i] = 1; break
    scores = inhib_gc.reindex(steps).fillna(0.0).values

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    auc = float('nan'); fpr = []; tpr = []
    if labels.sum() > 0 and labels.sum() < len(labels):
        fpr, tpr, _ = roc_curve(labels, scores)
        auc = roc_auc_score(labels, scores)
        ax.plot(fpr, tpr, color='#cc0066', linewidth=2,
                 label=f'ROC (AUC = {auc:.3f})')
    ax.plot([0, 1], [0, 1], color='#888', linewidth=1.0, linestyle='--',
             label='chance (AUC = 0.5)')
    ax.set_xlabel('false positive rate')
    ax.set_ylabel('true positive rate')
    ax.set_title('§5.10 ROC: inhibitor as pruning classifier (pre-window = '
                  f'{pre_window} frames)')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.25)
    n_pos = int(labels.sum()); n_neg = int(len(labels) - n_pos)
    txt = (f'positive frames: {n_pos}   negative: {n_neg}\n'
           f'inhibitor range: [{scores.min():.3g}, {scores.max():.3g}]')
    if not (labels.sum() > 0 and labels.sum() < len(labels)):
        txt += '\n(no class variation — ROC undefined)'
    ax.text(0.02, 0.98, txt, transform=ax.transAxes,
             ha='left', va='top', fontsize=8,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       edgecolor='#666', alpha=0.85))

    fig.suptitle('CH5 §5.10 — ROC (inhibitor is effectively zero in this '
                  'dataset; figure is scaffolding)', fontsize=10)
    add_exploratory_footer(fig, csv_path,
                            extra=f'pre-window = {pre_window} frames; '
                                  f'AUC near 0.5 expected when I is flat')
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_C_roc_inhibitor_pruning.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)
    return dict(auc=auc, n_pos=n_pos, n_neg=n_neg, events=events)


def figure_event_triggered(df, out_dir, csv_path, window=5):
    """
    §5.11 event-triggered average over GC concentrations, aligned to
    detected pruning events. With ≤2 events in D7.5_V3 the CI is wide;
    we still produce the layout.
    """
    elong = tip_velocity_proxy(df)
    events = detect_pruning_events(elong, eps=0.005, min_run=3)
    if not events:
        # No events — emit a placeholder figure and return early.
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5,
                 'no pruning events detected\n(scaffolding figure)',
                 ha='center', va='center', transform=ax.transAxes,
                 fontsize=14, color='#666')
        ax.set_xticks([]); ax.set_yticks([])
        add_exploratory_footer(fig, csv_path,
                                extra='no events; placeholder')
        out = Path(out_dir)
        for ext in ('png', 'svg'):
            fig.savefig(out / f'fig_C_event_triggered_average.{ext}',
                        dpi=160, bbox_inches='tight')
        plt.close(fig)
        return dict(n_events=0)

    species = [('u_t', '#d62728'), ('u_b', '#1f77b4'),
                ('u_ca_cyt', '#9467bd'), ('inhibitor', '#17becf')]
    series_by_sp = {}
    for sp, _ in species:
        series_by_sp[sp] = gc_timeseries(df, sp)['mean']

    lags = np.arange(-window, window + 1)
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    for sp, color in species:
        s = series_by_sp[sp]
        # z-score across the full series
        if s.std() > 1e-12:
            zs = (s - s.mean()) / s.std()
        else:
            zs = s * 0.0
        # collect aligned slices per event
        slices = []
        for ev_start, _ in events:
            row = np.full(len(lags), np.nan)
            for k, lag in enumerate(lags):
                t_target = ev_start + lag
                if t_target in zs.index:
                    row[k] = zs.loc[t_target]
            slices.append(row)
        slices = np.array(slices)
        mean = np.nanmean(slices, axis=0)
        # bootstrap percentile CI across events (very wide for small N)
        n_boot = 400; rng = np.random.default_rng(0)
        boot = np.empty((n_boot, len(lags)))
        for b in range(n_boot):
            idx = rng.integers(0, len(slices), size=len(slices))
            boot[b] = np.nanmean(slices[idx], axis=0)
        lo = np.nanpercentile(boot, 2.5, axis=0)
        hi = np.nanpercentile(boot, 97.5, axis=0)
        ax.plot(lags * DT_SECONDS, mean, color=color, linewidth=2,
                 marker='o', markersize=4, label=sp)
        ax.fill_between(lags * DT_SECONDS, lo, hi, color=color, alpha=0.18)

    ax.axvline(0, color='#888', linewidth=1.0, linestyle='--')
    ax.set_xlabel('Δt  (s; 0 = pruning-event start)')
    ax.set_ylabel('z-scored GC concentration')
    ax.set_title(f'§5.11 event-triggered average  ({len(events)} event(s))')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.25)
    add_exploratory_footer(fig, csv_path,
                            extra=('bootstrap CI across events; '
                                   'N events small ⇒ CI is wide'))
    fig.tight_layout()
    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_C_event_triggered_average.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)
    return dict(n_events=len(events), events=events)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    df = load_enriched(args.csv)
    km   = figure_km_pruning(df, out_dir, args.csv)
    roc  = figure_roc_inhibitor(df, out_dir, args.csv)
    eta  = figure_event_triggered(df, out_dir, args.csv)
    with open(out_dir / 'fig_C_summary.json', 'w') as fh:
        json.dump(dict(csv=args.csv, km=km, roc=roc, event_triggered=eta),
                  fh, indent=2, default=float)
    print(f'fig_C done. summary -> {out_dir / "fig_C_summary.json"}')


if __name__ == '__main__':
    main()
