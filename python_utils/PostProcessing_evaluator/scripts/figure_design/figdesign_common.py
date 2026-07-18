"""
figdesign_common.py
===================
Shared helpers for the second-generation figure-design iteration.

These figures are CANDIDATE DESIGNS for discussion with Nicole. They are
NOT publication-final. Each figure is saved with a banner footer that
includes timestamp, source CSV, and an "exploratory" tag, so prints can
never be mistaken for headline figures.
"""

from pathlib import Path
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FIELD_LABELS = {
    'u_t':       'Tubulin (u_t)',
    'u_u':       'MAP2 free (u_u)',
    'u_b':       'MAP2 bound (u_b)',
    'u_p':       'MAP2 phospho (u_p)',
    'u_ca_cyt':  'Cytosolic Ca (u_ca_cyt)',
    'inhibitor': 'Inhibitor',
    'lsf':       'Level-set (lsf)',
}

FIELD_COLORS = {
    'u_t':       '#d62728',
    'u_u':       '#2ca02c',
    'u_b':       '#1f77b4',
    'u_p':       '#ff7f0e',
    'u_ca_cyt':  '#9467bd',
    'inhibitor': '#17becf',
    'lsf':       '#7f7f7f',
}

MAIN_FIELDS = ['u_t', 'u_u', 'u_b', 'u_p', 'u_ca_cyt', 'inhibitor']
MAP2_GROUP  = ['u_u', 'u_b', 'u_p']

EXPLORATORY_FOOTER = ('Second-generation candidate design — exploratory; '
                      'not publication-final. PostProcessing_evaluator.')


def load_enriched(csv_path):
    df = pd.read_csv(csv_path)
    if 'valid_point' in df.columns:
        df = df[df['valid_point'] != 0].copy()
    for c in df.columns:
        if c in ('case_id', 'node_type'):
            continue
        df[c] = pd.to_numeric(df[c], errors='ignore')
    return df


def detect_meaningful_fields(df, fields, min_relative_range=1e-4):
    """
    Return a list of fields whose dynamic range across all rows is at least
    `min_relative_range * max(|values|)`. Used to skip flat/zero fields.
    Returns (kept_fields, dropped_fields, reasons).
    """
    kept, dropped, reasons = [], [], {}
    for f in fields:
        if f not in df.columns:
            dropped.append(f); reasons[f] = 'column missing'
            continue
        v = pd.to_numeric(df[f], errors='coerce').dropna().values
        if v.size == 0:
            dropped.append(f); reasons[f] = 'no finite values'
            continue
        rng = float(np.nanmax(v) - np.nanmin(v))
        scale = float(np.nanmax(np.abs(v)) + 1e-300)
        if rng <= min_relative_range * scale:
            dropped.append(f); reasons[f] = f'flat: range={rng:.3g} scale={scale:.3g}'
            continue
        kept.append(f)
    return kept, dropped, reasons


def add_exploratory_footer(fig, csv_path, extra=''):
    txt = EXPLORATORY_FOOTER + (f'  ·  {extra}' if extra else '') + \
          f'\nsource CSV: {csv_path}'
    fig.text(0.02, 0.005, txt, ha='left', va='bottom',
             fontsize=6, color='#555555', style='italic')


def save_figure(fig, out_dir, stem, dpi=160):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'svg'):
        p = out_dir / f'{stem}.{ext}'
        fig.savefig(p, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return [out_dir / f'{stem}.png', out_dir / f'{stem}.svg']


def lmax_per_step(df):
    """Return per-step Lmax(t) Series indexed by time_step (sorted ascending)."""
    g = df.groupby('time_step')['path_length_from_soma'].max()
    return g.sort_index()


def l0_value(df):
    """Initial maximum path length at the earliest available timestep."""
    if 'L0' in df.columns:
        v = df['L0'].dropna()
        if len(v):
            return float(v.iloc[0])
    return float(df.loc[df['time_step'] == df['time_step'].min(),
                        'path_length_from_soma'].max())


def time_values(df):
    """Return per-step time value (1:1 with time_step)."""
    if 'time' in df.columns:
        g = df.groupby('time_step')['time'].first().sort_index()
        return g
    return pd.Series(sorted(df['time_step'].unique()),
                     index=sorted(df['time_step'].unique()))


def percentile_clip(values, lo=2.0, hi=98.0):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0, 1.0
    return float(np.percentile(v, lo)), float(np.percentile(v, hi))
