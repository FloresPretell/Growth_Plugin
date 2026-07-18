"""
utils_ch5.py — shared helpers for CH5 post-analysis figures.

Implements only what is needed by the fig_A through fig_D modules:
loaders, GC-region time-series, stationary block bootstrap, plateau
detection, and the §5.12 spatial-correlation per-arclength-bin helper.

No heavy computation lives here; modules import these and orchestrate.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# -----------------------------------------------------------------------------
# Model parameters known from the Lua control file
# -----------------------------------------------------------------------------

# /project/.../Scripts_lua/test_parameter_fail.lua  (lines 80, 110, 138)
#   elongationRate = 225  (nondimensional)
#   = 225 × (10 µm)/(500 s · 40 µM · 6.25 µM) = 0.018 µm/(s · µM²)
ALPHA_MODEL = 0.018           # µm / (s · µM²)
DELTA_MODEL = 0.0             # retraction baseline (none in this run)
DT_SECONDS  = 6.0             # frame interval

# No I_avoid / I_prune thresholds appear in Scripts_lua/. The inhibitor field
# is essentially zero (max ≲ 5e-7) across D7.5_V3, so a sensible operational
# threshold is just "any nonzero inhibitor" or "above the 95-percentile".
# Modules will detect this and emit a clear caption if applicable.


EXPLORATORY_FOOTER = (
    'CH5 post-analysis — exploratory; not publication-final. '
    'PostProcessing_evaluator.')


# -----------------------------------------------------------------------------
# Loaders
# -----------------------------------------------------------------------------

def load_enriched(csv_path: str) -> pd.DataFrame:
    """Load enriched CSV; cast numerics; filter valid_point if present."""
    df = pd.read_csv(csv_path)
    if 'valid_point' in df.columns:
        df = df[df['valid_point'] != 0].copy()
    for c in df.columns:
        if c in ('case_id', 'node_type'):
            continue
        # pd.to_numeric with errors='ignore' is deprecated; do it manually
        try:
            df[c] = pd.to_numeric(df[c])
        except (ValueError, TypeError):
            pass
    return df


def detect_active_elongation_phase(elong: pd.Series,
                                   threshold: float = 0.03) -> list:
    """
    Return the list of time_steps that constitute the 'active elongation phase':
    the contiguous initial run where elongation_rate > threshold. Stops at the
    first step where the rate falls below threshold for ≥1 step.
    """
    e = elong.sort_index()
    active = []
    for t, v in e.items():
        if pd.isna(v) or v <= threshold:
            break
        active.append(int(t))
    return active


# -----------------------------------------------------------------------------
# Time-series helpers
# -----------------------------------------------------------------------------

def gc_timeseries(df: pd.DataFrame, field: str) -> pd.DataFrame:
    """
    Mean / SD / N of `field` over growth-cone-region nodes per time_step.
    Returns DataFrame indexed by time_step with columns ['mean', 'std', 'n'].
    """
    gc = df[df['is_growth_cone_region'] == 1]
    g = gc.groupby('time_step')[field].agg(['mean', 'std', 'count'])
    g = g.rename(columns={'count': 'n'}).sort_index()
    g['std'] = g['std'].fillna(0.0)
    return g


def tip_timeseries(df: pd.DataFrame, field: str) -> pd.DataFrame:
    """
    Like gc_timeseries but restricted to nodes with `is_tip == 1`.

    Use this for fields that are only well-defined at the moving interface
    (e.g. `norm_vel_int`, `curvature`, `GrowthVel`). For diffusive
    concentration fields, gc_timeseries is fine.
    """
    tips = df[df['is_tip'] == 1]
    if len(tips) == 0:
        return pd.DataFrame(columns=['mean', 'std', 'n'])
    g = tips.groupby('time_step')[field].agg(['mean', 'std', 'count'])
    g = g.rename(columns={'count': 'n'}).sort_index()
    g['std'] = g['std'].fillna(0.0)
    return g


def per_step_value(df: pd.DataFrame, field: str) -> pd.Series:
    """Return a Series of `field` indexed by time_step (first value per step)."""
    return df.groupby('time_step')[field].first().sort_index()


def tip_velocity_proxy(df: pd.DataFrame) -> pd.Series:
    """v_tip(t) ≈ elongation_rate (already dLmax/dt per step)."""
    return per_step_value(df, 'elongation_rate')


def lmax_per_step(df: pd.DataFrame) -> pd.Series:
    return df.groupby('time_step')['path_length_from_soma'].max().sort_index()


def time_values(df: pd.DataFrame) -> pd.Series:
    if 'time' in df.columns:
        return per_step_value(df, 'time')
    return pd.Series(sorted(df['time_step'].unique()),
                     index=sorted(df['time_step'].unique()))


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------

def stationary_block_bootstrap_corr(x, y, n_boot=1000, block_len=5,
                                    ci=0.95, seed=0):
    """
    Politis–Romano stationary block bootstrap CI for Pearson correlation.

    x, y: 1-D arrays of equal length.
    block_len: mean block length b̄ = 1/p; should be a few × τ_int.
    Returns (rho_point, rho_lo, rho_hi).
    """
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    n = len(x)
    if n < 4:
        return (float('nan'), float('nan'), float('nan'))
    rng = np.random.default_rng(seed)
    p = 1.0 / max(block_len, 1)
    rhos = np.empty(n_boot)
    for b in range(n_boot):
        idx = np.empty(n, dtype=int)
        i = 0
        cur = int(rng.integers(0, n))
        while i < n:
            idx[i] = cur
            i += 1
            if rng.random() < p:
                cur = int(rng.integers(0, n))
            else:
                cur = (cur + 1) % n
        xb = x[idx]; yb = y[idx]
        sx = xb.std(); sy = yb.std()
        if sx < 1e-30 or sy < 1e-30:
            rhos[b] = np.nan
        else:
            rhos[b] = ((xb - xb.mean()) * (yb - yb.mean())).mean() / (sx * sy)
    rho_point = ((x - x.mean()) * (y - y.mean())).mean() / (x.std() * y.std())
    alpha = (1 - ci) / 2
    lo, hi = np.nanpercentile(rhos, [100 * alpha, 100 * (1 - alpha)])
    return float(rho_point), float(lo), float(hi)


def cross_correlation(x: np.ndarray, y: np.ndarray, max_lag: int) -> np.ndarray:
    """ρ(τ) = corr(x(t), y(t+τ)) for τ in [-max_lag, +max_lag]."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    n = len(x)
    out = np.full(2 * max_lag + 1, np.nan)
    for k, tau in enumerate(range(-max_lag, max_lag + 1)):
        if tau >= 0:
            xs = x[:n - tau]; ys = y[tau:]
        else:
            xs = x[-tau:]; ys = y[:n + tau]
        m = np.isfinite(xs) & np.isfinite(ys)
        if m.sum() < 4: continue
        xs = xs[m]; ys = ys[m]
        sx = xs.std(); sy = ys.std()
        if sx < 1e-30 or sy < 1e-30: continue
        out[k] = ((xs - xs.mean()) * (ys - ys.mean())).mean() / (sx * sy)
    return out


def integrated_autocorr_time(series: np.ndarray, max_lag: int = None) -> float:
    """
    τ_int = 1 + 2 Σ_{k=1}^{cut} r_k with Sokal's automated windowing
    (cut = smallest k where k ≥ 5·τ_running).
    """
    x = np.asarray(series, dtype=float)
    x = x[np.isfinite(x)]
    n = len(x)
    if n < 4:
        return 1.0
    x = x - x.mean()
    var = x.var()
    if var < 1e-30:
        return 1.0
    if max_lag is None:
        max_lag = min(n - 2, 100)
    tau = 1.0
    for k in range(1, max_lag + 1):
        rk = (x[:n - k] * x[k:]).mean() / var
        tau += 2 * rk
        if k >= 5 * tau:
            break
    return max(1.0, tau)


def hac_se_ols(x: np.ndarray, y: np.ndarray, max_lag: int = None):
    """
    OLS y = β₀ + β₁ x + ε  with HAC (Newey–West) standard errors.
    Returns (slope, intercept, se_slope, se_intercept, r2, n).
    """
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    n = len(x)
    if n < 4:
        return (np.nan,) * 4 + (np.nan, n)
    X = np.column_stack([np.ones(n), x])
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    ss_res = (resid ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    if max_lag is None:
        max_lag = max(1, int(np.floor(4 * (n / 100) ** (2 / 9))))
    # Newey-West Bartlett kernel
    u = resid
    S = np.zeros((2, 2))
    for k in range(0, max_lag + 1):
        w = 1.0 - k / (max_lag + 1)
        gamma = np.zeros((2, 2))
        for t in range(k, n):
            xt = X[t]; xtk = X[t - k]
            ut = u[t]; utk = u[t - k]
            gamma += np.outer(xt * ut, xtk * utk)
        if k == 0:
            S += w * gamma
        else:
            S += w * (gamma + gamma.T)
    cov = XtX_inv @ S @ XtX_inv
    se = np.sqrt(np.diag(cov))
    return float(beta[1]), float(beta[0]), float(se[1]), float(se[0]), \
           float(r2), int(n)


def detect_plateau(series: pd.Series, eps: float = 0.1):
    """Return time_step at which |dL/dt| first falls below eps·peak."""
    v = series.dropna()
    if len(v) == 0:
        return None
    a = np.abs(v.values)
    peak = a.max()
    if peak <= 0:
        return None
    active = np.where(a > eps * peak)[0]
    if len(active) == 0:
        return None
    last_active = int(v.index[active[-1]])
    return last_active


# -----------------------------------------------------------------------------
# Spatial correlation per arclength bin (§5.12)
# -----------------------------------------------------------------------------

def spatial_corr_per_bin(df: pd.DataFrame, field_x: str, field_y: str,
                         n_bins: int = 12) -> tuple:
    """
    For each arclength bin s, compute Pearson rho of (mean field_x at s, t)
    vs (mean field_y at s, t) over time. Returns (centers, rho_array, n_array).
    """
    pmax = float(df['path_length_from_soma'].max())
    edges = np.linspace(0.0, pmax, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    steps = sorted(df['time_step'].unique())
    rho = np.full(n_bins, np.nan)
    nused = np.zeros(n_bins, dtype=int)
    for b in range(n_bins):
        xseries = []; yseries = []
        for t in steps:
            sub = df[(df['time_step'] == t) &
                     (df['path_length_from_soma'] >= edges[b]) &
                     (df['path_length_from_soma'] <  edges[b + 1])]
            if len(sub) == 0:
                xseries.append(np.nan); yseries.append(np.nan); continue
            xseries.append(pd.to_numeric(sub[field_x], errors='coerce').mean())
            yseries.append(pd.to_numeric(sub[field_y], errors='coerce').mean())
        x = np.array(xseries); y = np.array(yseries)
        m = np.isfinite(x) & np.isfinite(y)
        nused[b] = int(m.sum())
        if nused[b] < 4: continue
        sx = x[m].std(); sy = y[m].std()
        if sx < 1e-30 or sy < 1e-30: continue
        rho[b] = ((x[m] - x[m].mean()) * (y[m] - y[m].mean())).mean() / (sx * sy)
    return centers, rho, nused


# -----------------------------------------------------------------------------
# Adaptive helpers (case-agnostic; never hard-code physical ranges)
# -----------------------------------------------------------------------------

def make_correlation_norm(rho_values, pad: float = 0.05):
    """
    Build a TwoSlopeNorm for correlation data that is ALWAYS anchored at 0
    and uses symmetric vmin = -vmax so blue and red are comparable.

    Adapts to the actual data range with a small pad; clips to [-1, +1].
    Works whether the data spans [-0.1, +0.9] or [-0.8, +0.8].
    """
    from matplotlib.colors import TwoSlopeNorm
    import numpy as np
    finite = np.asarray(rho_values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return TwoSlopeNorm(vcenter=0.0, vmin=-1.0, vmax=1.0)
    half = max(abs(float(finite.min())), abs(float(finite.max()))) + pad
    half = min(max(half, 0.05), 1.0)
    return TwoSlopeNorm(vcenter=0.0, vmin=-half, vmax=half)


def adaptive_n_bins(values, lo: int = 10, hi: int = 50):
    """
    Sqrt rule with floor/ceil. Returns an int in [lo, hi].
    """
    import numpy as np
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return lo
    return int(max(lo, min(hi, round(v.size ** 0.5))))


def adaptive_percentile_clip(values, lo: float = 2.0, hi: float = 98.0,
                              pad_frac: float = 0.05):
    """
    Return (vmin, vmax) clipped to percentile range with a small pad.
    Replaces ad-hoc `vmin=0`, `vmax=ymax` constructions.
    """
    import numpy as np
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0, 1.0
    a, b = float(np.percentile(v, lo)), float(np.percentile(v, hi))
    pad = (b - a) * pad_frac
    return a - pad, b + pad


def adaptive_tip_band(df, coverage_nodes: int = 5,
                      min_fraction_of_range: float = 0.15,
                      absolute_min: float = 0.3):
    """
    Compute the tip-band half-width (µm) used by §5.4 from the data itself.

    Uses the larger of:
       (i) coverage_nodes × median inter-node spacing in the tip region
           (catches the clean-skeleton case where nodes are sparse), AND
      (ii) min_fraction_of_range × (path_length range)
           (catches the pixel-skeleton case where nodes cluster densely),
      and a fixed absolute_min as a final floor.

    Returns a positive float in µm.
    """
    import numpy as np
    if 'path_length_from_soma' not in df.columns:
        return absolute_min
    pl = df['path_length_from_soma']
    pmax = float(pl.max())
    pmin = float(pl.min())
    rng = max(pmax - pmin, 0.0)
    if pmax <= 0:
        return absolute_min
    tip_nodes = pl[pl > 0.8 * pmax]
    spacing_band = 0.0
    if tip_nodes.size >= 3:
        diffs = np.diff(np.sort(tip_nodes.unique()))
        diffs = diffs[diffs > 0]
        if diffs.size > 0:
            spacing_band = coverage_nodes * float(np.median(diffs))
    range_band = min_fraction_of_range * rng
    return max(absolute_min, spacing_band, range_band)


def adaptive_time_windows(df, n_windows: int = 3,
                           rate_threshold: float = 0.0):
    """
    Split the simulation into n_windows equal segments of the ACTIVE
    elongation phase (frames where elongation_rate > rate_threshold).
    Falls back to splitting the full time range if no active phase.
    Returns a list of n_windows lists of time_step values (ints).
    """
    import numpy as np
    if 'time_step' not in df.columns:
        return [[]]
    steps = sorted(df['time_step'].unique())
    if 'elongation_rate' in df.columns:
        elong = df.groupby('time_step')['elongation_rate'].mean()
        active = [int(s) for s in steps if elong.get(s, 0.0) > rate_threshold]
    else:
        active = []
    if len(active) < n_windows:
        active = [int(s) for s in steps]
    if len(active) == 0:
        return [[]]
    return [list(seg) for seg in np.array_split(np.array(active), n_windows)
            if len(seg) > 0]


def fit_lambda_T(xi_vals, T_vals):
    """
    Fit T(ξ) = T_inf + (T_tip − T_inf)·exp(ξ/λ_T) with adaptive bounds:
        T_inf  ∈ [0, mean(T)]
        T_tip  ∈ [mean(T), 2·max(T)]
        λ_T    ∈ [median inter-ξ spacing, 3·(max|ξ| − min|ξ|)]
    Returns (lam, lam_se, T_inf, T_tip, success_bool).
    """
    import numpy as np
    from scipy.optimize import curve_fit
    xi = np.asarray(xi_vals, dtype=float)
    T = np.asarray(T_vals, dtype=float)
    m = np.isfinite(xi) & np.isfinite(T)
    xi = xi[m]; T = T[m]
    if xi.size < 5:
        return (float('nan'),) * 4 + (False,)
    T_min = float(T.min()); T_max = float(T.max())
    T_mean = float(T.mean())
    xi_range = float(np.abs(xi).max() - np.abs(xi).min()) + 1e-6
    diffs = np.diff(np.sort(np.unique(xi)))
    diffs = diffs[diffs > 0]
    min_spacing = float(diffs.min()) + 1e-9 if diffs.size else 1e-3

    def model(xi_in, T_inf, T_tip, lam):
        return T_inf + (T_tip - T_inf) * np.exp(xi_in / lam)

    p0 = [max(T_min, 0.0), T_max, xi_range / 3]
    lo = [0.0, T_mean, min_spacing]
    hi = [max(T_mean, 1e-9), max(T_max * 2.0, T_min + 1e-9), xi_range * 3.0]
    # Clamp p0 into [lo, hi]
    p0 = [min(max(p0[i], lo[i]), hi[i]) for i in range(3)]
    try:
        popt, pcov = curve_fit(model, xi, T, p0=p0, bounds=(lo, hi),
                               maxfev=5000)
        perr = np.sqrt(np.diag(pcov))
        return (float(popt[2]), float(perr[2]),
                float(popt[0]), float(popt[1]), True)
    except Exception:
        return (float('nan'),) * 4 + (False,)


def select_velocity_column(df):
    """
    Prefer the real interface normal velocity (`norm_vel_int`) over the
    `elongation_rate` proxy. Returns (column_name, is_proxy, description).
    Raises ValueError if neither is present.
    """
    if 'norm_vel_int' in df.columns:
        v = df['norm_vel_int'].dropna()
        if v.size > 0 and float(v.std()) > 1e-9:
            return 'norm_vel_int', False, \
                   'local normal velocity from level-set (v_n)'
    if 'growth_vel_n' in df.columns:
        v = df['growth_vel_n'].dropna()
        if v.size > 0 and float(v.std()) > 1e-9:
            return 'growth_vel_n', False, \
                   'GrowthVel projected onto SWC tangent'
    if 'elongation_rate' in df.columns:
        return 'elongation_rate', True, \
               'dLmax/dt proxy (use only if norm_vel_int absent)'
    raise ValueError("No velocity column found in CSV.")


def available_fields(df, candidates, min_finite: int = 10,
                    min_std: float = 1e-12):
    """
    Return the subset of `candidates` that:
        - exists as a column
        - has at least `min_finite` finite values
        - has std > min_std (avoids trivially-constant fields)
    """
    out = []
    skipped = {}
    for f in candidates:
        if f not in df.columns:
            skipped[f] = 'missing column'
            continue
        v = df[f].dropna()
        if v.size < min_finite:
            skipped[f] = f'only {v.size} finite values'
            continue
        try:
            std = float(v.std())
        except Exception:
            std = 0.0
        if std < min_std:
            skipped[f] = f'std={std:.3g}<{min_std:g}'
            continue
        out.append(f)
    return out, skipped


# -----------------------------------------------------------------------------
# Plot footer helper
# -----------------------------------------------------------------------------

def add_exploratory_footer(fig, csv_path: str, extra: str = ''):
    txt = EXPLORATORY_FOOTER + (f'  ·  {extra}' if extra else '') + \
          f'\nsource CSV: {csv_path}'
    fig.text(0.02, 0.005, txt, ha='left', va='bottom',
             fontsize=6, color='#555555', style='italic')
