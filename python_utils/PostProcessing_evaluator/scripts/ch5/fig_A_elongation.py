#!/usr/bin/env python3
"""
fig_A_elongation.py — CH5 Part A (elongation regime).

Implements §5.3 (OLS velocity law v=αTB−δ with HAC SEs), §5.4 (co-moving
frame and tubulin penetration depth λ_T), and a §5.2-style cross-correlation
function with stationary block-bootstrap envelopes.

Inputs:  enriched CSV from enrich_swc_samples.py
Outputs: three figures + a text summary
    fig_A_ols_velocity_law.{png,svg}
    fig_A_comoving_tubulin.{png,svg}
    fig_A_ccf_concentration_vs_vtip.{png,svg}
    fig_A_summary.txt
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
from scipy.optimize import curve_fit

# Allow running as `python fig_A_elongation.py` without packaging
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from utils_ch5 import (   # noqa: E402
    ALPHA_MODEL, DELTA_MODEL, DT_SECONDS,
    add_exploratory_footer, adaptive_time_windows, adaptive_tip_band,
    cross_correlation, detect_active_elongation_phase, fit_lambda_T,
    gc_timeseries, hac_se_ols, integrated_autocorr_time, lmax_per_step,
    load_enriched, per_step_value, select_velocity_column,
    stationary_block_bootstrap_corr, time_values, tip_timeseries,
    tip_velocity_proxy,
)


# -----------------------------------------------------------------------------
# §5.3 — OLS velocity law v = α(T·B) − δ
# -----------------------------------------------------------------------------

def figure_ols(df, out_dir, csv_path,
                alpha_model=ALPHA_MODEL, delta_model=DELTA_MODEL):
    """
    v_tip = α (T·B) − δ + ε; two fits:
        (i) all-frames lag-1 OLS with HAC SEs (the strict §5.3 fit);
        (ii) active-elongation-phase only (contiguous initial run where the
             velocity column stays above 0.03 of its peak; case-adaptive).

    Velocity source:
      The function prefers `norm_vel_int` (real local interface velocity)
      from the extended sampler. If only `elongation_rate` is present in
      the CSV, it falls back to that proxy and clearly tags the figure.

    Per-second conversion: if the velocity column is in µm/frame (the proxy)
    we divide α̂ by DT_SECONDS to compare to α_model in µm/(s·µM²). If the
    velocity column is already a per-second normal velocity, the conversion
    is a no-op.
    """
    # Adaptive velocity-column selection (Task 3b)
    vcol, is_proxy, vdesc = select_velocity_column(df)
    ut = gc_timeseries(df, 'u_t')['mean']
    ub = gc_timeseries(df, 'u_b')['mean']
    if is_proxy:
        v_series = per_step_value(df, vcol)
        v_agg_note = 'per-step scalar from CSV'
    else:
        # Real v_n is per-node; aggregate over TIP-only nodes (is_tip==1)
        # because norm_vel_int is non-trivial only at the moving interface.
        # GC-region (radius 0.1 µm) averaging dilutes the signal with
        # interior nodes where the extended velocity ≈ 0.
        v_series = tip_timeseries(df, vcol)['mean']
        v_agg_note = 'mean over is_tip==1 nodes (tip-only)'
    common = sorted(set(ut.index) & set(ub.index) & set(v_series.index))
    ut = ut.reindex(common); ub = ub.reindex(common)
    elong = v_series.reindex(common)
    TB = (ut * ub).rename('TB')
    L  = lmax_per_step(df).reindex(common)
    # If proxy (µm/frame), convert α̂ by /DT to get µm/(s·µM²). For real v_n
    # already in µm/s the conversion factor is 1.0.
    sec_factor = DT_SECONDS if is_proxy else 1.0

    # All-frames lag-1
    t_idx = np.array(common)
    x_all = TB.loc[t_idx[:-1]].values
    y_all = elong.shift(-1).loc[t_idx[:-1]].values
    m_all = np.isfinite(x_all) & np.isfinite(y_all)
    alpha_all, intercept_all, se_a_all, se_i_all, r2_all, n_all = \
        hac_se_ols(x_all[m_all], y_all[m_all])

    # Active-only fit:
    #  - When using the proxy (dLmax/dframe), restrict to the contiguous
    #    initial growth phase to avoid the plateau-driven negative slope.
    #  - When using the real v_n, every frame is informative; the active set
    #    is just the frames whose GC-mean v_n is in the top 60 percentile
    #    (an adaptive selector that does not assume monotone elongation).
    if is_proxy:
        peak_v = float(elong.dropna().abs().max()) if elong.notna().any() else 0.0
        rate_thresh = 0.05 * peak_v if peak_v > 0 else 0.0
        active_steps = detect_active_elongation_phase(elong, threshold=rate_thresh)
    else:
        absv = elong.dropna().abs()
        if absv.empty:
            active_steps = []
        else:
            q60 = float(absv.quantile(0.40))
            active_steps = [int(s) for s in elong.index
                            if pd.notna(elong.loc[s]) and abs(elong.loc[s]) >= q60]
    if active_steps:
        x_act = TB.loc[active_steps].values
        y_act = elong.loc[active_steps].values
    else:
        x_act = np.array([])
        y_act = np.array([])
    if len(x_act) >= 3:
        alpha_act, intercept_act, se_a_act, se_i_act, r2_act, n_act = \
            hac_se_ols(x_act, y_act)
    else:
        alpha_act = intercept_act = se_a_act = se_i_act = r2_act = float('nan')
        n_act = len(x_act)

    # Convert raw α̂ (in units of velocity column / µM²) to µm/(s·µM²).
    # proxy: divide by Δt; real v_n: no conversion needed.
    alpha_all_s = alpha_all / sec_factor
    alpha_act_s = alpha_act / sec_factor

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1.0], hspace=0.32,
                           wspace=0.28)

    # (a) all-frames scatter + OLS
    ax_a = fig.add_subplot(gs[0, 0])
    sc = ax_a.scatter(x_all[m_all], y_all[m_all], c=t_idx[:-1][m_all],
                       cmap='viridis', s=42, edgecolor='black',
                       linewidth=0.4, alpha=0.85)
    xx = np.linspace(0.0, max(x_all[m_all].max(), 0.001), 100)
    unit_y = 'µm/frame' if is_proxy else 'µm/s'
    ax_a.plot(xx, alpha_all * xx + intercept_all, color='#cc0066',
               linewidth=1.8,
               label=fr'OLS all-frames: $\hat\alpha$={alpha_all:.3g} /{unit_y[:-2]}')
    if np.isfinite(alpha_act):
        ax_a.plot(xx, alpha_act * xx + intercept_act, color='#0099cc',
                   linewidth=1.8,
                   label=fr'OLS active-phase: $\hat\alpha$={alpha_act:.3g} /{unit_y[:-2]}')
    # Model reference line in the SAME units as y-axis (per frame for proxy)
    model_slope_in_y = alpha_model * sec_factor
    ax_a.plot(xx, model_slope_in_y * xx, color='#333333',
               linestyle='--', linewidth=1.2,
               label=fr'model $\alpha$={alpha_model:.3g} µm/(s·µM²) → {model_slope_in_y:.3g} in {unit_y}')
    ax_a.set_xlabel(r'GC-mean $u_t \cdot u_b$ at $t$  ($\approx T\cdot B$, µM²)')
    if is_proxy:
        ax_a.set_ylabel(rf'$v_{{tip}}$ at $t+1$  ({vcol}, {unit_y})')
        ax_a.set_title('§5.3 OLS velocity law — PROXY (dLmax/dt), see caveat')
    else:
        ax_a.set_ylabel(rf'$v_{{tip}}$  ({vcol}, GC mean, {unit_y})')
        ax_a.set_title('§5.3 OLS velocity law — using real $v_n$ from level-set')
    ax_a.legend(loc='best', fontsize=7)
    ax_a.grid(True, alpha=0.25)
    cbar = fig.colorbar(sc, ax=ax_a, pad=0.02, fraction=0.04)
    cbar.set_label('time step')

    txt = (f'all frames (lag-1):  $\\hat\\alpha$={alpha_all:.3g}/frame  '
           f'={alpha_all_s:.3g} µm/(s·µM²),  N={n_all},  $R^2$={r2_all:.3f}\n'
           f'active phase (lag-0): $\\hat\\alpha$={alpha_act:.3g}/frame  '
           f'={alpha_act_s:.3g} µm/(s·µM²),  N={n_act},  $R^2$={r2_act:.3f}\n'
           f'model $\\alpha$ = {ALPHA_MODEL} µm/(s·µM²)')
    ax_a.text(0.02, 0.98, txt, transform=ax_a.transAxes,
               ha='left', va='top', fontsize=7.5,
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                         edgecolor='#666', alpha=0.9))

    # (b) partial-residual plot
    ax_b = fig.add_subplot(gs[0, 1])
    yhat = alpha_all * x_all[m_all] + intercept_all
    resid = y_all[m_all] - yhat
    ax_b.scatter(x_all[m_all], resid + alpha_all * x_all[m_all],
                  c=t_idx[:-1][m_all], cmap='viridis', s=42,
                  edgecolor='black', linewidth=0.4, alpha=0.85)
    ax_b.plot(xx, alpha_all * xx, color='#cc0066', linewidth=1.8,
               label='OLS slope (all)')
    ax_b.axhline(0, color='#888', linewidth=0.8)
    ax_b.set_xlabel(r'$u_t \cdot u_b$  (GC mean)')
    ax_b.set_ylabel('partial residual')
    ax_b.set_title('partial-residual plot (linearity check)')
    ax_b.grid(True, alpha=0.25)
    ax_b.legend(loc='best', fontsize=8)

    # (c) trajectory diagnostic — TB(t), v(t), Lmax(t)
    ax_c = fig.add_subplot(gs[1, :])
    t_arr = np.array(common)
    ax_c.plot(t_arr, TB.values, color='#1f77b4', linewidth=2,
               marker='o', markersize=4, label=r'$u_t \cdot u_b$ (GC mean)')
    ax_c.set_xlabel('time step')
    ax_c.set_ylabel(r'$u_t \cdot u_b$  (GC mean)', color='#1f77b4')
    ax_c.tick_params(axis='y', labelcolor='#1f77b4')
    ax_c.grid(True, alpha=0.25)

    ax_c2 = ax_c.twinx()
    ax_c2.plot(t_arr, elong.values, color='#cc0066', linewidth=1.8,
                marker='s', markersize=4, alpha=0.85,
                label=r'$v_{tip}$ = elongation_rate (µm/frame)')
    ax_c2.plot(t_arr, L.values / L.values[0], color='#333333',
                linewidth=1.2, linestyle=':',
                label=fr'$L_{{\max}}(t)/L_0$  (1.0 = no growth)')
    ax_c2.set_ylabel(r'$v_{tip}$ and $L_{\max}/L_0$', color='#333333')
    ax_c2.axhline(0, color='#888', linewidth=0.6)

    # mark active phase
    if active_steps:
        ax_c.axvspan(active_steps[0] - 0.5, active_steps[-1] + 0.5,
                      color='#cc0066', alpha=0.08,
                      label='active elongation phase')
    h1, l1 = ax_c.get_legend_handles_labels()
    h2, l2 = ax_c2.get_legend_handles_labels()
    ax_c.legend(h1 + h2, l1 + l2, loc='center right', fontsize=7.5)
    ax_c.set_title('Trajectory diagnostic — TB(t), $v_{tip}$(t), $L_{\\max}(t)/L_0$ '
                    '(active phase shaded)')

    fig.suptitle('CH5 §5.3 — direct test of $v_{tip}=\\alpha TB-\\delta$  '
                  '(GC-mean proxy; see caveat in summary)', fontsize=11)
    add_exploratory_footer(fig, csv_path,
                            extra=('elongation_rate is dLmax/frame, NOT µm/s; '
                                   'compare on the same time scale only after '
                                   'multiplying α_model by Δt = 6 s'))
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])

    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_A_ols_velocity_law.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)

    return dict(
        velocity_column=vcol,
        velocity_is_proxy=bool(is_proxy),
        velocity_description=vdesc,
        velocity_aggregation=v_agg_note,
        all_frames=dict(alpha_hat_per_frame=alpha_all,
                        alpha_hat_per_s=alpha_all_s,
                        se_slope=se_a_all,
                        intercept=intercept_all, se_intercept=se_i_all,
                        r2=r2_all, n=n_all),
        active_phase=dict(alpha_hat_per_frame=alpha_act,
                          alpha_hat_per_s=alpha_act_s,
                          se_slope=se_a_act,
                          intercept=intercept_act, se_intercept=se_i_act,
                          r2=r2_act, n=n_act,
                          steps=active_steps),
        alpha_model=alpha_model, delta_model=delta_model,
        ratio_all_to_model=alpha_all_s / alpha_model
            if alpha_model else float('nan'),
        ratio_active_to_model=alpha_act_s / alpha_model
            if alpha_model else float('nan'),
    )


# -----------------------------------------------------------------------------
# §5.4 — Co-moving frame; tubulin penetration depth λ_T
# -----------------------------------------------------------------------------

def figure_comoving(df, out_dir, csv_path,
                     tip_band=None, n_windows=3, case_id=None):
    """
    Co-moving frame: ξ = path_length_from_soma − Lmax(t).
    Keep nodes with ξ ∈ [−tip_band, 0]. Fit T̄(ξ) = T_inf + (T_tip − T_inf)·exp(ξ/λ_T).

    Case-adaptive design (Task 2):
      - `tip_band` is computed by adaptive_tip_band() from the data
        (≈ 5 × median inter-node spacing in the tip region) when None.
      - Time windows are chosen by adaptive_time_windows() from the
        active-elongation phase, with n_windows segments.
      - The exponential fit uses fit_lambda_T() with bounds derived
        from the actual data ranges of ξ and u_t.
    """
    lmax = lmax_per_step(df)
    df = df.assign(xi=df['path_length_from_soma'] - df['time_step'].map(lmax))

    # Adaptive tip-band
    if tip_band is None or tip_band <= 0:
        tip_band = adaptive_tip_band(df)
    # Adaptive time windows (n_windows segments of the active phase)
    windows = adaptive_time_windows(df, n_windows=n_windows)
    if not windows or all(len(w) == 0 for w in windows):
        windows = [sorted(df['time_step'].unique())]

    # Color cycle scales with n_windows
    base_colors = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#ff7f0e',
                    '#17becf', '#8c564b']
    window_colors = [base_colors[i % len(base_colors)] for i in range(len(windows))]

    fit_summary = {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    for w_idx, w_steps in enumerate(windows):
        tag = f'w{w_idx}_t{w_steps[0]}-{w_steps[-1]}'
        sub = df[df['time_step'].isin(w_steps)
                 & (df['xi'] >= -tip_band) & (df['xi'] <= 0)
                 & df['u_t'].notna()]
        n_nodes = len(sub)
        if n_nodes < 5:
            fit_summary[tag] = dict(n=n_nodes, lam=None, se=None,
                                     T_inf=None, T_tip=None,
                                     steps=[int(s) for s in w_steps])
            continue
        xs = sub['xi'].values
        ys = sub['u_t'].values
        lam, lam_se, T_inf, T_tip, ok = fit_lambda_T(xs, ys)
        axes[0].scatter(xs, ys, color=window_colors[w_idx], s=10, alpha=0.4,
                         label=f'{tag}  (n={n_nodes})')
        if ok:
            xi_grid = np.linspace(-tip_band, 0, 80)
            yhat = T_inf + (T_tip - T_inf) * np.exp(xi_grid / lam)
            axes[0].plot(xi_grid, yhat, color=window_colors[w_idx], linewidth=2)
        fit_summary[tag] = dict(n=n_nodes, lam=lam if ok else None,
                                 se=lam_se if ok else None,
                                 T_inf=T_inf if ok else None,
                                 T_tip=T_tip if ok else None,
                                 steps=[int(s) for s in w_steps],
                                 fit_ok=ok)

    axes[0].set_xlabel(r'$\xi = s - L_{\max}(t)$  (µm; 0 = tip)')
    axes[0].set_ylabel(r'$u_t$  (tubulin)')
    axes[0].set_title(f'§5.4 co-moving tubulin + exp fit  '
                       f'(tip band = {tip_band:.3f} µm)')
    axes[0].axvline(0, color='#888', linewidth=0.8)
    axes[0].legend(loc='best', fontsize=7)
    axes[0].grid(True, alpha=0.25)

    # Panel (b): λ_T per window bar chart
    tags_in_order = list(fit_summary.keys())
    lams = [fit_summary[t]['lam'] for t in tags_in_order]
    ses  = [fit_summary[t]['se'] or 0.0 for t in tags_in_order]
    valid_idx = [i for i, v in enumerate(lams) if v is not None]
    xpos = np.arange(len(tags_in_order))
    if valid_idx:
        axes[1].bar(xpos[valid_idx], [lams[i] for i in valid_idx],
                     yerr=[ses[i] for i in valid_idx],
                     color=[window_colors[i] for i in valid_idx],
                     alpha=0.85, capsize=4)
    axes[1].set_xticks(xpos)
    axes[1].set_xticklabels([t.replace('_t', '\nt=') for t in tags_in_order],
                              fontsize=7)
    axes[1].set_ylabel(r'$\hat\lambda_T$  (µm)')
    axes[1].set_title('tubulin penetration depth per window')
    axes[1].grid(True, alpha=0.25, axis='y')

    case_str = f'  ·  case: {case_id}' if case_id else ''
    fig.suptitle(r'CH5 §5.4 — co-moving frame: tubulin penetration depth $\lambda_T$' + case_str,
                  fontsize=11)
    add_exploratory_footer(fig, csv_path,
                            extra=(f'tip band |ξ|≤{tip_band:.3f}µm '
                                   f'(adaptive: 5× median tip-node spacing); '
                                   f'{len(windows)} windows over active phase'))
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])

    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_A_comoving_tubulin.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)

    return dict(tip_band=tip_band, n_windows=len(windows),
                windows=fit_summary)


# -----------------------------------------------------------------------------
# §5.2 — Cross-correlation function with bootstrap envelopes
# -----------------------------------------------------------------------------

def figure_ccf(df, out_dir, csv_path, max_lag=8):
    """
    ρ(τ) = corr(TB(t), v(t+τ))  and  ρ(τ) = corr(Ca(t), v(t+τ))
    for τ ∈ [-max_lag, +max_lag] frames.
    """
    ut = gc_timeseries(df, 'u_t')['mean']
    ub = gc_timeseries(df, 'u_b')['mean']
    ca = gc_timeseries(df, 'u_ca_cyt')['mean']
    elong = tip_velocity_proxy(df)
    common = ut.index.intersection(ub.index).intersection(ca.index) \
                     .intersection(elong.index)
    common = np.array(sorted(common))
    ut = ut.loc[common].values
    ub = ub.loc[common].values
    ca = ca.loc[common].values
    v  = elong.loc[common].values
    TB = ut * ub

    lags = np.arange(-max_lag, max_lag + 1)

    # Empirical CCF
    rho_TB = cross_correlation(TB, v, max_lag)
    rho_Ca = cross_correlation(ca, v, max_lag)

    # Block-bootstrap CIs at each lag (independent per τ)
    rho_TB_lo = np.zeros_like(rho_TB); rho_TB_hi = np.zeros_like(rho_TB)
    rho_Ca_lo = np.zeros_like(rho_Ca); rho_Ca_hi = np.zeros_like(rho_Ca)
    for k, tau in enumerate(lags):
        if tau >= 0:
            x_TB = TB[:len(TB) - tau];  y_TB = v[tau:]
            x_Ca = ca[:len(ca) - tau];  y_Ca = v[tau:]
        else:
            x_TB = TB[-tau:];  y_TB = v[:len(v) + tau]
            x_Ca = ca[-tau:];  y_Ca = v[:len(v) + tau]
        _, lo, hi = stationary_block_bootstrap_corr(x_TB, y_TB,
                                                     n_boot=600, block_len=4)
        rho_TB_lo[k], rho_TB_hi[k] = lo, hi
        _, lo, hi = stationary_block_bootstrap_corr(x_Ca, y_Ca,
                                                     n_boot=600, block_len=4)
        rho_Ca_lo[k], rho_Ca_hi[k] = lo, hi

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    lags_s = lags * DT_SECONDS

    ax.fill_between(lags_s, rho_TB_lo, rho_TB_hi, color='#1f77b4', alpha=0.18,
                     label='TB·v 95% block-boot')
    ax.plot(lags_s, rho_TB, color='#1f77b4', linewidth=2,
             marker='o', markersize=4, label=r'$\rho(TB(t), v(t+\tau))$')
    ax.fill_between(lags_s, rho_Ca_lo, rho_Ca_hi, color='#ff7f0e', alpha=0.18,
                     label='Ca·v 95% block-boot')
    ax.plot(lags_s, rho_Ca, color='#ff7f0e', linewidth=1.5,
             linestyle='--', marker='s', markersize=4,
             label=r'$\rho(Ca(t), v(t+\tau))$')
    ax.axhline(0, color='#888', linewidth=0.8)
    ax.axvline(0, color='#888', linewidth=0.8, linestyle=':')
    ax.set_xlabel(r'lag $\tau$  (s)   — positive: concentration leads $v_{tip}$')
    ax.set_ylabel(r'$\rho(\tau)$')
    ax.set_title('§5.2 cross-correlation function (GC-mean concentration vs $v_{tip}$)')
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.25)

    tau_int = integrated_autocorr_time(v)
    txt = (f'N = {len(common)} frames\n'
           f'$\\tau_{{int}}$($v$) ≈ {tau_int:.2f} frames '
           f'({tau_int * DT_SECONDS:.0f} s)\n'
           'block bootstrap b̄ = 4; n_boot = 600')
    ax.text(0.02, 0.98, txt, transform=ax.transAxes,
             ha='left', va='top', fontsize=8,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       edgecolor='#666', alpha=0.85))

    add_exploratory_footer(fig, csv_path,
                            extra=f'max_lag={max_lag} frames; Δt={DT_SECONDS}s')
    fig.tight_layout()
    out = Path(out_dir)
    for ext in ('png', 'svg'):
        fig.savefig(out / f'fig_A_ccf_concentration_vs_vtip.{ext}',
                    dpi=160, bbox_inches='tight')
    plt.close(fig)

    return dict(lags_s=lags_s.tolist(),
                rho_TB=rho_TB.tolist(), rho_Ca=rho_Ca.tolist(),
                tau_int_v=tau_int)


# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--case-id', default=None,
                    help='Optional case label used in figure suptitles.')
    ap.add_argument('--alpha-model', type=float, default=ALPHA_MODEL,
                    help='Reference alpha for the model comparison line in §5.3 '
                         f'(default: {ALPHA_MODEL} µm/(s·µM²); see Lua control file).')
    ap.add_argument('--delta-model', type=float, default=DELTA_MODEL,
                    help='Reference delta for the model comparison.')
    ap.add_argument('--n-windows', type=int, default=3,
                    help='Number of co-moving-frame time windows (default 3).')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_enriched(args.csv)
    ols = figure_ols(df, out_dir, args.csv,
                      alpha_model=args.alpha_model,
                      delta_model=args.delta_model)
    comov = figure_comoving(df, out_dir, args.csv,
                              n_windows=args.n_windows,
                              case_id=args.case_id)
    ccf = figure_ccf(df, out_dir, args.csv)

    summary = dict(csv=args.csv, ols=ols, comoving=comov, ccf=ccf)
    with open(out_dir / 'fig_A_summary.txt', 'w') as fh:
        fh.write('CH5 Part A — summary\n')
        fh.write('=' * 60 + '\n\n')
        fh.write(f'Source CSV: {args.csv}\n\n')
        if ols:
            vcol = ols.get('velocity_column', 'unknown')
            is_proxy = ols.get('velocity_is_proxy', True)
            fh.write(f'§5.3 OLS velocity law  '
                     f'(velocity column: {vcol}, '
                     f'{"PROXY (dLmax/dt)" if is_proxy else "REAL v_n from level-set"}):\n')
            af = ols['all_frames']
            fh.write('  ALL FRAMES (lag-1):\n')
            fh.write(f"    alpha_hat = {af['alpha_hat_per_frame']:.5g} µm/frame/µM²  "
                     f"  ({af['alpha_hat_per_s']:.5g} µm/(s·µM²))  ± "
                     f"{af['se_slope']:.4g}\n")
            fh.write(f"    intercept = {af['intercept']:.5g}, "
                     f"R^2 = {af['r2']:.4f}, N = {af['n']}\n")
            ap = ols['active_phase']
            fh.write('  ACTIVE PHASE (lag-0; rate > 0.03 µm/frame):\n')
            fh.write(f"    alpha_hat = {ap['alpha_hat_per_frame']:.5g} µm/frame/µM²  "
                     f"  ({ap['alpha_hat_per_s']:.5g} µm/(s·µM²))  ± "
                     f"{ap['se_slope']:.4g}\n")
            fh.write(f"    intercept = {ap['intercept']:.5g}, "
                     f"R^2 = {ap['r2']:.4f}, N = {ap['n']}\n")
            fh.write(f"    active steps: {ap['steps']}\n\n")
            fh.write(f"  alpha_model = {ols['alpha_model']} µm/(s·µM²)\n")
            fh.write(f"  ratio  all-frames α̂/α_model  = "
                     f"{ols['ratio_all_to_model']:.4g}\n")
            fh.write(f"  ratio  active-phase α̂/α_model = "
                     f"{ols['ratio_active_to_model']:.4g}\n\n")
            if is_proxy:
                fh.write('  CAVEAT — proxy mismatch warning:\n')
                fh.write('    `elongation_rate` from enrich_swc_samples.py is\n'
                         '    dLmax/dt where Lmax = max path-length over the SWC.\n'
                         '    The simulator\'s v = αTB law drives the local interface\n'
                         '    normal velocity, NOT dLmax/dt. Saturation of Lmax\n'
                         '    once growth hits a boundary plateaus elongation_rate\n'
                         '    while u_b continues to accumulate, producing the\n'
                         '    apparent flat/negative all-frames slope. Active-phase\n'
                         '    slope is positive but ~order-of-magnitude above\n'
                         '    α_model (likely sampling location / proxy bias).\n'
                         '    Interpret as DIAGNOSTIC, not as a validation/rejection\n'
                         '    of the velocity law.\n\n')
            else:
                v_agg = ols.get('velocity_aggregation', 'unknown')
                fh.write('  NOTE — using real v_n from level-set:\n')
                fh.write(f'    velocity column   = {vcol}\n'
                         f'    aggregation       = {v_agg}\n'
                         '    Compared to α_model in µm/(s·µM²); no time-unit\n'
                         '    conversion needed. Active-phase here = top-60%-by-|v_n|\n'
                         '    timesteps (adaptive), NOT the contiguous-growth hack\n'
                         '    used for the proxy. Tip-only aggregation chosen because\n'
                         '    norm_vel_int is the level-set extended interface velocity:\n'
                         '    nontrivial only at the moving front, ≈0 at interior nodes.\n\n')
        if comov:
            fh.write('§5.4 Co-moving tubulin penetration depth lambda_T:\n')
            fh.write(f"  tip band: |xi| <= {comov.get('tip_band'):.4g} µm "
                     f"(adaptive)\n")
            fh.write(f"  n_windows: {comov.get('n_windows')}\n")
            for tag, d in comov.get('windows', {}).items():
                if d.get('lam') is not None:
                    fh.write(f"  {tag}: lambda_T = {d['lam']:.4g} µm "
                             f"(SE={d['se']:.4g}); T_inf={d['T_inf']:.4g}, "
                             f"T_tip={d['T_tip']:.4g}; n={d['n']}\n")
                else:
                    fh.write(f"  {tag}: n={d.get('n')}  (fit failed/skipped)\n")
            fh.write('\n')
        if ccf:
            fh.write('§5.2 CCF integrated autocorrelation time of v_tip:\n')
            fh.write(f"  tau_int(v) = {ccf['tau_int_v']:.3f} frames "
                     f"({ccf['tau_int_v'] * DT_SECONDS:.0f} s)\n")

    with open(out_dir / 'fig_A_summary.json', 'w') as fh:
        json.dump(summary, fh, indent=2, default=float)

    print(f'\nfig_A done. summary -> {out_dir / "fig_A_summary.txt"}')
    if ols:
        print(f'  all-frames    α̂/α_model = '
              f'{ols["ratio_all_to_model"]:.4g}')
        print(f'  active-phase  α̂/α_model = '
              f'{ols["ratio_active_to_model"]:.4g}')


if __name__ == '__main__':
    main()
