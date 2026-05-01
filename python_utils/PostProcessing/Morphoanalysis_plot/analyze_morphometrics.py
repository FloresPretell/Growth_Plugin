#!/usr/bin/env python3
"""
Morphology analysis — dynamically adapts to all metric columns in morphometrics.csv.

Figures
  fig1   final-state D×V heatmaps       (key metrics, 2 rows × 4 cols)
  fig2   time-course by D               (key metrics)
  fig3   time-course by V               (key metrics)
  fig4   Spearman correlations with D/V (all metrics, sorted)
  fig5   pairwise metric correlation matrix
  fig6   PCA scores + loadings
  fig7   D×V interaction lines          (4 metrics)
  fig8   marginal effects D and V       (4 metrics)
  fig9   simulation length per case
  fig10  information gain scatter       (D/V sensitivity vs size redundancy)

CSVs
  summary_final_state.csv
  correlation_metrics_DV.csv
  pca_variance.csv  /  pca_loadings.csv
  info_gain_table.csv
"""

import argparse
import math
import re
import warnings
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ── style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        9,
    "axes.linewidth":   0.8,
    "axes.spines.top":  False,
    "axes.spines.right": False,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size":  3,
    "ytick.major.size":  3,
    "pdf.fonttype":     42,
    "svg.fonttype":     "none",
    "figure.dpi":       150,
})

# ── constants ─────────────────────────────────────────────────────────────────
IDENTIFIER_COLS = {"case", "timestep", "swc_file", "D", "V"}

METRIC_LABELS = {
    "n_nodes":                       "Nodes",
    "n_tips":                        "Terminal tips",
    "n_branchpoints":                "Branch points",
    "n_sections":                    "Sections",
    "total_length_px":               "Total length (px)",
    "mean_section_length_px":        "Mean section len (px)",
    "max_path_length_px":            "Max path len (px)",
    "mean_path_length_px":           "Mean path len (px)",
    "median_path_length_px":         "Median path len (px)",
    "terminal_section_length_mean_px": "Terminal section len (px)",
    "width_px":                      "Width (px)",
    "height_px":                     "Height (px)",
    "aspect_ratio":                  "Aspect ratio",
    "convex_hull_area_px2":          "Convex hull area (px²)",
    "occupancy_density":             "Occupancy density",
    "mean_tortuosity":               "Mean tortuosity",
    "max_tortuosity":                "Max tortuosity",
    "path_length_cv":                "Path length CV",
    "section_length_cv":             "Section length CV",
    "tip_branch_ratio":              "Tip/branch ratio",
    "left_right_asymmetry_length":   "LR asymmetry (length)",
    "left_right_asymmetry_tips":     "LR asymmetry (tips)",
    "subtree_asymmetry_mean":        "Subtree asymmetry",
    "mean_branch_angle_deg":         "Mean branch angle (°)",
    "std_branch_angle_deg":          "Branch angle SD (°)",
}

def mlabel(m):
    return METRIC_LABELS.get(m, m)

# One representative per conceptual group; order determines priority in figures
KEY_METRICS_PRIORITY = [
    "total_length_px",           # size
    "n_tips",                    # count
    "max_path_length_px",        # reach
    "mean_section_length_px",    # section detail
    "aspect_ratio",              # shape
    "occupancy_density",         # packing density
    "mean_tortuosity",           # path geometry
    "subtree_asymmetry_mean",    # topology
    "mean_branch_angle_deg",     # branching geometry
    "section_length_cv",         # variability
]

# Fixed 4-metric set for interaction and marginal figures
INTERACTION_METRICS = [
    "total_length_px",
    "n_tips",
    "mean_section_length_px",
    "occupancy_density",
]


# ── data helpers ──────────────────────────────────────────────────────────────
def load_data(path="morphometrics.csv"):
    df = pd.read_csv(path)

    def parse_case(s):
        m = re.match(r"D(?P<D>[0-9.]+)_V(?P<V>[0-9]+(?:p[0-9]+)?|[0-9.]+)$", s)
        return float(m.group("D")), float(m.group("V").replace("p", "."))

    df[["D", "V"]] = pd.DataFrame(df["case"].apply(parse_case).tolist(), index=df.index)
    df["timestep"] = df["timestep"].astype(int)
    return df


def detect_metrics(df):
    return [c for c in df.columns if c not in IDENTIFIER_COLS]


def select_key_metrics(available, n=8):
    km = [m for m in KEY_METRICS_PRIORITY if m in available]
    if len(km) < n:
        extras = [m for m in available if m not in km]
        km = km + extras[: n - len(km)]
    return km[:n]


def final_state(df):
    return df.loc[df.groupby("case")["timestep"].idxmax()].copy()


def save_fig(fig, stem, outdir="."):
    fig.savefig(f"{outdir}/{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{outdir}/{stem}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {stem}")


def save_placeholder_fig(stem, title, message, outdir="."):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.axis("off")
    ax.text(0.5, 0.62, title, ha="center", va="center", fontsize=11, weight="bold")
    ax.text(0.5, 0.40, message, ha="center", va="center", fontsize=9, wrap=True)
    save_fig(fig, stem, outdir)


def scalar_norm(values):
    vmin = min(values)
    vmax = max(values)
    if vmin == vmax:
        pad = 0.5 if vmin == 0 else max(abs(vmin) * 0.05, 0.5)
        return mcolors.Normalize(vmin=vmin - pad, vmax=vmax + pad)
    return mcolors.Normalize(vmin=vmin, vmax=vmax)


def palette_from_values(values, cmap):
    if len(values) == 1:
        return [cmap(0.5)]
    denom = max(len(values) - 1, 1)
    return [cmap(i / denom) for i in range(len(values))]


# ── fig1: final-state heatmaps ────────────────────────────────────────────────
def plot_final_heatmaps(final, metrics, outdir="."):
    D_vals = sorted(final["D"].unique())
    V_vals = sorted(final["V"].unique())
    n = len(metrics)
    ncols = min(n, 4)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(3.9 * ncols, 3.5 * nrows),
                             squeeze=False)
    for idx, metric in enumerate(metrics):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        grid = (
            final.pivot_table(index="V", columns="D", values=metric, aggfunc="mean")
            .reindex(index=V_vals[::-1], columns=D_vals)
        )
        im = ax.imshow(grid.values, aspect="auto", cmap="viridis",
                       extent=[-0.5, len(D_vals) - 0.5, -0.5, len(V_vals) - 0.5])
        ax.set_xticks(range(len(D_vals)))
        ax.set_xticklabels([f"{v:g}" for v in D_vals], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(V_vals)))
        ax.set_yticklabels([f"{v:g}" for v in V_vals[::-1]], fontsize=7)
        ax.set_xlabel("D", fontsize=8)
        if col == 0:
            ax.set_ylabel("V", fontsize=8)
        ax.set_title(mlabel(metric), fontsize=8, pad=4)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=6)

    for idx in range(n, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.suptitle("Final-state morphology — D × V parameter space", fontsize=10, y=1.01)
    fig.tight_layout()
    save_fig(fig, "fig1_final_heatmaps", outdir)


# ── fig2/3: time-course ───────────────────────────────────────────────────────
def plot_timecourse(df, metrics, by, outdir="."):
    other = "V" if by == "D" else "D"
    vals = sorted(df[by].unique())
    cmap = plt.cm.plasma if by == "D" else plt.cm.coolwarm
    colors = palette_from_values(vals, cmap)
    n = len(metrics)

    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 3.4), constrained_layout=True,
                             squeeze=False)
    axes = axes[0]

    for ax, metric in zip(axes, metrics):
        for val, color in zip(vals, colors):
            ts_mean = df[df[by] == val].groupby("timestep")[metric].mean()
            ax.plot(ts_mean.index, ts_mean.values, color=color, linewidth=1.0, alpha=0.9)
        ax.set_xlabel("Timestep", fontsize=9)
        ax.set_title(mlabel(metric), fontsize=9, pad=4)

    axes[0].set_ylabel(f"Mean over {other}", fontsize=9)

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=scalar_norm(vals))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=list(axes), fraction=0.015, pad=0.02, label=by)
    cb.ax.tick_params(labelsize=7)

    num = "2" if by == "D" else "3"
    fig.suptitle(f"Time-course by {by} (averaged over {other})", fontsize=10)
    save_fig(fig, f"fig{num}_timecourse_by_{by}", outdir)


# ── fig4: Spearman correlations with D and V ──────────────────────────────────
def plot_correlations_DV(final, metrics, outdir="."):
    records = []
    for m in metrics:
        col = final[m].dropna()
        if len(col) < 10:
            continue
        idx = col.index
        r_D, _ = stats.spearmanr(final.loc[idx, "D"], col)
        r_V, _ = stats.spearmanr(final.loc[idx, "V"], col)
        records.append({"metric": m, "rho_D": r_D, "rho_V": r_V})

    if not records:
        pd.DataFrame(columns=["metric", "rho_D", "rho_V"]).to_csv(
            f"{outdir}/correlation_metrics_DV.csv", index=False
        )
        save_placeholder_fig(
            "fig4_correlations_DV",
            "Spearman correlations with D and V",
            "Not enough final-state cases to compute stable D/V correlations.",
            outdir,
        )
        return pd.DataFrame(columns=["rho_D", "rho_V"])

    corr_df = (
        pd.DataFrame(records)
        .set_index("metric")
        .assign(max_abs=lambda x: x[["rho_D", "rho_V"]].abs().max(axis=1))
        .sort_values("max_abs", ascending=True)      # ascending → most sensitive at top of barh
    )
    corr_df.drop(columns="max_abs").to_csv(f"{outdir}/correlation_metrics_DV.csv")

    n = len(corr_df)
    fig, axes = plt.subplots(1, 2, figsize=(12, max(4.0, n * 0.42 + 1.2)),
                             constrained_layout=True)
    y = np.arange(n)
    labels = [mlabel(m) for m in corr_df.index]

    for ax, col_name, base_color, title in [
        (axes[0], "rho_D", "#2c7bb6", "Spearman ρ with D"),
        (axes[1], "rho_V", "#d7191c", "Spearman ρ with V"),
    ]:
        vals = corr_df[col_name].values
        bar_colors = [base_color if v >= 0 else "#aaaaaa" for v in vals]
        ax.barh(y, vals, color=bar_colors, alpha=0.85, height=0.65)
        ax.axvline(0, color="k", linewidth=0.6)
        ax.axvline( 0.3, color="gray", linewidth=0.5, linestyle="--", alpha=0.6)
        ax.axvline(-0.3, color="gray", linewidth=0.5, linestyle="--", alpha=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(-1.05, 1.05)
        ax.set_xlabel(title, fontsize=9)
        ax.set_title(title, fontsize=9)

    fig.suptitle("Spearman correlations with D and V — final state (sorted by |ρ|)",
                 fontsize=10)
    save_fig(fig, "fig4_correlations_DV", outdir)
    return corr_df.drop(columns=[], errors="ignore")


# ── fig5: pairwise metric correlation matrix ──────────────────────────────────
def plot_metric_correlation_matrix(final, metrics, outdir="."):
    if len(final) < 2:
        save_placeholder_fig(
            "fig5_metric_correlation_matrix",
            "Pairwise metric correlation",
            "At least two final-state cases are required to compute a correlation matrix.",
            outdir,
        )
        return

    corr = final[metrics].corr(method="spearman")
    n = len(metrics)
    sz = max(6, n * 0.55)
    fig, ax = plt.subplots(figsize=(sz, sz * 0.88))

    mask = np.eye(n, dtype=bool)          # hide the diagonal (all 1.0)
    sns.heatmap(
        corr, mask=mask, ax=ax,
        cmap="RdBu_r", center=0, vmin=-1, vmax=1,
        square=True, linewidths=0.25, linecolor="white",
        annot=(n <= 12), fmt=".2f", annot_kws={"size": 7},
        cbar_kws={"shrink": 0.65, "label": "Spearman ρ"},
        xticklabels=[mlabel(m) for m in metrics],
        yticklabels=[mlabel(m) for m in metrics],
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)
    ax.set_title("Pairwise Spearman correlation — all metrics (final state)", fontsize=10, pad=8)
    fig.tight_layout()
    save_fig(fig, "fig5_metric_correlation_matrix", outdir)


# ── fig6: PCA ─────────────────────────────────────────────────────────────────
def plot_pca(final, metrics, outdir="."):
    if len(final) < 2 or len(metrics) < 2:
        save_placeholder_fig(
            "fig6_pca",
            "PCA of final-state morphology",
            "PCA requires at least two final-state cases and two metric columns.",
            outdir,
        )
        pd.DataFrame(columns=["PC", "explained_variance_ratio", "cumulative"]).to_csv(
            f"{outdir}/pca_variance.csv", index=False
        )
        pd.DataFrame(columns=metrics).to_csv(f"{outdir}/pca_loadings.csv")
        print("  PCA skipped: insufficient final-state cases or metrics")
        return None, None, np.array([])

    X = final[metrics].values
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)

    n_comp = min(len(metrics), len(Xs), 6)
    pca = PCA(n_components=n_comp)
    scores = pca.fit_transform(Xs)
    loadings = pca.components_
    ev = pca.explained_variance_ratio_

    D_vals = final["D"].values
    V_vals = final["V"].values

    fig = plt.figure(figsize=(14, 4.8))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.5, 1.5, 1.4], wspace=0.35)

    for panel, (param_vals, param_name, cmap_name) in enumerate([
        (D_vals, "D", "plasma"),
        (V_vals, "V", "coolwarm"),
    ]):
        ax = fig.add_subplot(gs[panel])
        sc = ax.scatter(scores[:, 0], scores[:, 1],
                        c=param_vals, cmap=cmap_name, s=20, alpha=0.85, linewidths=0)
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label=param_name)
        cb.ax.tick_params(labelsize=7)
        ax.set_xlabel(f"PC1 ({ev[0]*100:.1f}%)", fontsize=9)
        ax.set_ylabel(f"PC2 ({ev[1]*100:.1f}%)", fontsize=9)
        ax.set_title(f"PCA scores — coloured by {param_name}", fontsize=9)

    # Loadings: PC1 and PC2 only
    ax3 = fig.add_subplot(gs[2])
    yl = np.arange(len(metrics))
    ax3.barh(yl - 0.18, loadings[0], height=0.34,
             label=f"PC1 ({ev[0]*100:.1f}%)", color="#2c7bb6", alpha=0.85)
    ax3.barh(yl + 0.18, loadings[1], height=0.34,
             label=f"PC2 ({ev[1]*100:.1f}%)", color="#d7191c", alpha=0.85)
    ax3.set_yticks(yl)
    ax3.set_yticklabels([mlabel(m) for m in metrics], fontsize=7)
    ax3.axvline(0, color="k", linewidth=0.6)
    ax3.set_xlabel("Loading", fontsize=9)
    ax3.set_title("PC1 / PC2 loadings", fontsize=9)
    ax3.legend(frameon=False, fontsize=7)

    fig.suptitle("PCA of final-state morphology (all metrics, StandardScaler)", fontsize=10, y=1.02)
    save_fig(fig, "fig6_pca", outdir)

    pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(len(ev))],
        "explained_variance_ratio": ev,
        "cumulative": np.cumsum(ev),
    }).to_csv(f"{outdir}/pca_variance.csv", index=False)

    pd.DataFrame(
        loadings, columns=metrics,
        index=[f"PC{i+1}" for i in range(n_comp)]
    ).to_csv(f"{outdir}/pca_loadings.csv")

    print(f"  PCA variance explained: {[f'{v*100:.1f}%' for v in ev]}")
    return pca, scores, ev


# ── fig7: D×V interaction lines ───────────────────────────────────────────────
def plot_interaction(final, metrics, outdir="."):
    n = len(metrics)
    D_vals = sorted(final["D"].unique())
    cmap = plt.cm.plasma
    colors = palette_from_values(D_vals, cmap)

    fig, axes = plt.subplots(1, n, figsize=(4.3 * n, 3.8),
                             constrained_layout=True, squeeze=False)
    axes = axes[0]

    for ax, metric in zip(axes, metrics):
        for d, color in zip(D_vals, colors):
            sub = final[final["D"] == d].sort_values("V")
            ax.plot(sub["V"], sub[metric], color=color, linewidth=1.0,
                    marker="o", markersize=2.5, alpha=0.85)
        ax.set_xlabel("V", fontsize=9)
        ax.set_title(mlabel(metric), fontsize=9, pad=4)

    axes[0].set_ylabel("Final-state value", fontsize=9)

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=scalar_norm(D_vals))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=list(axes), fraction=0.015, pad=0.02, label="D")
    cb.ax.tick_params(labelsize=7)

    fig.suptitle("D × V interaction at final state (each line = one D level)", fontsize=10)
    save_fig(fig, "fig7_interaction", outdir)


# ── fig8: marginal effects ────────────────────────────────────────────────────
def plot_marginal_effects(final, metrics, outdir="."):
    n = len(metrics)
    fig, axes = plt.subplots(2, n, figsize=(4.1 * n, 6.0),
                             constrained_layout=True, squeeze=False)

    for col, metric in enumerate(metrics):
        for row, (param, color) in enumerate([("D", "#2c7bb6"), ("V", "#d7191c")]):
            ax = axes[row][col]
            agg = final.groupby(param)[metric].agg(["mean", "std"])
            ax.errorbar(agg.index, agg["mean"], yerr=agg["std"],
                        fmt="o-", color=color, linewidth=1.2, markersize=4,
                        capsize=3, elinewidth=0.8)
            ax.set_xlabel(param, fontsize=9)
            if col == 0:
                ax.set_ylabel(f"Mean ± SD over {'V' if param == 'D' else 'D'}", fontsize=9)
        axes[0][col].set_title(mlabel(metric), fontsize=9)

    fig.suptitle("Marginal effects of D and V on key metrics (final state)", fontsize=10)
    save_fig(fig, "fig8_marginal_effects", outdir)


# ── fig9: simulation length heatmap ──────────────────────────────────────────
def plot_series_lengths(df, outdir="."):
    lengths = df.groupby(["D", "V"])["timestep"].max().reset_index()
    lengths.columns = ["D", "V", "max_timestep"]
    D_vals = sorted(lengths["D"].unique())
    V_vals = sorted(lengths["V"].unique())
    grid = (
        lengths.pivot(index="V", columns="D", values="max_timestep")
        .reindex(index=V_vals[::-1], columns=D_vals)
    )

    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    im = ax.imshow(grid.values, aspect="auto", cmap="YlOrRd",
                   extent=[-0.5, len(D_vals) - 0.5, -0.5, len(V_vals) - 0.5])
    ax.set_xticks(range(len(D_vals)))
    ax.set_xticklabels([f"{v:g}" for v in D_vals], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(V_vals)))
    ax.set_yticklabels([f"{v:g}" for v in V_vals[::-1]], fontsize=8)
    ax.set_xlabel("D", fontsize=9)
    ax.set_ylabel("V", fontsize=9)
    ax.set_title("Simulation length — max timestep per case", fontsize=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Max timestep")
    cb.ax.tick_params(labelsize=8)
    fig.tight_layout()
    save_fig(fig, "fig9_series_lengths", outdir)


# ── fig10: information gain scatter ──────────────────────────────────────────
def plot_info_gain(final, metrics, outdir="."):
    """
    For each metric (excluding total_length_px as size reference):
      x = |Spearman(metric, total_length_px)|  →  redundancy with bulk size
      y = max(|Spearman(metric, D)|, |Spearman(metric, V)|)  →  D/V sensitivity

    Metrics in the top-left quadrant capture D/V variation independently of size.
    """
    size_ref = "total_length_px"
    if len(final) < 2:
        pd.DataFrame(
            columns=["metric", "rho_D", "rho_V", "abs_rho_size", "max_abs_DV", "size_redundant", "informative"]
        ).to_csv(f"{outdir}/info_gain_table.csv", index=False)
        save_placeholder_fig(
            "fig10_info_gain",
            "Information gain beyond size",
            "At least two final-state cases are required to compare D/V sensitivity against size redundancy.",
            outdir,
        )
        print("\n  Information-gain analysis skipped: insufficient final-state cases")
        return pd.DataFrame()

    records = []
    for m in metrics:
        if m == size_ref:
            continue
        r_D, _ = stats.spearmanr(final["D"], final[m])
        r_V, _ = stats.spearmanr(final["V"], final[m])
        r_sz, _ = stats.spearmanr(final[size_ref], final[m])
        records.append({
            "metric":        m,
            "rho_D":         round(r_D, 4),
            "rho_V":         round(r_V, 4),
            "abs_rho_size":  round(abs(r_sz), 4),
            "max_abs_DV":    round(max(abs(r_D), abs(r_V)), 4),
        })

    info = pd.DataFrame(records)
    info["size_redundant"] = info["abs_rho_size"] > 0.80
    info["informative"]    = (info["max_abs_DV"] > 0.30) & ~info["size_redundant"]
    info.to_csv(f"{outdir}/info_gain_table.csv", index=False)

    # color scheme
    def pt_color(row):
        if row["informative"]:
            return "#d7191c"       # informative beyond size
        if row["size_redundant"]:
            return "#2c7bb6"       # size-correlated
        return "#888888"           # neither

    fig, ax = plt.subplots(figsize=(8, 6))

    for _, row in info.iterrows():
        c = pt_color(row)
        ax.scatter(row["abs_rho_size"], row["max_abs_DV"],
                   color=c, s=60, zorder=3, alpha=0.9)
        ax.annotate(mlabel(row["metric"]),
                    (row["abs_rho_size"], row["max_abs_DV"]),
                    fontsize=7, xytext=(5, 3), textcoords="offset points")

    ax.axvline(0.80, color="#2c7bb6", linestyle="--", linewidth=0.9, alpha=0.65)
    ax.axhline(0.30, color="#888888", linestyle="--", linewidth=0.9, alpha=0.65)
    ax.set_xlabel(f"|Spearman ρ| with {mlabel(size_ref)}   (size redundancy)", fontsize=9)
    ax.set_ylabel("max(|ρ with D|, |ρ with V|)   (D/V sensitivity)", fontsize=9)
    ax.set_title(
        "Information gain beyond size\n"
        "Red = informative  |  Blue = size-redundant  |  Grey = low D/V sensitivity",
        fontsize=9,
    )
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.text(0.88, 0.05, "size-\nredundant", color="#2c7bb6", fontsize=7,
            ha="center", va="bottom", alpha=0.8)
    ax.text(0.08, 0.92, "informative\nbeyond size", color="#d7191c", fontsize=7,
            ha="center", va="top", alpha=0.8)
    fig.tight_layout()
    save_fig(fig, "fig10_info_gain", outdir)

    print("\n  Informative metrics (max|ρ_DV| > 0.30 and |ρ_size| < 0.80):")
    subset = info[info["informative"]].sort_values("max_abs_DV", ascending=False)
    if len(subset):
        for _, r in subset.iterrows():
            print(f"    {mlabel(r['metric']):35s}  ρ_D={r['rho_D']:+.3f}  "
                  f"ρ_V={r['rho_V']:+.3f}  ρ_size={r['abs_rho_size']:.3f}")
    else:
        print("    (none at current thresholds)")

    print("\n  Size-redundant metrics (|ρ_size| > 0.80):")
    redundant = info[info["size_redundant"]].sort_values("abs_rho_size", ascending=False)
    for _, r in redundant.iterrows():
        print(f"    {mlabel(r['metric']):35s}  ρ_size={r['abs_rho_size']:.3f}")

    return info


# ── main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze morphometrics and generate figures.")
    parser.add_argument("--input", default="morphometrics.csv", help="Input morphometrics CSV path")
    parser.add_argument("--output-dir", default=".", help="Directory for generated figures and tables")
    cli_args = parser.parse_args()

    input_csv = Path(cli_args.input).expanduser().resolve()
    outdir = Path(cli_args.output_dir).expanduser().resolve()

    if not input_csv.is_file():
        raise SystemExit(f"ERROR: morphometrics CSV not found: {input_csv}")

    outdir.mkdir(parents=True, exist_ok=True)
    outdir_str = str(outdir)

    print("Loading data...")
    df    = load_data(str(input_csv))
    final = final_state(df)

    all_metrics = detect_metrics(df)
    km  = select_key_metrics(all_metrics, n=8)
    km5 = km[:5]                                    # time-course: 5 panels
    im  = [m for m in INTERACTION_METRICS if m in all_metrics]

    print(f"  {len(df)} rows | {df['case'].nunique()} cases | "
          f"{len(all_metrics)} metric columns")
    print(f"  D: {sorted(df['D'].unique())}")
    print(f"  V: {sorted(df['V'].unique())}")
    print(f"  Key metrics (heatmaps/PCA): {km}")
    print(f"  Time-course metrics       : {km5}")
    print(f"  Interaction metrics       : {im}")
    print()

    # summary CSV
    cols_out = ["case", "D", "V", "timestep"] + all_metrics
    (final[cols_out]
     .rename(columns={"timestep": "final_timestep"})
     .sort_values(["D", "V"])
     .to_csv(f"{outdir_str}/summary_final_state.csv", index=False))
    print("  Saved summary_final_state.csv")

    print("\nGenerating figures...")
    plot_final_heatmaps(final, km, outdir_str)
    plot_timecourse(df, km5, by="D", outdir=outdir_str)
    plot_timecourse(df, km5, by="V", outdir=outdir_str)
    plot_correlations_DV(final, all_metrics, outdir_str)
    plot_metric_correlation_matrix(final, all_metrics, outdir_str)
    plot_pca(final, all_metrics, outdir_str)
    plot_interaction(final, im, outdir_str)
    plot_marginal_effects(final, im, outdir_str)
    plot_series_lengths(df, outdir_str)
    plot_info_gain(final, all_metrics, outdir_str)

    print("\nDone.")
