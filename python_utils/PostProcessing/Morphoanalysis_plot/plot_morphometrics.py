
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

df = pd.read_csv("morphometrics.csv")

def parse_case(s):
    m = re.match(r"D(?P<D>[0-9.]+)_V(?P<V>[0-9]+(?:p[0-9]+)?|[0-9.]+)$", s)
    return float(m.group("D")), float(m.group("V").replace("p", "."))

df[["D", "V"]] = pd.DataFrame(df["case"].apply(parse_case).tolist(), index=df.index)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

fig, axes = plt.subplots(1, 3, figsize=(11, 3.2), constrained_layout=True)
for ax, (key, title) in zip(
    axes,
    [("total_length_px", "Total length (px)"),
     ("n_tips", "Terminal tips"),
     ("n_branchpoints", "Branch points")]
):
    for d_value, sub in df.groupby("D"):
        summary = sub.groupby("timestep")[key].mean()
        ax.plot(summary.index, summary.values, linewidth=1.1, label=f"D={d_value:g}")
    ax.set_title(title)
    ax.set_xlabel("Timestep")
axes[0].set_ylabel("Mean value across V")
axes[-1].legend(frameon=False, ncol=3, fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1.0))
fig.savefig("morphology_timecourse_by_D.pdf", bbox_inches="tight")
fig.savefig("morphology_timecourse_by_D.png", dpi=300, bbox_inches="tight")
