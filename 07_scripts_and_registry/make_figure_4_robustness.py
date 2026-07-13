from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle


ROOT = Path(r"C:\Users\刘子轩\Desktop\nature")
INPUT = ROOT / "03_exposure_closure" / "chao_phraya_robustness_grid" / "chao_phraya_robustness_grid.csv"
OUTDIR = ROOT / "09_figures_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    }
)


def save_pub(fig: mpl.figure.Figure, stem: str) -> None:
    fig.savefig(OUTDIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.tiff", dpi=600, bbox_inches="tight")


def fmt(x: float) -> str:
    return f"{x:.2f}"


def build_heatmap(ax: plt.Axes, df: pd.DataFrame) -> None:
    pivot = df.pivot_table(index="strong_threshold", columns="threshold", values="strong_vs_nonstrong_not_majority_odds_ratio")
    strong_levels = sorted(pivot.index.tolist())
    thresh_levels = sorted(pivot.columns.tolist())
    values = pivot.loc[strong_levels, thresh_levels].to_numpy()
    norm = TwoSlopeNorm(vmin=np.nanmin(values), vcenter=1.0, vmax=np.nanmax(values))
    im = ax.imshow(values, cmap="RdYlBu_r", norm=norm, aspect="auto")
    ax.set_xticks(np.arange(len(thresh_levels)))
    ax.set_xticklabels([fmt(x) for x in thresh_levels])
    ax.set_yticks(np.arange(len(strong_levels)))
    ax.set_yticklabels([fmt(x) for x in strong_levels])
    ax.set_xlabel("LiCSAR observability threshold")
    ax.set_ylabel("Strong-subsidence threshold (mm/yr)")
    ax.set_title("A  Threshold sensitivity of censoring signal", loc="left", fontweight="bold")
    for i, strong in enumerate(strong_levels):
        for j, thr in enumerate(thresh_levels):
            row = df[(df["strong_threshold"] == strong) & (df["threshold"] == thr)].iloc[0]
            val = row["strong_vs_nonstrong_not_majority_odds_ratio"]
            txt_color = "white" if val > (np.nanmin(values) + np.nanmax(values)) / 2 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color=txt_color)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Odds ratio")


def build_block_intervals(ax: plt.Axes, df: pd.DataFrame) -> None:
    colors = {
        0.2: "#5E81AC",
        0.3: "#A3BE8C",
        0.4: "#D08770",
    }
    x_pos = np.arange(len(sorted(df["strong_threshold"].unique())))
    labels = [fmt(x) for x in sorted(df["strong_threshold"].unique())]
    width = 0.22
    for idx, thr in enumerate(sorted(df["threshold"].unique())):
        sub = df[(df["threshold"] == thr) & (df["block_size"] == 10)].sort_values("strong_threshold")
        if sub.empty:
            continue
        xs = x_pos + (idx - 1) * width
        y = sub["strong_vs_nonstrong_not_majority_odds_ratio"].to_numpy()
        lo = sub["strong_vs_nonstrong_not_majority_odds_ratio_q025"].to_numpy()
        hi = sub["strong_vs_nonstrong_not_majority_odds_ratio_q975"].to_numpy()
        ax.errorbar(
            xs,
            y,
            yerr=[y - lo, hi - y],
            fmt="o-",
            lw=1.2,
            ms=4,
            color=colors[thr],
            capsize=2.8,
            label=f"threshold {thr}",
        )
    ax.axhline(1.0, color="#444444", lw=0.9, ls="--")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Strong-subsidence threshold (mm/yr)")
    ax.set_ylabel("Odds ratio")
    ax.set_title("B  Block bootstrap intervals at 10x10 blocks", loc="left", fontweight="bold")
    ax.legend(title="Observability threshold", loc="upper left", fontsize=6, title_fontsize=6)


def build_line_plot(ax: plt.Axes, df: pd.DataFrame) -> None:
    colors = {
        0.2: "#5E81AC",
        0.3: "#A3BE8C",
        0.4: "#D08770",
    }
    block_size = 10
    for thr in sorted(df["threshold"].unique()):
        sub = df[(df["threshold"] == thr) & (df["block_size"] == block_size)].sort_values("strong_threshold")
        x = sub["strong_threshold"].to_numpy()
        y = sub["strong_not_majority_fraction"].to_numpy()
        lo = sub["strong_not_majority_fraction_q025"].to_numpy()
        hi = sub["strong_not_majority_fraction_q975"].to_numpy()
        ax.plot(x, y, marker="o", lw=1.4, color=colors[thr], label=f"threshold {thr}")
        ax.fill_between(x, lo, hi, color=colors[thr], alpha=0.18)
    ax.set_xlabel("Strong-subsidence threshold (mm/yr)")
    ax.set_ylabel("Strong cells not-majority observable")
    ax.set_ylim(0.45, 0.9)
    ax.set_title("C  Strong-motion censoring survives moderate and strong thresholds", loc="left", fontweight="bold")
    ax.legend(title="Observability threshold", loc="upper right", fontsize=6, title_fontsize=6)


def add_callout(fig: mpl.figure.Figure, ax: plt.Axes) -> None:
    bbox = ax.get_position()
    x = bbox.x1 - 0.27
    y = bbox.y0 + 0.035
    box = Rectangle(
        (x, y),
        0.245,
        0.085,
        transform=fig.transFigure,
        facecolor="#F6F7FB",
        edgecolor="#C9CEDA",
        linewidth=0.8,
        zorder=3,
    )
    fig.patches.append(box)
    fig.text(
        x + 0.010,
        y + 0.056,
        "Interpretation",
        fontsize=7.2,
        fontweight="bold",
        ha="left",
        va="center",
    )
    fig.text(
        x + 0.010,
        y + 0.032,
        "The censoring signal remains above 1 for the lead case",
        fontsize=6.0,
        ha="left",
        va="center",
    )
    fig.text(
        x + 0.010,
        y + 0.013,
        "at 5 mm/yr and survives block resampling.",
        fontsize=6.0,
        ha="left",
        va="center",
    )


def main() -> int:
    df = pd.read_csv(INPUT)
    fig = plt.figure(figsize=(7.2, 5.5), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.95], width_ratios=[1.0, 1.0], hspace=0.38, wspace=0.28)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    build_heatmap(ax_a, df)
    build_block_intervals(ax_b, df)
    build_line_plot(ax_c, df)
    add_callout(fig, ax_c)

    fig.suptitle(
        "Chao Phraya robustness grid: observability censoring survives threshold and block-size sensitivity",
        y=0.995,
        fontsize=9.2,
        fontweight="bold",
    )
    save_pub(fig, "fig4_chao_phraya_robustness_grid")

    meta = {
        "input": str(INPUT),
        "output_dir": str(OUTDIR),
        "rows": int(len(df)),
    }
    (OUTDIR / "fig4_chao_phraya_robustness_grid_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    plt.close(fig)
    print(json.dumps(meta, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
