from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\刘子轩\Desktop\nature")
FIG_DIR = ROOT / "09_figures_v1"
FIG_DIR.mkdir(parents=True, exist_ok=True)

BENCHMARK = ROOT / "02_benchmark_v0_1" / "benchmark_region_evidence_v0_1.csv"
AREA_SUMMARY = ROOT / "03_exposure_closure" / "chao_phraya_area_weighted_exposure_censoring" / "chao_phraya_area_weighted_exposure_summary.csv"
OSM_SUMMARY = ROOT / "03_exposure_closure" / "chao_phraya_osm_exposure_censoring" / "chao_phraya_osm_exposure_censoring_summary.csv"


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
    fig.savefig(FIG_DIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.tiff", dpi=600, bbox_inches="tight")


def add_label(fig: mpl.figure.Figure, ax: mpl.Axes, label: str) -> None:
    pos = ax.get_position()
    fig.text(pos.x0 - 0.02, pos.y1 + 0.01, label, fontsize=9.5, fontweight="bold")


def main() -> int:
    benchmark = pd.read_csv(BENCHMARK).sort_values("landcover_adjusted_or", ascending=False).copy()
    area = pd.read_csv(AREA_SUMMARY).iloc[0]
    osm = pd.read_csv(OSM_SUMMARY)

    transport = osm[osm["exposure_category"] == "transport_total"].copy()
    transport["year"] = transport["pair"].str.slice(0, 4).astype(int)
    transport = transport.sort_values("year")
    transport["visible_km"] = transport["transport_km_strong_subsidence"] - transport["transport_km_hidden_strong_subsidence"]

    fig = plt.figure(figsize=(7.6, 5.8), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.02, 1.0], hspace=0.34, wspace=0.28)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # Panel A: stacked exposure bars
    labels = ["Population (M)", "Built-up (km$^2$)"]
    total_vals = [
        float(area["total_population_vlm_grid"]) / 1_000_000.0,
        float(area["total_builtup_km2_vlm_grid"]),
    ]
    hidden_vals = [
        float(area["strong_sub_5mm_population_not_majority"]) / 1_000_000.0,
        float(area["strong_sub_5mm_builtup_not_majority_km2"]),
    ]
    visible_vals = [t - h for t, h in zip(total_vals, hidden_vals)]
    x = np.arange(len(labels))
    ax_a.bar(x, visible_vals, color="#5E81AC", width=0.62, label="Visible")
    ax_a.bar(x, hidden_vals, bottom=visible_vals, color="#D08770", width=0.62, label="Hidden")
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(labels)
    ax_a.set_ylabel("Exposure")
    ax_a.set_title("A  Strong-motion exposure split into visible and hidden parts", loc="left", fontweight="bold")
    ax_a.grid(axis="y", color="#E1E5EE", linewidth=0.7, zorder=0)
    ax_a.legend(loc="upper right", fontsize=6.5)
    ax_a.text(
        0,
        total_vals[0] * 1.02,
        f"hidden {float(area['strong_sub_5mm_population_not_majority_fraction']):.2f}",
        ha="center",
        va="bottom",
        fontsize=6.2,
    )
    ax_a.text(
        1,
        total_vals[1] * 1.02,
        f"hidden {float(area['strong_sub_5mm_builtup_not_majority_fraction']):.2f}",
        ha="center",
        va="bottom",
        fontsize=6.2,
    )

    # Panel B: regional underestimation factor
    signal_order = {"strong_positive": 0, "positive": 1, "inconclusive": 2}
    benchmark["order"] = benchmark["observability_bias_signal"].map(signal_order).fillna(3).astype(int)
    benchmark = benchmark.sort_values(["order", "landcover_adjusted_or"], ascending=[True, False])
    colors = benchmark["observability_bias_signal"].map(
        {"strong_positive": "#D08770", "positive": "#A3BE8C", "inconclusive": "#8FA3BF"}
    )
    y = np.arange(len(benchmark))
    ax_b.barh(y, benchmark["landcover_adjusted_or"], color=colors, height=0.65)
    ax_b.errorbar(
        benchmark["landcover_adjusted_or"],
        y,
        xerr=[
            benchmark["landcover_adjusted_or"] - benchmark["landcover_adjusted_boot_q025"],
            benchmark["landcover_adjusted_boot_q975"] - benchmark["landcover_adjusted_or"],
        ],
        fmt="none",
        ecolor="#3B4252",
        elinewidth=0.8,
        capsize=2.5,
    )
    ax_b.axvline(1.0, color="#444444", lw=0.9, ls="--")
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(benchmark["region"])
    ax_b.set_xlabel("Land-cover-adjusted OR")
    ax_b.set_title("B  Regional bias factors are not near unity", loc="left", fontweight="bold")
    ax_b.set_xlim(0, max(benchmark["landcover_adjusted_boot_q975"]) * 1.1)
    ax_b.grid(axis="x", color="#E1E5EE", linewidth=0.7)
    ax_b.invert_yaxis()

    # Panel C: hidden fractions by exposure layer
    layers = ["Population", "Built-up", "Transport"]
    not_majority = [
        float(area["strong_sub_5mm_population_not_majority_fraction"]),
        float(area["strong_sub_5mm_builtup_not_majority_fraction"]),
        float(transport["hidden_strong_subsidence_transport_fraction"].mean()),
    ]
    x2 = np.arange(len(layers))
    ax_c.bar(x2, not_majority, width=0.55, color="#5E81AC")
    ax_c.set_xticks(x2)
    ax_c.set_xticklabels(layers)
    ax_c.set_ylim(0, 1.0)
    ax_c.set_ylabel("Fraction")
    ax_c.set_title("C  Hidden exposure remains substantial across layers", loc="left", fontweight="bold")
    ax_c.grid(axis="y", color="#E1E5EE", linewidth=0.7)
    for xi, val in zip(x2, not_majority):
        ax_c.text(xi, val + 0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=6.3)

    # Panel D: transport hidden by year
    years = transport["year"].to_numpy()
    visible = transport["visible_km"].to_numpy()
    hidden = transport["transport_km_hidden_strong_subsidence"].to_numpy()
    ax_d.bar(years, visible, color="#4C566A", width=0.62, label="Visible")
    ax_d.bar(years, hidden, bottom=visible, color="#A3BE8C", width=0.62, label="Hidden")
    ax_d.set_xlabel("Year")
    ax_d.set_ylabel("Transport length (km)")
    ax_d.set_title("D  Transport exposure hidden year by year", loc="left", fontweight="bold")
    ax_d.grid(axis="y", color="#E1E5EE", linewidth=0.7)
    ax_d.legend(loc="upper right", fontsize=6.5)

    fig.text(
        0.5,
        0.01,
        f"Chao Phraya hidden strong-motion fractions: population {float(area['strong_sub_5mm_population_not_majority_fraction']):.3f}, "
        f"built-up {float(area['strong_sub_5mm_builtup_not_majority_fraction']):.3f}.",
        ha="center",
        fontsize=6.2,
        color="#3B4252",
    )
    fig.suptitle(
        "From observability censoring to exposure underestimation",
        y=0.995,
        fontsize=9.2,
        fontweight="bold",
    )

    save_pub(fig, "fig3_exposure_translation")
    meta = {
        "benchmark": str(BENCHMARK),
        "area_summary": str(AREA_SUMMARY),
        "osm_summary": str(OSM_SUMMARY),
        "output_dir": str(FIG_DIR),
        "regions": int(len(benchmark)),
    }
    (FIG_DIR / "fig3_exposure_translation_meta.json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    plt.close(fig)
    print(json.dumps(meta, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
