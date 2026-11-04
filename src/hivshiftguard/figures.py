"""Generate publication-style, fully reproducible figures after real-data analysis."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

COLORS = {"ink": "#15263E", "blue": "#186AA9", "teal": "#0C8D83", "orange": "#DC8C37", "red": "#BF5353"}


def setup():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 220, "figure.dpi": 120})


def save_both(fig, folder, basename):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / (basename + ".png"), bbox_inches="tight", facecolor="white")
    fig.savefig(folder / (basename + ".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def architecture(folder):
    setup()
    fig, ax = plt.subplots(figsize=(13.8, 7))
    fig.patch.set_facecolor("#FFFFFF")
    ax.set_xlim(0, 13.8)
    ax.set_ylim(0, 7)
    ax.axis("off")
    ax.text(6.9, 6.72, "HIV-ShiftGuard | research architecture", ha="center", fontsize=17, weight="bold", color=COLORS["ink"])
    ax.text(6.9, 6.34, "In-vitro HIV-1 fold susceptibility; not treatment selection or cure discovery", ha="center", fontsize=9, color=COLORS["red"])
    boxes = [
        (0.35, 4.9, 3.25, 1.0, "Stanford HIVDB", "PI  |  NRTI  |  NNRTI  |  INI", COLORS["blue"]),
        (4.05, 4.9, 3.25, 1.0, "Data integrity", "Parse, deduplicate, hash & verify", COLORS["teal"]),
        (7.75, 4.9, 5.7, 1.0, "Leakage-audited grouped partitions", "Training  /  validation  /  calibration  /  untouched test", COLORS["blue"]),
        (0.35, 2.98, 3.25, 1.0, "Mutation representation", "Residues, positions & mutation sets", COLORS["teal"]),
        (4.05, 2.98, 3.25, 1.0, "Baseline and ablation", "Median, ridge, position-only trees", COLORS["orange"]),
        (7.75, 2.98, 5.7, 1.0, "Primary ensemble", "Mutation forest + validation-tuned ridge fusion", COLORS["blue"]),
        (0.35, 1.08, 3.25, 1.0, "Novelty estimation", "Nearest training genotype (Jaccard)", COLORS["teal"]),
        (4.05, 1.08, 3.25, 1.0, "Conformal uncertainty", "Fixed target; standard vs novelty-scaled", COLORS["orange"]),
        (7.75, 1.08, 5.7, 1.0, "Evaluation & reporting", "Unseen genotypes, coverage, abstention, tables, SVG", COLORS["blue"]),
    ]
    for x, y, w, h, head, sub, col in boxes:
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.18",
                               linewidth=1.4, edgecolor=col, facecolor="#F3F8FC")
        ax.add_patch(patch)
        ax.text(x + w / 2, y + .65, head, ha="center", va="center", fontsize=11.5, weight="bold", color=COLORS["ink"])
        ax.text(x + w / 2, y + .30, sub, ha="center", va="center", fontsize=8.2, color="#4A5B6B")
    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=17,
                                     linewidth=1.65, color="#667E91"))
    for y in (5.4, 3.48, 1.58):
        arrow(3.65, y, 3.96, y)
        arrow(7.35, y, 7.65, y)
    arrow(9.9, 4.8, 9.9, 4.1)
    arrow(9.9, 2.88, 9.9, 2.16)
    arrow(1.97, 4.8, 1.97, 4.1)
    arrow(1.97, 2.88, 1.97, 2.16)
    save_both(fig, folder, "Figure_01_Research_Architecture")


def build_figures(folder, metrics_df, preds_df):
    setup()
    architecture(folder)
    if metrics_df.empty or preds_df.empty:
        return
    order = ["median_baseline", "linear_mutations", "position_only_ablation", "mutation_forest", "shiftguard_ensemble"]
    labels = ["Median", "Mutation ridge", "Position only", "Mutation forest", "ShiftGuard"]
    means = [metrics_df.loc[metrics_df.model == k, "rmse"].mean() for k in order]
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    bars = ax.barh(labels, means, color=["#A4B3C0", "#7E9DB5", "#D1AA68", "#7EA6B6", COLORS["teal"]])
    ax.invert_yaxis(); ax.set_xlabel("Mean RMSE of log2 fold susceptibility (lower is better)")
    ax.set_title("Model and ablation comparison — real held-out predictions")
    for bar, value in zip(bars, means):
        ax.text(bar.get_width() + max(means) * .008, bar.get_y() + bar.get_height() / 2, f"{value:.3f}", va="center", fontsize=9)
    ax.set_xlim(0, max(means) * 1.17)
    save_both(fig, folder, "Figure_02_Model_Ablation_RMSE")

    summary = (metrics_df[metrics_df.model == "shiftguard_ensemble"].groupby(["gene", "drug"])
               .agg(rmse=("rmse", "mean"), experiments=("seed", "count")).reset_index()
               .sort_values("rmse", ascending=False).head(17))
    if not summary.empty:
        fig, ax = plt.subplots(figsize=(10, max(4.2, .33 * len(summary) + 1.5)))
        ax.barh(summary["gene"] + ":" + summary["drug"], summary["rmse"], color=COLORS["blue"])
        ax.invert_yaxis(); ax.set_xlabel("Test RMSE, log2 fold susceptibility")
        ax.set_title("Drug-level variation (largest errors shown)")
        save_both(fig, folder, "Figure_03_Drugwise_Test_RMSE")

    fig, ax = plt.subplots(figsize=(7, 6.2))
    actual = preds_df["observed_log2fc"].to_numpy()
    pred = preds_df["predicted_log2fc"].to_numpy()
    ax.scatter(actual, pred, s=11, alpha=.20, rasterized=True, color=COLORS["blue"])
    bounds = [float(min(actual.min(), pred.min())), float(max(actual.max(), pred.max()))]
    ax.plot(bounds, bounds, ls="--", color=COLORS["red"], linewidth=1.2, label="Perfect agreement")
    ax.set_xlabel("Observed log2 fold susceptibility")
    ax.set_ylabel("Predicted log2 fold susceptibility")
    ax.set_title("Untouched test: prediction vs observation (all seeds)")
    ax.legend(loc="upper left")
    save_both(fig, folder, "Figure_04_Predicted_vs_Observed")

    fig, ax = plt.subplots(figsize=(8, 4.8))
    y = preds_df["observed_log2fc"].to_numpy()
    methods = [("Standard", preds_df.interval_lower_unscaled.to_numpy(), preds_df.interval_upper_unscaled.to_numpy(), COLORS["blue"]),
               ("Novelty-scaled", preds_df.interval_lower_novelty.to_numpy(), preds_df.interval_upper_novelty.to_numpy(), COLORS["teal"])]
    cover = [float(np.mean((y >= lower) & (y <= upper))) for _, lower, upper, _ in methods]
    ax.bar([x[0] for x in methods], cover, color=[x[3] for x in methods])
    ax.axhline(.90, color=COLORS["red"], ls="--", linewidth=1.4, label="Nominal target 90%")
    for i, value in enumerate(cover):
        ax.text(i, value + .012, f"{value:.1%}", ha="center", fontsize=11, weight="bold")
    ax.set_ylim(0, min(1.08, max(.94, max(cover) + .08)))
    ax.set_ylabel("Observed test-interval coverage")
    ax.set_title("Conformal interval calibration — no test tuning")
    ax.legend(loc="lower right")
    save_both(fig, folder, "Figure_05_Interval_Coverage")

    residuals = np.abs(y - preds_df["predicted_log2fc"].to_numpy())
    novelty = preds_df["novelty_jaccard"].to_numpy()
    fig, ax = plt.subplots(figsize=(8.3, 4.8))
    ax.scatter(novelty, residuals, s=14, alpha=.2, color=COLORS["teal"], rasterized=True)
    ax.set_xlabel("Nearest-training-genotype Jaccard novelty")
    ax.set_ylabel("Absolute test error (log2 fold)")
    ax.set_title("Mutation novelty and generalization error")
    save_both(fig, folder, "Figure_06_Error_vs_Novelty")

    fig, ax = plt.subplots(figsize=(8.5, 5.1))
    ranks = preds_df.sort_values("novelty_jaccard").copy()
    ranks["bin"] = pd.qcut(ranks["novelty_jaccard"].rank(method="first"), q=5, labels=False, duplicates="drop")
    across = ranks.groupby("bin", as_index=False).apply(
        lambda d: pd.Series({"mae": float(np.mean(abs(d.observed_log2fc - d.predicted_log2fc))),
                             "novelty": float(d.novelty_jaccard.mean())}), include_groups=False)
    ax.plot(across["novelty"], across["mae"], "o-", lw=2, color=COLORS["orange"])
    ax.set_xlabel("Mean Jaccard novelty by equal-size quintile")
    ax.set_ylabel("Mean absolute error (log2 fold)")
    ax.set_title("Performance across genotype novelty strata")
    save_both(fig, folder, "Figure_07_Novelty_Quintiles")
