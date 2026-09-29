#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from common import load_config, resolve

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    root = Path.cwd()
    cfg = load_config(args.config)
    sample_table = resolve(root, cfg["samples"])
    color_table = resolve(root, cfg["colors"])
    results = resolve(root, cfg["results"])
    cand_dir = results / "03_mpileup_candidates"
    summary_dir = results / "04_summary"
    fig_dir = results / "05_figures"
    summary_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    samples = pd.read_csv(sample_table, sep=r"\s+", engine="python")[["sample", "group"]]
    samples["sample"] = samples["sample"].astype(str)
    samples["group"] = samples["group"].astype(str)
    sample_order = samples["sample"].tolist()
    group_order = list(dict.fromkeys(samples["group"].tolist()))
    sample_to_group = dict(zip(samples["sample"], samples["group"]))

    colors = pd.read_csv(color_table, sep=r"\s+", engine="python")
    group_colors = dict(zip(colors["group"].astype(str), colors["color"].astype(str)))

    all_data = []
    burden = []
    for sample in sample_order:
        path = cand_dir / f"{sample}.ati_candidates.tsv"
        df = pd.read_csv(path, sep="\t")
        df["sample"] = sample
        df["group"] = sample_to_group[sample]
        all_data.append(df)
        burden.append({"sample": sample, "group": sample_to_group[sample], "n_A_to_I_like_candidates": len(df)})

    all_df = pd.concat(all_data, ignore_index=True)
    burden_df = pd.DataFrame(burden)
    all_df.to_csv(summary_dir / "all_A_to_I_like_candidates.tsv", sep="\t", index=False)
    burden_df.to_csv(summary_dir / "per_sample_candidate_burden.tsv", sep="\t", index=False)

    max_n = max((sum(burden_df["group"] == g) for g in group_order), default=0)
    wide = {}
    for g in group_order:
        vals = burden_df.loc[burden_df["group"] == g, "n_A_to_I_like_candidates"].tolist()
        wide[g] = vals + [np.nan] * (max_n - len(vals))
    pd.DataFrame(wide).to_csv(summary_dir / "per_sample_candidate_burden_GraphPad_wide.tsv", sep="\t", index=False)

    group_summary = burden_df.groupby("group")["n_A_to_I_like_candidates"].agg(["count", "mean", "median", "std", "min", "max"]).reset_index()
    group_summary.to_csv(summary_dir / "group_candidate_burden_summary.tsv", sep="\t", index=False)

    if len(group_order) == 2:
        x = burden_df.loc[burden_df["group"] == group_order[0], "n_A_to_I_like_candidates"].to_numpy()
        y = burden_df.loc[burden_df["group"] == group_order[1], "n_A_to_I_like_candidates"].to_numpy()
        if len(x) >= 2 and len(y) >= 2:
            stat = mannwhitneyu(x, y, alternative="two-sided", method="exact")
            pd.DataFrame([{
                "comparison": f"{group_order[0]}_vs_{group_order[1]}",
                "test": "two-sided exact Mann-Whitney U test",
                "U_statistic": stat.statistic,
                "p_value": stat.pvalue,
            }]).to_csv(summary_dir / "group_comparison.tsv", sep="\t", index=False)

    x_pos = np.arange(len(sample_order))
    bar_values = burden_df.set_index("sample").loc[sample_order, "n_A_to_I_like_candidates"].to_numpy()
    bar_colors = [group_colors[sample_to_group[s]] for s in sample_order]
    fig = plt.figure(figsize=(max(5.2, 0.55 * len(sample_order)), 4.8))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 3.1], hspace=0.05)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ax1.bar(x_pos, bar_values, color=bar_colors, width=0.78, edgecolor="none")
    ax1.set_ylabel("Number of\nA-to-I-like sites", fontsize=9)
    ax1.tick_params(axis="x", labelbottom=False)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["bottom"].set_visible(False)

    rng = np.random.default_rng(1)
    for i, sample in enumerate(sample_order):
        values = all_df.loc[all_df["sample"] == sample, "frequency"].astype(float).to_numpy() * 100
        jitter = rng.uniform(-0.34, 0.34, len(values))
        ax2.scatter(np.full(len(values), i) + jitter, values, s=1.2, color=group_colors[sample_to_group[sample]], alpha=0.25, linewidths=0)

    ax2.set_ylabel("RNA A-to-I-like editing (%)", fontsize=9)
    ax2.set_ylim(0, 100)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(sample_order, rotation=45, ha="right", fontsize=8)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    f = cfg["filters"]
    ax2.text(0.01, 0.98, f"Filter: coverage >={f['min_coverage']}, alt reads >={f['min_alt_reads']}, alternative-base fraction >={f['min_alt_fraction']}", transform=ax2.transAxes, ha="left", va="top", fontsize=8)

    for ext in ["pdf", "svg", "png"]:
        out = fig_dir / f"per_sample_A_to_I_like_candidates.{ext}"
        fig.savefig(out, dpi=600 if ext == "png" else None, transparent=True, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
