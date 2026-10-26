from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PROJECT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_DIR / "group_comparison"
FEATURES = ["Relative theta power", "Dominant frequency 4-13 Hz"]


def main():
    data = pd.read_csv(PROJECT_DIR / "EEG_features.csv")
    OUTPUT_DIR.mkdir(exist_ok=True)
    results, summaries = [], []
    group_counts = data["Group"].value_counts()

    for feature in FEATURES:
        pd_values = data.loc[data["Group"] == "PD", feature].to_numpy()
        control_values = data.loc[data["Group"] == "Control", feature].to_numpy()

        pd_mean = pd_values.mean()
        control_mean = control_values.mean()

        for group, values, mean in [("PD", pd_values, pd_mean),
                                    ("Control", control_values, control_mean)]:
            summaries.append({
                "Feature": feature, "Group": group, "N": len(values),
                "Mean": mean, "SD": values.std(ddof=1),
                "Median": np.median(values),
                "Q1": np.quantile(values, 0.25), "Q3": np.quantile(values, 0.75),
            })

        # Welch compares means without assuming equal group variances.
        test = stats.ttest_ind(pd_values, control_values, equal_var=False,
                               alternative="two-sided", nan_policy="raise")
        interval = test.confidence_interval(confidence_level=0.95)
        difference = pd_mean - control_mean

        # Hedges' g expresses the difference in pooled standard-deviation units.
        # The correction reduces the small-sample bias of Cohen's d.
        n_pd, n_control = len(pd_values), len(control_values)
        pooled_df = n_pd + n_control - 2
        pooled_sd = np.sqrt(
            ((n_pd - 1) * pd_values.var(ddof=1)
             + (n_control - 1) * control_values.var(ddof=1)) / pooled_df
        )
        hedges_g = (1 - 3 / (4 * pooled_df - 1)) * difference / pooled_sd
        results.append({
            "Feature": feature, "PD mean": pd_mean,
            "Control mean": control_mean,
            "Difference PD-Control": difference,
            "95% CI low (unadjusted)": interval.low,
            "95% CI high (unadjusted)": interval.high,
            "Hedges g (pooled SD)": hedges_g,
            "Welch t": test.statistic, "Welch df": test.df,
            "P raw": test.pvalue,
            "P Bonferroni": min(test.pvalue * len(FEATURES), 1.0),
        })

    summary = pd.DataFrame(summaries)
    results = pd.DataFrame(results)
    summary.to_csv(OUTPUT_DIR / "group_summary.csv", index=False)
    results.to_csv(OUTPUT_DIR / "group_tests.csv", index=False)

    # Dots are participants; boxes show the middle 50%; diamonds show means.
    rng = np.random.default_rng(42)
    colors = ["#24748C", "#C5653C"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    for ax, feature, title, ylabel, scale in zip(
        axes, FEATURES, ["Relative theta power", "Dominant frequency"],
        ["Theta power / 1–30 Hz power (%)", "Peak frequency in 4–13 Hz (Hz)"],
        [100, 1],
    ):
        for position, (group, color) in enumerate(zip(["Control", "PD"], colors)):
            values = data.loc[data["Group"] == group, feature].to_numpy() * scale
            ax.boxplot(values, positions=[position], widths=0.42,
                       patch_artist=True, showfliers=False,
                       boxprops={"facecolor": color, "alpha": 0.18},
                       medianprops={"color": color, "linewidth": 2})
            ax.scatter(position + rng.uniform(-0.16, 0.16, len(values)), values,
                       s=19, color=color, alpha=0.55, edgecolors="none", zorder=3)
            ax.scatter(position, values.mean(), marker="D", color="black",
                       s=42, zorder=4)
        corrected = results.loc[results["Feature"] == feature, "P Bonferroni"].iloc[0]
        ax.set_title(f"{title}\nBonferroni-adjusted p = {corrected:.3g}", fontsize=12)
        ax.set_xticks([0, 1], [f"Control (n={group_counts['Control']})",
                               f"PD (n={group_counts['PD']})"])
        ax.set_ylabel(ylabel)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
    fig.suptitle("EEG slowing features: PD and controls", fontsize=16)
    fig.text(0.5, 0.025,
             "Each dot = one participant · Box = middle 50% · Line = median · Black diamond = mean",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    fig.savefig(OUTPUT_DIR / "slowing_comparison.png", dpi=200)
    fig.savefig(OUTPUT_DIR / "slowing_comparison.svg")
    plt.close(fig)

    boundaries = pd.crosstab(data["Group"], data["Peak at search boundary"])
    print(summary.to_string(index=False))
    print(results.to_string(index=False))
    print("\nPeak boundary counts:\n", boundaries)
    print(f"\nOutput: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
