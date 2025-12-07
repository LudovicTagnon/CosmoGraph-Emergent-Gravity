import json
import os

import matplotlib.pyplot as plt


def main():
    stats_path = os.path.join(os.path.dirname(__file__), "..", "outputs", "sensitivity_stats.json")
    if not os.path.exists(stats_path):
        print(f"No stats file found at {stats_path}. Run sensitivity_analysis.py first.")
        return

    with open(stats_path, "r") as f:
        data = json.load(f)

    plt.figure(figsize=(8, 5))
    for cname, rows in data.items():
        ks = [r["k"] for r in rows]
        means = [r["corr_mean"] for r in rows]
        stds = [r["corr_std"] for r in rows]
        plt.errorbar(ks, means, yerr=stds, marker="o", capsize=3, label=cname)

    plt.xlabel("k (nearest neighbors)")
    plt.ylabel("corr(strength, centrality)")
    plt.title("Robustesse : corr vs k (bootstrap)")
    plt.legend()
    out_path = os.path.join(os.path.dirname(__file__), "..", "outputs", "sensitivity_corr_vs_k.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
