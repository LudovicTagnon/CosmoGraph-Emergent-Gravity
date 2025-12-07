from __future__ import annotations

from pathlib import Path
import math

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def shannon_entropy_degree(G: nx.Graph) -> float:
    degrees = [d for _, d in G.degree()]
    if not degrees:
        return 0.0
    values, counts = np.unique(degrees, return_counts=True)
    probs = counts / counts.sum()
    return float(-(probs * np.log(probs)).sum())


def von_neumann_entropy(G: nx.Graph) -> float:
    if G.number_of_nodes() == 0:
        return 0.0
    L = nx.laplacian_matrix(G).astype(float)
    trace = L.diagonal().sum()
    if trace <= 0:
        return 0.0
    rho = L / trace
    # eigvalsh for symmetric matrix
    vals = np.linalg.eigvalsh(rho.toarray())
    vals = vals[vals > 1e-12]
    return float(-(vals * np.log(vals)).sum()) if len(vals) else 0.0


def summary(name: str, G: nx.Graph) -> dict:
    degrees = [d for _, d in G.degree()]
    return {
        "name": name,
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "avg_degree": float(np.mean(degrees)) if degrees else 0.0,
        "max_degree": int(max(degrees)) if degrees else 0,
        "shannon_degree_entropy": shannon_entropy_degree(G),
        "von_neumann_entropy": von_neumann_entropy(G),
    }


def main():
    out_dir = Path(__file__).resolve().parents[1] / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "thermodynamics_report.txt"
    fig_path = out_dir / "entropy_evolution.png"

    print("--- PHASE 6 : ANALYSE THERMODYNAMIQUE ---")

    # Univers jeune (peu structuré)
    G_young = nx.erdos_renyi_graph(30, 0.1, seed=42)
    # Univers mûr (structuré avec hubs, expansion)
    G_mature = nx.barabasi_albert_graph(150, 3, seed=42)

    stats_y = summary("Univers jeune (ER 30, p=0.1)", G_young)
    stats_m = summary("Univers mûr (BA 150, m=3)", G_mature)

    # Rapport texte
    with report_path.open("w", encoding="utf-8") as f:
        for stats in (stats_y, stats_m):
            f.write(f"{stats['name']}\n")
            f.write(f"  Noeuds      : {stats['nodes']}\n")
            f.write(f"  Liens       : {stats['edges']}\n")
            f.write(f"  Degre moyen : {stats['avg_degree']:.3f}\n")
            f.write(f"  Degre max   : {stats['max_degree']}\n")
            f.write(f"  Entropie de Shannon (degres) : {stats['shannon_degree_entropy']:.6f}\n")
            f.write(f"  Entropie de Von Neumann      : {stats['von_neumann_entropy']:.6f}\n")
            f.write("\n")

        delta_shannon = stats_m["shannon_degree_entropy"] - stats_y["shannon_degree_entropy"]
        delta_vn = stats_m["von_neumann_entropy"] - stats_y["von_neumann_entropy"]
        f.write("Comparatif\n")
        f.write(f"  Δ Entropie Shannon (mûr - jeune) : {delta_shannon:.6f}\n")
        f.write(f"  Δ Entropie Von Neumann           : {delta_vn:.6f}\n")

    print(f"Rapport écrit : {report_path}")

    # Figure comparative
    labels = ["Shannon deg", "Von Neumann"]
    young_vals = [stats_y["shannon_degree_entropy"], stats_y["von_neumann_entropy"]]
    mature_vals = [stats_m["shannon_degree_entropy"], stats_m["von_neumann_entropy"]]

    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - width / 2, young_vals, width, label="Jeune")
    ax.bar(x + width / 2, mature_vals, width, label="Mûr")
    ax.set_ylabel("Entropie")
    ax.set_title("Entropie structurale et quantique\nUnivers jeune vs univers mûr")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(fig_path)
    print(f"Graphique généré : {fig_path}")


if __name__ == "__main__":
    main()
