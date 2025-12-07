from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path


def calculate_sigma(G: nx.Graph) -> float:
    if len(G) < 10 or not nx.is_connected(G):
        return 0.0
    try:
        C = nx.average_clustering(G)
        L = nx.average_shortest_path_length(G)
        n = len(G)
        k = np.mean([d for _, d in G.degree()])
        if k <= 1:
            return 0.0
        Cr = k / n
        Lr = np.log(n) / np.log(k) if k > 1 else n
        if Cr == 0 or L == 0:
            return 0.0
        sigma = (C / Cr) / (L / Lr)
        return sigma
    except Exception:
        return 0.0


def scan_multiverse():
    results = []
    m_values = [1, 2, 3, 4, 5, 6]
    sizes = [50, 100, 200]

    print("🌌 Lancement du Scan Multivers (Sigma Metric)...")
    for n in sizes:
        for m in m_values:
            if m >= n:
                continue
            G = nx.barabasi_albert_graph(n, m)
            sigma = calculate_sigma(G)
            density = nx.density(G)
            diameter = nx.diameter(G) if nx.is_connected(G) else 0

            status = "MORT"
            if sigma > 1.0:
                status = "VIABLE"
            if sigma > 2.5:
                status = "OPTIMAL"

            results.append(
                {
                    "Nodes": n,
                    "Connectivity_m": m,
                    "Sigma_Complexity": round(sigma, 4),
                    "Density": round(density, 4),
                    "Diameter": diameter,
                    "Status": status,
                }
            )
            print(f"   Univers (N={n}, m={m}) -> Sigma={sigma:.2f} [{status}]")

    df = pd.DataFrame(results)
    outputs = Path(__file__).resolve().parents[1] / "outputs"
    outputs.mkdir(exist_ok=True, parents=True)
    csv_path = outputs / "multiverse_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"✅ Données sauvegardées : {csv_path}")

    max_n = max(sizes)
    subset = df[df["Nodes"] == max_n]
    if not subset.empty:
        plt.figure(figsize=(10, 6))
        sns.barplot(data=subset, x="Connectivity_m", y="Sigma_Complexity", palette="viridis")
        plt.axhline(1.0, color="red", linestyle="--", label="Seuil de Viabilité (Sigma=1)")
        plt.title(f"Paysage Anthropique : Complexité selon la Connectivité (N={max_n})")
        plt.ylabel("Complexité Sigma (Small-Worldness)")
        plt.xlabel("Paramètre de Connectivité (m)")
        plt.legend()
        plt.grid(axis="y", alpha=0.3)
        img_path = outputs / "multiverse_landscape.png"
        plt.tight_layout()
        plt.savefig(img_path)
        print(f"✅ Visuel généré : {img_path}")


if __name__ == "__main__":
    scan_multiverse()
