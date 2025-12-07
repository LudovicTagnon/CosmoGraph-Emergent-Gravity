from __future__ import annotations

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def generate_quantum_network(n_nodes: int = 1000) -> nx.Graph:
    return nx.barabasi_albert_graph(n_nodes, 3, seed=42)


def find_dark_matter_candidates(G: nx.Graph):
    print("🕵️‍♂️ Analyse de la topologie du réseau...")
    degree = dict(G.degree())
    max_degree = max(degree.values()) if degree else 1
    try:
        centrality = nx.eigenvector_centrality(G, max_iter=1000)
    except Exception:
        centrality = nx.degree_centrality(G)

    vals_degree = np.array(list(degree.values()), dtype=float)
    vals_centrality = np.array(list(centrality.values()), dtype=float)

    norm_degree = vals_degree / max(vals_degree) if len(vals_degree) else np.array([])
    norm_centrality = vals_centrality / max(vals_centrality) if len(vals_centrality) else np.array([])
    discrepancy = norm_centrality - norm_degree

    threshold = np.percentile(discrepancy, 95) if len(discrepancy) else 0.0
    print(f"📊 Seuil de détection (95e percentile) : {threshold:.4f}")

    candidate_indices = np.where(discrepancy > threshold)[0]
    return candidate_indices, discrepancy


def main():
    print("🌌 Initialisation du scanner de Matière Noire...")
    G = generate_quantum_network(1000)

    candidates, scores = find_dark_matter_candidates(G)
    print(f"⚠️ DÉTECTION : {len(candidates)} amas de matière noire trouvés.")

    pos = nx.spring_layout(G, seed=42, k=0.15)
    out_path = Path(__file__).resolve().parents[1] / "outputs" / "dark_matter_map.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(12, 8))
    nx.draw_networkx_nodes(G, pos, node_size=10, node_color="#333333", alpha=0.3, label="Espace-Temps")

    dark_nodes = [list(G.nodes())[i] for i in candidates]
    if dark_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=dark_nodes, node_size=100, node_color="red", alpha=0.9, label="Matière Noire (Anomalies)")

    plt.title("Carte de la Matière Noire Émergente\n(Nœuds à haute influence / faible connectivité)", color="white")
    plt.legend()
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_path, facecolor="black")
    print(f"📸 Carte générée : {out_path}")


if __name__ == "__main__":
    main()
