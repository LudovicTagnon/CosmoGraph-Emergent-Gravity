from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from pathlib import Path


def measure_expansion():
    print("💥 Initialisation du Big Bang...")
    G = nx.barabasi_albert_graph(20, 2, seed=42)

    history_nodes = []
    history_avg_path = []

    epochs = 50
    nodes_per_epoch = 20

    print(f"🚀 Simulation de l'expansion sur {epochs} époques...")
    for t in range(epochs):
        current_nodes = list(G.nodes())
        new_nodes = range(len(current_nodes), len(current_nodes) + nodes_per_epoch)

        for new_node in new_nodes:
            target = np.random.choice(current_nodes)
            G.add_edge(new_node, target)

        if nx.is_connected(G):
            avg_path = nx.average_shortest_path_length(G)
        else:
            avg_path = 0

        history_nodes.append(len(G.nodes()))
        history_avg_path.append(avg_path)

        if t % 10 == 0:
            print(f"   Époque {t}: {len(G.nodes())} nœuds, Distance Moyenne={avg_path:.2f}")

    return history_nodes, history_avg_path


def main():
    nodes, distances = measure_expansion()

    plt.figure(figsize=(10, 6))
    plt.plot(nodes, distances, marker="o", color="purple", linestyle="-")
    plt.title("Loi de Hubble Émergente : Expansion du Réseau")
    plt.xlabel("Temps (Nombre de Nœuds / Volume)")
    plt.ylabel("Distance Moyenne (Métrique de l'Espace)")
    plt.grid(True, linestyle="--", alpha=0.5)

    out_path = Path(__file__).resolve().parents[1] / "outputs" / "expansion_hubble.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"🌌 Graphique d'expansion généré : {out_path}")


if __name__ == "__main__":
    main()
