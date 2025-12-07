from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx


def main():
    print("--- PHASE 5 : SIMULATION DE L'EXPANSION COSMIQUE (LOI DE HUBBLE) ---")

    initial_nodes = 20
    m_parameter = 2
    G = nx.barabasi_albert_graph(initial_nodes, m_parameter, seed=42)

    galaxy_A = 0
    galaxy_B = 1
    initial_dist = nx.shortest_path_length(G, galaxy_A, galaxy_B)

    print(f"1. Univers Primordial né ({initial_nodes} nœuds).")
    print(f"   Distance initiale entre Galaxie A et B : {initial_dist}")

    time_steps = []
    avg_distances = []
    galaxy_distances = []

    target_nodes = 300
    current_nodes = initial_nodes

    print("2. Démarrage de l'Inflation...")

    while current_nodes < target_nodes:
        current_nodes += 10
        G_expanded = nx.barabasi_albert_graph(current_nodes, m_parameter, seed=42)

        avg_dist = nx.average_shortest_path_length(G_expanded)
        try:
            g_dist = nx.shortest_path_length(G_expanded, galaxy_A, galaxy_B)
        except nx.NetworkXNoPath:
            g_dist = 0

        time_steps.append(current_nodes)
        avg_distances.append(avg_dist)
        galaxy_distances.append(g_dist)

    print(f"3. Inflation terminée. Taille finale : {current_nodes} nœuds.")
    print(f"   Distance finale moyenne : {avg_distances[-1]:.4f}")
    print(f"   Distance finale A-B : {galaxy_distances[-1]}")

    fig, ax1 = plt.subplots(figsize=(10, 6))

    color = "tab:blue"
    ax1.set_xlabel("Temps (Nombre de Nœuds)")
    ax1.set_ylabel("Distance Moyenne (Métrique)", color=color)
    ax1.plot(time_steps, avg_distances, color=color, linewidth=2, label="Expansion Moyenne")
    ax1.tick_params(axis="y", labelcolor=color)
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    color = "tab:red"
    ax2.set_ylabel("Distance Galaxie A-B (Sauts)", color=color)
    ax2.plot(time_steps, galaxy_distances, color=color, linestyle="--", label="Distance Galaxie A-B")
    ax2.tick_params(axis="y", labelcolor=color)

    plt.title("Loi de Hubble Algorithmique : Expansion de l'Espace-Temps")
    fig.tight_layout()

    out_path = (
        Path(__file__).resolve().parents[1] / "outputs" / "cosmic_expansion.png"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    print(f"Figure enregistrée : {out_path}")

    if avg_distances[-1] > avg_distances[0]:
        print(">>> CONCLUSION : L'univers est en expansion (le chemin moyen s'allonge).")
    else:
        print(">>> CONCLUSION : L'univers est stable ou se contracte (Small World dominant).")


if __name__ == "__main__":
    from pathlib import Path

    main()
