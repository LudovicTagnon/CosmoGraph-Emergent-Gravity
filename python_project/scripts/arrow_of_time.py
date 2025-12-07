from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path


def von_neumann_entropy(G: nx.Graph) -> float:
    if G.number_of_nodes() == 0:
        return 0.0
    L = nx.laplacian_matrix(G).astype(float)
    trace = L.diagonal().sum()
    if trace <= 0:
        return 0.0
    rho = L / trace
    vals = np.linalg.eigvalsh(rho.toarray())
    vals = vals[vals > 1e-12]
    return float(-(vals * np.log(vals)).sum()) if len(vals) else 0.0


def run_arrow_of_time(steps: int = 100, m_attach: int = 2):
    # Big Bang initial
    G = nx.barabasi_albert_graph(10, max(1, m_attach - 1), seed=42)
    entropies = []
    times = []
    for t in range(steps):
        entropies.append(von_neumann_entropy(G))
        times.append(t)

        # Expansion : ajouter un nœud, attachement préférentiel simple
        new_idx = G.number_of_nodes()
        targets = list(G.nodes())
        degrees = np.array([deg for _, deg in G.degree()], dtype=float)
        probs = degrees / degrees.sum() if degrees.sum() > 0 else np.ones_like(degrees) / len(degrees)
        attach = np.random.choice(targets, size=max(1, m_attach), replace=False, p=probs)
        for tgt in attach:
            G.add_edge(new_idx, tgt)

    entropies = np.array(entropies)
    times = np.array(times)
    dS = np.diff(entropies)
    monotonic = np.all(dS >= -1e-6)

    return times, entropies, dS, monotonic


def main():
    outputs = Path(__file__).resolve().parents[1] / "outputs"
    outputs.mkdir(parents=True, exist_ok=True)
    img_path = outputs / "arrow_of_time.png"
    verdict_path = outputs / "thermo_law.txt"

    times, entropies, dS, monotonic = run_arrow_of_time()

    plt.figure(figsize=(10, 6))
    plt.plot(times, entropies, color="orange", label="S_VN(t)")
    plt.xlabel("Temps cosmologique (itérations)")
    plt.ylabel("Entropie de Von Neumann")
    plt.title("Flèche du Temps : Evolution de l'Entropie")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(img_path)

    verdict = "Respecté" if monotonic else "Violé"
    verdict_text = (
        f"Flèche du Temps : {verdict}\n"
        f"Entropie finale : {entropies[-1]:.4f}\n"
        f"dS/dt min       : {dS.min():.6f}\n"
    )
    verdict_path.write_text(verdict_text, encoding="utf-8")
    print(verdict_text.strip())
    print(f"Graphique : {img_path}")


if __name__ == "__main__":
    main()
