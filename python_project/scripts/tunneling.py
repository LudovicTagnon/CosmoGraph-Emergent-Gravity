import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from matplotlib.colors import LogNorm
from pathlib import Path

# Parameters
Lx, Ly = 100, 40
wall_start, wall_end = 40, 45  # wall spans x in [40,45)
V0 = 2.0
x0, y0 = 20.0, 20.0
sigma = 2.0
kx = 1.5
k_energy = kx ** 2 / 2  # ~1.125 < V0
t_final = 40.0


def build_potential():
    V = np.zeros((Lx, Ly), dtype=float)
    V[wall_start:wall_end, :] = V0
    return V


def build_graph():
    G = nx.grid_2d_graph(Lx, Ly)
    mapping = {node: i for i, node in enumerate(G.nodes())}
    G = nx.relabel_nodes(G, mapping)
    pos = {mapping[node]: node for node in nx.grid_2d_graph(Lx, Ly).nodes()}
    return G, pos


def hamiltonian(G, V, pos):
    L = nx.laplacian_matrix(G).astype(complex)  # sparse
    diag_V = diags([V[x, y] for _, (x, y) in sorted(pos.items())], offsets=0, dtype=complex, shape=L.shape)
    return L + diag_V


def initial_packet(pos):
    psi0 = np.zeros(len(pos), dtype=complex)
    for idx, (x, y) in pos.items():
        amp = np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * sigma ** 2))
        phase = np.exp(1j * kx * x)
        psi0[idx] = amp * phase
    norm = np.linalg.norm(psi0)
    if norm > 0:
        psi0 /= norm
    return psi0


def main():
    out_dir = Path("outputs")
    out_dir.mkdir(exist_ok=True)

    V = build_potential()
    G, pos = build_graph()
    H = hamiltonian(G, V, pos)
    psi0 = initial_packet(pos)

    psi_t = expm_multiply(-1j * H * t_final, psi0)
    prob = np.abs(psi_t) ** 2

    grid = np.zeros((Lx, Ly))
    for idx, (x, y) in pos.items():
        grid[x, y] = prob[idx]

    prob_right = sum(prob[idx] for idx, (x, y) in pos.items() if x > wall_end)

    plt.figure(figsize=(12, 5))
    plt.imshow(
        grid.T,
        origin="lower",
        cmap="inferno",
        extent=[0, Lx, 0, Ly],
        norm=LogNorm(vmin=1e-6, vmax=1.0),
    )
    plt.axvspan(wall_start, wall_end, color="white", alpha=0.3, label=f"Barrière V={V0}")
    plt.title("Effet tunnel : paquet face à une barrière (E < V0)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.colorbar(label="Probabilité", fraction=0.046, pad=0.04)
    plt.legend(loc="upper right")
    plt.tight_layout()
    out_path = out_dir / "tunneling.png"
    plt.savefig(out_path, dpi=200)

    print(f"Image enregistrée : {out_path}")
    print(f"Énergie cinétique ~ {k_energy:.3f} (doit être < V0={V0})")
    print(f"Probabilité à droite du mur (x>{wall_end}) : {prob_right:.3e}")


if __name__ == "__main__":
    main()
