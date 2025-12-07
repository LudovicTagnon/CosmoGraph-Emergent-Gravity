import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import expm_multiply
from matplotlib.colors import LogNorm

# Parameters
Lx, Ly = 80, 40
wall_x = 25
wall_thickness = 2
V_wall = 100.0
slit_width = 3
slit_y1 = 14
slit_y2 = 26
sigma = 4.0
kx = 3.0
t_final = 30.0
vmax_plot = 0.01


def build_potential():
    V = np.zeros((Lx, Ly), dtype=float)
    # wall region
    V[wall_x : wall_x + wall_thickness, :] = V_wall
    # carve slits (set back to 0)
    def apply_slit(center_y):
        half = slit_width // 2
        y_start = max(0, center_y - half)
        y_end = min(Ly, center_y + half + 1)
        V[wall_x : wall_x + wall_thickness, y_start:y_end] = 0.0
    apply_slit(slit_y1)
    apply_slit(slit_y2)
    return V


def build_graph():
    G = nx.grid_2d_graph(Lx, Ly)
    mapping = {node: i for i, node in enumerate(G.nodes())}
    G = nx.relabel_nodes(G, mapping)
    pos = {mapping[node]: node for node in nx.grid_2d_graph(Lx, Ly).nodes()}
    return G, pos


def hamiltonian_sparse(G, V, pos):
    L = nx.laplacian_matrix(G).astype(complex)  # sparse
    diag_V = diags([V[x, y] for _, (x, y) in sorted(pos.items())], offsets=0, dtype=complex, shape=L.shape)
    return L + diag_V


def initial_packet(pos):
    psi0 = np.zeros(len(pos), dtype=complex)
    x0, y0 = 10.0, 20.0
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
    H = hamiltonian_sparse(G, V, pos)
    psi0 = initial_packet(pos)

    psi_t = expm_multiply(-1j * H * t_final, psi0)
    prob = np.abs(psi_t) ** 2

    # Heatmap grid
    grid = np.zeros((Lx, Ly))
    for idx, (x, y) in pos.items():
        grid[x, y] = prob[idx]

    plt.figure(figsize=(10, 5))
    plt.imshow(
        grid.T,
        origin="lower",
        cmap="inferno",
        extent=[0, Lx, 0, Ly],
        norm=LogNorm(vmin=1e-6, vmax=1.0),
    )
    plt.axvspan(wall_x, wall_x + wall_thickness, color="white", alpha=0.2, label="Mur V=100")
    plt.title("Double fente (Young) | t={:.1f}".format(t_final))
    plt.xlabel("x")
    plt.ylabel("y")
    plt.colorbar(label="Probabilité", fraction=0.046, pad=0.04)
    plt.legend(loc="upper right")
    plt.tight_layout()
    out_path = out_dir / "double_slit_fixed.png"
    plt.savefig(out_path, dpi=200)

    # screen profile at rightmost column
    x_max = Lx - 1
    screen_vals = []
    screen_y = []
    for idx, (x, y) in pos.items():
        if x == x_max:
            screen_vals.append(prob[idx])
            screen_y.append(y)
    if screen_vals:
        screen_y, screen_vals = zip(*sorted(zip(screen_y, screen_vals)))
        screen_vals = np.array(screen_vals)
        screen_total = screen_vals.sum()
        norm_vals = screen_vals / screen_vals.max() if screen_vals.max() > 0 else screen_vals
        mean_x = sum(prob[idx] * pos[idx][0] for idx in range(len(prob)))
        # print metrics
        print(f"Probabilité totale sur l'écran (x={x_max}) : {screen_total:.3e}")
        print(f"<x> moyen ≈ {mean_x:.2f} / {Lx}")
        # find peaks
        peaks = []
        for i in range(1, len(norm_vals) - 1):
            if norm_vals[i] > norm_vals[i - 1] and norm_vals[i] > norm_vals[i + 1]:
                peaks.append((screen_y[i], norm_vals[i]))
        peaks_sorted = sorted(peaks, key=lambda p: p[1], reverse=True)
        print(f"Nb de pics (franges) détectés : {len(peaks_sorted)}")
        print("Top 5 pics (y, amplitude norm.):", peaks_sorted[:5])
    print(f"Heatmap: {out_path}")


if __name__ == "__main__":
    main()
