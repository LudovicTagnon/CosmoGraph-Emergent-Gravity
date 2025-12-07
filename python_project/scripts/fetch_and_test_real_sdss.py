import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.stats import pearsonr
from sklearn.neighbors import kneighbors_graph

URL = "https://raw.githubusercontent.com/astroML/astroML-data/master/datasets/SDSS_specgals.csv"
ALT_FILE = Path(__file__).resolve().parent.parent / "data" / "sdss_100k_galaxy_form_burst.csv"
LOCAL_FILE = Path("data/galaxy_subset.csv")
OUT_IMG = Path("outputs/real_sdss_result.png")
SENS_IMG = Path("outputs/sensitivity_k.png")


def fetch_data(limit=20000):
    for path in [ALT_FILE, LOCAL_FILE]:
        if path.exists():
            try:
                df_all = pd.read_csv(path, comment="#")
                red_col = "redshift" if "redshift" in df_all.columns else "z" if "z" in df_all.columns else None
                if red_col is None or not {"ra", "dec"}.issubset(df_all.columns):
                    continue
                cols = ["ra", "dec", red_col]
                if "petroMag_r" in df_all.columns:
                    cols.append("petroMag_r")
                df = df_all[cols].rename(columns={red_col: "z"})
                if limit:
                    df = df.head(limit)
                return df, False
            except Exception:
                continue
    try:
        df = pd.read_csv(URL, usecols=["ra", "dec", "z"], nrows=limit)
        return df, True
    except Exception as e:
        print("Download failed:", e)
        print("Provide data/galaxy_subset.csv with ra, dec, z")
        sys.exit(1)


def preprocess(df):
    df = df[
        (df["z"] > 0.04)
        & (df["z"] < 0.12)
        & (df["ra"] > 130)
        & (df["ra"] < 240)
        & (df["dec"] > -5)
        & (df["dec"] < 60)
    ].copy()
    if len(df) > 20000:
        df = df.sort_values(by="ra").head(20000)
    if "petroMag_r" in df.columns:
        mass = 10 ** (-0.4 * df["petroMag_r"].values)
        mass = mass / (mass.mean() + 1e-12)
    else:
        mass = np.ones(len(df))
    ra = np.deg2rad(df["ra"].values)
    dec = np.deg2rad(df["dec"].values)
    z = df["z"].values
    x = z * np.cos(dec) * np.cos(ra)
    y = z * np.cos(dec) * np.sin(ra)
    zc = z * np.sin(dec)
    coords = np.stack([x, y, zc], axis=1)
    return coords, mass, df.reset_index(drop=True)


def metrics_on_coords(coords, masses, k=20):
    dist_graph = kneighbors_graph(coords, n_neighbors=k, mode="distance", include_self=False)
    coo = dist_graph.tocoo()
    w_data = masses[coo.row] * masses[coo.col] / (coo.data + 1e-8)
    weights = coo_matrix((w_data, (coo.row, coo.col)), shape=dist_graph.shape)
    G_full = nx.from_scipy_sparse_array(weights)
    comps = sorted(nx.connected_components(G_full), key=len, reverse=True)
    G = G_full.subgraph(comps[0]).copy()
    lcc_nodes = np.array(sorted(G.nodes()))
    coords_lcc = coords[lcc_nodes]
    strength = np.array([d for _, d in G.degree(weight="weight")], dtype=float)
    try:
        cent_dict = nx.eigenvector_centrality(G, max_iter=2000, tol=1e-06, weight="weight")
        cent = np.array([cent_dict[i] for i in lcc_nodes])
    except Exception:
        vals, vecs = np.linalg.eig(nx.to_numpy_array(G))
        idx = np.argmax(vals.real)
        cent = np.abs(vecs[:, idx].real)
        cent /= cent.max()
    return strength, cent, G, coords_lcc, lcc_nodes


def plot_web(coords, centrality, path):
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=centrality, cmap="plasma", s=5)
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.set_title("SDSS Web (color = centralité)")
    fig.colorbar(sc, ax=ax, label="Eigenvector centrality")
    plt.tight_layout()
    path.parent.mkdir(exist_ok=True)
    plt.savefig(path, dpi=200)
    plt.close(fig)


def main():
    out_dir = Path("outputs")
    out_dir.mkdir(exist_ok=True)

    df, downloaded = fetch_data(limit=20000)
    coords, masses, df_slice = preprocess(df)
    if len(coords) == 0:
        print("No galaxies after redshift filter; adjust z range.")
        sys.exit(1)

    k_values = [10, 20, 30, 50, 75, 100]
    results = []
    best = None
    # buffer mask on original slice
    mask_inner = (
        (df_slice["ra"] > 135)
        & (df_slice["ra"] < 235)
        & (df_slice["dec"] > 0)
        & (df_slice["dec"] < 55)
        & (df_slice["z"] > 0.05)
        & (df_slice["z"] < 0.11)
    ).values

    for k in k_values:
        strength, cent, G, coords_lcc, lcc_nodes = metrics_on_coords(coords, masses, k=k)
        mask_lcc = mask_inner[lcc_nodes]
        strength_f = strength[mask_lcc]
        cent_f = cent[mask_lcc]
        corr, _ = pearsonr(strength_f, cent_f)
        results.append((k, corr))
        if best is None or corr > best[1]:
            best = (k, corr, strength_f, cent_f, coords_lcc[mask_lcc], G, mask_lcc)
        print(f"k={k}: corr(strength, centrality)={corr:.3f} | nodes LCC={G.number_of_nodes()} | kept {mask_lcc.sum()}")

    ks = [r[0] for r in results]
    corrs = [r[1] for r in results]
    plt.figure(figsize=(6, 4))
    plt.plot(ks, corrs, marker="o")
    plt.xlabel("Nombre de voisins (k)")
    plt.ylabel("Corrélation force vs centralité")
    plt.title("Impact of Graph Connectivity on Physical Correlation")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(SENS_IMG, dpi=200)

    best_k, best_corr, strength_b, cent_b, coords_b, G_b, mask_b = best
    plot_web(coords_b, cent_b, OUT_IMG)

    print(f"Best k={best_k} corr={best_corr:.3f}")
    print(f"Plot saved: {OUT_IMG}")
    print(f"Sensitivity plot: {SENS_IMG}")


if __name__ == "__main__":
    main()
