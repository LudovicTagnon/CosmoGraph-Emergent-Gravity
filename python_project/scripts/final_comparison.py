import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.stats import pearsonr
from sklearn.datasets import make_blobs
from sklearn.neighbors import kneighbors_graph

ALT_FILE = Path(__file__).resolve().parent.parent / "data" / "sdss_100k_galaxy_form_burst.csv"
LOCAL_FILE = Path("data/galaxy_subset.csv")
OUT_IMG = Path("outputs/final_theory_vs_reality.png")


def build_weighted_graph(coords, masses, k):
    dist_graph = kneighbors_graph(coords, n_neighbors=k, mode="distance", include_self=False)
    coo = dist_graph.tocoo()
    w_data = masses[coo.row] * masses[coo.col] / (coo.data + 1e-8)
    weights = coo_matrix((w_data, (coo.row, coo.col)), shape=dist_graph.shape)
    G_full = nx.from_scipy_sparse_array(weights)
    # largest component
    comps = sorted(nx.connected_components(G_full), key=len, reverse=True)
    G = G_full.subgraph(comps[0]).copy()
    lcc_nodes = np.array(sorted(G.nodes()))
    return G, lcc_nodes


def corr_strength_centrality(coords, masses, k=20, margin=0.1):
    G, lcc_nodes = build_weighted_graph(coords, masses, k)
    strength = np.array([d for _, d in G.degree(weight="weight")], dtype=float)
    try:
        cent_dict = nx.eigenvector_centrality(G, max_iter=2000, tol=1e-06, weight="weight")
        cent = np.array([cent_dict[i] for i in lcc_nodes])
    except Exception:
        vals, vecs = np.linalg.eig(nx.to_numpy_array(G))
        idx = np.argmax(vals.real)
        cent = np.abs(vecs[:, idx].real)
        cent /= cent.max()
    # inner mask (buffer) for fairness
    if margin and margin > 0:
        pts = coords[lcc_nodes]
        mins = pts.min(axis=0)
        maxs = pts.max(axis=0)
        rng = maxs - mins
        lo = mins + margin * rng
        hi = maxs - margin * rng
        mask = np.all((pts > lo) & (pts < hi), axis=1)
        strength = strength[mask]
        cent = cent[mask]
    corr, _ = pearsonr(strength, cent)
    return corr


def fetch_sdss(limit=20000):
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
                return df
            except Exception:
                continue
    print("SDSS data not found. Aborting.")
    sys.exit(1)


def preprocess_sdss(df):
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


def metrics_sdss(coords, masses, df_slice, k=10):
    G, lcc_nodes = build_weighted_graph(coords, masses, k)
    strength = np.array([d for _, d in G.degree(weight="weight")], dtype=float)
    try:
        cent_dict = nx.eigenvector_centrality(G, max_iter=2000, tol=1e-06, weight="weight")
        cent = np.array([cent_dict[i] for i in lcc_nodes])
    except Exception:
        vals, vecs = np.linalg.eig(nx.to_numpy_array(G))
        idx = np.argmax(vals.real)
        cent = np.abs(vecs[:, idx].real)
        cent /= cent.max()
    mask_inner = (
        (df_slice["ra"] > 135)
        & (df_slice["ra"] < 235)
        & (df_slice["dec"] > 0)
        & (df_slice["dec"] < 55)
        & (df_slice["z"] > 0.05)
        & (df_slice["z"] < 0.11)
    ).values
    mask_lcc = mask_inner[lcc_nodes]
    strength_f = strength[mask_lcc]
    cent_f = cent[mask_lcc]
    corr, _ = pearsonr(strength_f, cent_f)
    return corr, len(strength_f)


def main():
    out_dir = Path("outputs")
    out_dir.mkdir(exist_ok=True)

    # Ideal theory
    pts, _ = make_blobs(n_samples=3000, centers=3, cluster_std=1.5, n_features=3, random_state=42)
    masses = np.ones(len(pts))
    corr_ideal = corr_strength_centrality(pts, masses, k=20, margin=0.1)

    # Realistic theory: add RSD on z and mass scatter
    pts_real = pts.copy()
    pts_real[:, 2] += np.random.normal(0, 5.0, size=len(pts_real))
    masses_real = masses * np.random.lognormal(mean=0.0, sigma=0.5, size=len(masses))
    corr_realistic = corr_strength_centrality(pts_real, masses_real, k=20, margin=0.1)

    # Reality SDSS
    df = fetch_sdss()
    coords, masses_sdss, df_slice = preprocess_sdss(df)
    corr_sdss, n_sdss = metrics_sdss(coords, masses_sdss, df_slice, k=10)

    labels = ["Théorie idéale", "Théorie réaliste (RSD+bruit)", "Observation SDSS"]
    values = [corr_ideal, corr_realistic, corr_sdss]

    plt.figure(figsize=(7, 5))
    bars = plt.bar(labels, values, color=["tab:blue", "tab:purple", "tab:orange"], alpha=0.8)
    for bar, val in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, val + 0.01, f"{val:.3f}", ha="center", va="bottom")
    plt.ylim(0, max(values) + 0.1)
    plt.ylabel("Corrélation (force vs centralité)")
    plt.title("Théorie vs Réalité (avec RSD et bruit de masse)")
    plt.tight_layout()
    out_dir.mkdir(exist_ok=True)
    plt.savefig(OUT_IMG, dpi=200)

    print(f"Corr Théorie idéale: {corr_ideal:.3f}")
    print(f"Corr Théorie réaliste: {corr_realistic:.3f}")
    print(f"Corr Réalité SDSS (k=10, inner mask, n={n_sdss}): {corr_sdss:.3f}")
    print(f"Figure enregistrée : {OUT_IMG}")


if __name__ == "__main__":
    main()
