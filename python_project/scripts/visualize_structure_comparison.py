import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from sklearn.datasets import make_blobs
from sklearn.neighbors import kneighbors_graph
from scipy.sparse import coo_matrix

ALT_FILE = Path(__file__).resolve().parent.parent / "data" / "sdss_100k_galaxy_form_burst.csv"
LOCAL_FILE = Path("data/galaxy_subset.csv")
OUT_IMG = Path("outputs/structural_comparison_3d.png")


def toy_data(n=3000, centers=3):
    coords, _ = make_blobs(n_samples=n, centers=centers, cluster_std=1.5, n_features=3, random_state=42)
    masses = np.ones(len(coords))
    return coords, masses


def realistic_toy(coords, masses, z_noise=5.0, mass_sigma=0.5):
    coords_noisy = coords.copy()
    coords_noisy[:, 2] += np.random.normal(0, z_noise, size=len(coords_noisy))
    masses_noisy = masses * np.random.lognormal(mean=0.0, sigma=mass_sigma, size=len(masses))
    return coords_noisy, masses_noisy


def load_sdss(limit=20000):
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
    print("SDSS data not found")
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
    return coords, mass


def centrality_vals(coords, masses, k=20):
    dist_graph = kneighbors_graph(coords, n_neighbors=k, mode="distance", include_self=False)
    coo = dist_graph.tocoo()
    w_data = masses[coo.row] * masses[coo.col] / (coo.data + 1e-8)
    weights = coo_matrix((w_data, (coo.row, coo.col)), shape=dist_graph.shape)
    G = nx.from_scipy_sparse_array(weights)
    # use pagerank as stable centrality
    pr = nx.pagerank(G, weight="weight", alpha=0.85, max_iter=200)
    cent = np.array([pr[i] for i in range(len(coords))])
    return cent


def main():
    out_dir = Path("outputs")
    out_dir.mkdir(exist_ok=True)

    # Ideal Toy
    coords_toy, masses_toy = toy_data()
    cent_toy = centrality_vals(coords_toy, masses_toy, k=20)

    # Realistic Toy
    coords_toy_rsd, masses_toy_rsd = realistic_toy(coords_toy, masses_toy, z_noise=5.0, mass_sigma=0.5)
    cent_toy_rsd = centrality_vals(coords_toy_rsd, masses_toy_rsd, k=20)

    # SDSS real slice
    df_sdss = load_sdss()
    coords_sdss, masses_sdss = preprocess_sdss(df_sdss)
    # sub-sample for display
    if len(coords_sdss) > 2000:
        coords_sdss = coords_sdss[:2000]
        masses_sdss = masses_sdss[:2000]
    cent_sdss = centrality_vals(coords_sdss, masses_sdss, k=10)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    data_list = [
        (coords_toy, cent_toy, "Théorie idéale"),
        (coords_toy_rsd, cent_toy_rsd, "Théorie réaliste (RSD + bruit)") ,
        (coords_sdss, cent_sdss, "Observation SDSS"),
    ]
    sc = None
    for ax, (coords, cent, title) in zip(axes, data_list):
        sc = ax.scatter(coords[:, 0], coords[:, 2], c=cent, cmap="plasma", s=2, alpha=0.6)
        ax.set_title(title)
        ax.set_xlabel("X")
        ax.set_ylabel("Z (ligne de visée)")
    fig.colorbar(sc, ax=axes, label="Centralité (PageRank)")
    plt.tight_layout()
    out_path = Path("outputs/structural_comparison_2d_XZ.png")
    plt.savefig(out_path, dpi=200)
    print(f"Image enregistrée : {out_path}")


if __name__ == "__main__":
    main()
