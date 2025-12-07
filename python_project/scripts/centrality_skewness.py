import os
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.neighbors import kneighbors_graph
from scipy.stats import skew
import matplotlib.pyplot as plt


def load_sdss(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", low_memory=False)
    if {"ra", "dec", "z"}.issubset(df.columns):
        pass
    elif {"ra", "dec", "redshift"}.issubset(df.columns):
        df = df.rename(columns={"redshift": "z"})
    else:
        raise ValueError("CSV must contain ra, dec and either z or redshift")
    return df


def to_cartesian(df: pd.DataFrame) -> np.ndarray:
    ra = np.deg2rad(df["ra"].to_numpy())
    dec = np.deg2rad(df["dec"].to_numpy())
    r = df["z"].to_numpy()
    x = r * np.cos(dec) * np.cos(ra)
    y = r * np.cos(dec) * np.sin(ra)
    zc = r * np.sin(dec)
    return np.column_stack([x, y, zc])


def main():
    data_path = os.path.join(os.path.dirname(__file__), "..", "data", "sdss_100k_galaxy_form_burst.csv")
    df = load_sdss(data_path)
    df = df.head(5000)  # keep manageable for fast viz
    coords = to_cartesian(df)

    conn = kneighbors_graph(coords, n_neighbors=10, mode="distance", include_self=False)
    sources, targets = conn.nonzero()
    weights = conn.data
    G = nx.Graph()
    G.add_nodes_from(range(coords.shape[0]))
    for s, t, w in zip(sources, targets, weights):
        if s == t:
            continue
        G.add_edge(int(s), int(t), weight=float(1.0 / (w + 1e-9)))
    if not nx.is_connected(G):
        largest = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest).copy()

    central = np.array(list(nx.eigenvector_centrality(G, max_iter=1000, weight="weight").values()))
    strength = np.array([w for _, w in G.degree(weight="weight")], float)

    stats = {
        "N_LCC": len(G),
        "corr_strength_central": float(np.corrcoef(strength, central)[0, 1]),
        "central_skew": float(skew(central)),
        "central_min": float(central.min()),
        "central_median": float(np.median(central)),
        "central_p95": float(np.percentile(central, 95)),
        "central_max": float(central.max()),
        "strength_min": float(strength.min()),
        "strength_median": float(np.median(strength)),
        "strength_p95": float(np.percentile(strength, 95)),
        "strength_max": float(strength.max()),
    }
    print("Stats:", stats)

    # Histogram centrality
    plt.figure(figsize=(7, 4))
    plt.hist(central, bins=100, log=True, color="steelblue", alpha=0.8)
    plt.xlabel("Eigenvector centrality")
    plt.ylabel("Count (log)")
    plt.title(f"Centrality skewness={stats['central_skew']:.1f} (N={stats['N_LCC']})")
    out_path = os.path.join(os.path.dirname(__file__), "..", "outputs", "centrality_hist.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
