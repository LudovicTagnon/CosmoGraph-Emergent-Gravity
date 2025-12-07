import json
import os
from typing import Dict, List

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.neighbors import kneighbors_graph
from tqdm import tqdm


def load_sdss_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"ra", "dec", "z"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain columns {required}")
    return df


def to_cartesian(df: pd.DataFrame) -> np.ndarray:
    # Low-z approximation: D ~ z (or cz/H0); we just scale by z for relative geometry
    ra = np.deg2rad(df["ra"].to_numpy())
    dec = np.deg2rad(df["dec"].to_numpy())
    r = df["z"].to_numpy()
    x = r * np.cos(dec) * np.cos(ra)
    y = r * np.cos(dec) * np.sin(ra)
    zc = r * np.sin(dec)
    return np.column_stack([x, y, zc])


def build_graph(coords: np.ndarray, k: int) -> nx.Graph:
    conn = kneighbors_graph(coords, n_neighbors=k, mode="distance", include_self=False)
    sources, targets = conn.nonzero()
    weights = conn.data
    G = nx.Graph()
    G.add_nodes_from(range(coords.shape[0]))
    for s, t, w in zip(sources, targets, weights):
        if s == t:
            continue
        G.add_edge(int(s), int(t), weight=float(1.0 / (w + 1e-9)))
    # Largest Connected Component
    if not nx.is_connected(G):
        largest = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest).copy()
    return G


def compute_strength(G: nx.Graph) -> np.ndarray:
    return np.array([w for _, w in G.degree(weight="weight")], dtype=float)


def compute_centralities(G: nx.Graph) -> Dict[str, np.ndarray]:
    centralities = {}
    centralities["eigenvector"] = np.array(
        list(nx.eigenvector_centrality(G, max_iter=1000, weight="weight").values())
    )
    centralities["degree"] = np.array([v for _, v in G.degree(weight="weight")], dtype=float)
    # closeness with distance = 1/weight (so we need lengths)
    inv_len = {(u, v): 1.0 / data.get("weight", 1.0) for u, v, data in G.edges(data=True)}
    nx.set_edge_attributes(G, inv_len, name="length")
    centralities["closeness"] = np.array(
        list(nx.closeness_centrality(G, distance="length").values()), dtype=float
    )
    return centralities


def bootstrap_correlations(
    strength: np.ndarray, centrality: np.ndarray, n_boot: int = 50
) -> (float, float):
    n = len(strength)
    if n < 50:
        return np.nan, np.nan
    corrs = []
    for _ in range(n_boot):
        idx = np.random.randint(0, n, size=n)
        s = strength[idx]
        c = centrality[idx]
        if np.std(s) == 0 or np.std(c) == 0:
            continue
        corrs.append(np.corrcoef(s, c)[0, 1])
    if not corrs:
        return np.nan, np.nan
    return float(np.mean(corrs)), float(np.std(corrs))


def main():
    data_path = os.path.join(os.path.dirname(__file__), "..", "data", "sdss_100k_galaxy_form_burst.csv")
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data not found at {data_path}")

    df = load_sdss_csv(data_path)
    if len(df) < 200:
        raise ValueError("Dataset too small for sensitivity analysis.")

    coords = to_cartesian(df)
    k_values = [5, 10, 15, 20, 30]
    results: Dict[str, List[Dict[str, float]]] = {}

    for cent_name in ["eigenvector", "degree", "closeness"]:
        results[cent_name] = []

    for k in tqdm(k_values, desc="k-scan"):
        G = build_graph(coords, k)
        strength = compute_strength(G)
        centrals = compute_centralities(G)
        for cname, cvals in centrals.items():
            mean_corr, std_corr = bootstrap_correlations(strength, cvals)
            results[cname].append(
                {"k": k, "corr_mean": mean_corr, "corr_std": std_corr, "n_nodes": len(G)}
            )

    out_path = os.path.join(os.path.dirname(__file__), "..", "outputs", "sensitivity_stats.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Saved sensitivity stats to {out_path}")

    # Optional quick text summary
    for cname, arr in results.items():
        print(f"\nCentrality: {cname}")
        for row in arr:
            print(
                f"k={row['k']}: corr={row['corr_mean']:.3f}±{row['corr_std']:.3f} (N={row['n_nodes']})"
            )


if __name__ == "__main__":
    main()
