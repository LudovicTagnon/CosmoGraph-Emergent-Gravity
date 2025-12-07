import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import kneighbors_graph


DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "sdss_100k_galaxy_form_burst.csv"
OUT_JSON = Path(__file__).resolve().parent.parent / "outputs" / "jackknife_power.json"


def load_sdss_slice(max_rows=8000):
    """Load the SDSS slice used in the analysis and return a DataFrame.

    We keep a contiguous volume in (ra, dec, z) similar to the main pipeline,
    then truncate to `max_rows` (sorted by ra) to keep the runtime reasonable.
    """

    df = pd.read_csv(DATA_PATH, comment="#")
    zcol = "redshift" if "redshift" in df.columns else "z"
    if not {"ra", "dec", zcol}.issubset(df.columns):
        raise ValueError("SDSS file missing required columns ra, dec, z/redshift")

    df = df[["ra", "dec", zcol]].rename(columns={zcol: "z"})
    m = (df["z"].between(0.04, 0.12)) & (df["ra"].between(130, 240)) & (df["dec"].between(-5, 60))
    df = df[m].copy()
    if len(df) > max_rows:
        df = df.sort_values("ra").head(max_rows)
    return df.reset_index(drop=True)


def to_cartesian(df, H0=70.0, c=300000.0):
    dist = (df["z"].to_numpy() * c) / H0
    ra = np.deg2rad(df["ra"].to_numpy())
    dec = np.deg2rad(df["dec"].to_numpy())
    x = dist * np.cos(dec) * np.cos(ra)
    y = dist * np.cos(dec) * np.sin(ra)
    zc = dist * np.sin(dec)
    coords = np.vstack([x, y, zc]).T
    coords -= coords.mean(axis=0, keepdims=True)
    return coords


def build_strength(points, k=10):
    # k-NN graph; weights are inverse distances; strength = weighted degree
    A = kneighbors_graph(points, n_neighbors=k, mode="distance", include_self=False)
    A.data = 1.0 / (A.data + 1e-9)
    strength = np.array(A.sum(axis=1)).ravel()
    return strength


def power_spectrum(points, weights, grid=64):
    # Voxelise and compute isotropised P(k)
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    span = maxs - mins
    span[span == 0] = 1.0
    box_size = span.max() * 1.05
    norm = (points - mins) / box_size
    idx = np.floor(norm * grid).astype(int)
    idx = np.clip(idx, 0, grid - 1)

    rho = np.zeros((grid, grid, grid), dtype=np.float64)
    for (i, j, k), w in zip(idx, weights):
        rho[i, j, k] += w

    mean = rho.mean()
    if mean == 0:
        return None, None
    delta = (rho - mean) / mean
    delta_k = np.fft.fftn(delta)
    power = np.abs(delta_k) ** 2

    kfreq = np.fft.fftfreq(grid) * 2 * np.pi
    kx, ky, kz = np.meshgrid(kfreq, kfreq, kfreq, indexing="ij")
    kmag = np.sqrt(kx**2 + ky**2 + kz**2).ravel()
    power = power.ravel()
    mask = kmag > 0
    kmag = kmag[mask]
    power = power[mask]

    bins = np.logspace(np.log10(kmag.min()), np.log10(kmag.max()), 30)
    k_centers, pk_vals = [], []
    for a, b in zip(bins[:-1], bins[1:]):
        m = (kmag >= a) & (kmag < b)
        if np.any(m):
            k_centers.append(np.mean(kmag[m]))
            pk_vals.append(np.mean(power[m]))
    return np.array(k_centers), np.array(pk_vals)


def fit_slope(k, p, kmin=0.05, kmax=0.3):
    m = (k >= kmin) & (k <= kmax) & (p > 0)
    if m.sum() < 5:
        return np.nan
    coef = np.polyfit(np.log10(k[m]), np.log10(p[m]), 1)
    return coef[0]


def jackknife_slope(coords, k_neighbors=10, grid=64, kmin=0.05, kmax=0.3):
    """Jackknife on (ra,dec,z) sub-volumes defined by quantiles.

    Returns dict with mean slopes and 1-sigma jackknife errors for
    P_mass, P_topo, T_mass, T_topo.
    """

    N = len(coords)
    # Split RA/DEC/Z into two quantile bins each (2x2x2 = 8 chunks)
    # We work from the original ra/dec/z used to build coords
    # To avoid recomputing, keep a copy alongside coords
    raise NotImplementedError


def jackknife_from_dataframe(df, k_neighbors=10, grid=64, kmin=0.05, kmax=0.3):
    """Perform jackknife leaving out each RA/DEC/Z quantile chunk."""

    # Build quantile edges
    ra_edges = np.quantile(df["ra"], [0, 0.5, 1])
    dec_edges = np.quantile(df["dec"], [0, 0.5, 1])
    z_edges = np.quantile(df["z"], [0, 0.5, 1])

    chunks = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                mask = (
                    (df["ra"].between(ra_edges[i], ra_edges[i + 1], inclusive="right" if i == 0 else "both"))
                    & (df["dec"].between(dec_edges[j], dec_edges[j + 1], inclusive="right" if j == 0 else "both"))
                    & (df["z"].between(z_edges[k], z_edges[k + 1], inclusive="right" if k == 0 else "both"))
                )
                chunks.append(mask)

    slopes = []
    for idx, mask in enumerate(chunks):
        keep = ~mask
        df_sub = df[keep]
        if len(df_sub) < 300:
            slopes.append({"P_mass": np.nan, "P_topo": np.nan, "T_mass": np.nan, "T_topo": np.nan})
            continue

        pts = to_cartesian(df_sub)
        strength = build_strength(pts, k=k_neighbors)
        w_mass = np.ones(len(df_sub))
        w_topo = strength / strength.mean()

        k_m, p_m = power_spectrum(pts, w_mass, grid=grid)
        k_t, p_t = power_spectrum(pts, w_topo, grid=grid)

        # Random catalogue matching box
        mins = pts.min(axis=0)
        maxs = pts.max(axis=0)
        box_size = (maxs - mins).max() * 1.05
        rng = np.random.default_rng(1234 + idx)
        rand_pts = rng.uniform(-box_size / 2, box_size / 2, size=pts.shape)
        strength_r = build_strength(rand_pts, k=k_neighbors)
        w_m_r = np.ones(len(rand_pts))
        w_t_r = strength_r / strength_r.mean()
        k_m_r, p_m_r = power_spectrum(rand_pts, w_m_r, grid=grid)
        k_t_r, p_t_r = power_spectrum(rand_pts, w_t_r, grid=grid)

        if any(x is None for x in [k_m, p_m, k_t, p_t, k_m_r, p_m_r, k_t_r, p_t_r]):
            slopes.append({"P_mass": np.nan, "P_topo": np.nan, "T_mass": np.nan, "T_topo": np.nan})
            continue

        T_m = p_m / p_m_r
        T_t = p_t / p_t_r

        slopes.append(
            {
                "P_mass": fit_slope(k_m, p_m, kmin=kmin, kmax=kmax),
                "P_topo": fit_slope(k_t, p_t, kmin=kmin, kmax=kmax),
                "T_mass": fit_slope(k_m, T_m, kmin=kmin, kmax=kmax),
                "T_topo": fit_slope(k_t, T_t, kmin=kmin, kmax=kmax),
            }
        )

    # Jackknife combine
    def jackknife_mean_error(values):
        vals = np.array(values, dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            return np.nan, np.nan
        mean = vals.mean()
        err = np.sqrt(((len(vals) - 1) / len(vals)) * np.sum((vals - mean) ** 2))
        return mean, err

    result = {key: jackknife_mean_error([s[key] for s in slopes]) for key in ["P_mass", "P_topo", "T_mass", "T_topo"]}
    result["n_chunks"] = len(slopes)
    result["chunks_used"] = int(sum(np.isfinite([s["P_mass"] for s in slopes])))
    return result


def main():
    df = load_sdss_slice()
    res = jackknife_from_dataframe(df, k_neighbors=10, grid=64, kmin=0.05, kmax=0.3)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
