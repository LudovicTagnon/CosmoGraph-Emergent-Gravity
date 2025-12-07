import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr
from sklearn.neighbors import kneighbors_graph


DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "sdss_100k_galaxy_form_burst.csv"
OUT_IMG = Path(__file__).resolve().parent.parent / "outputs" / "centrality_vs_matter_spectrum.png"


def load_sdss_slice(max_rows=20000):
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Missing SDSS file at {DATA_PATH}")
    df = pd.read_csv(DATA_PATH, comment="#")
    zcol = "redshift" if "redshift" in df.columns else "z"
    if not {"ra", "dec", zcol}.issubset(df.columns):
        raise ValueError("SDSS file missing required columns ra, dec, z/redshift")
    df = df[["ra", "dec", zcol]].rename(columns={zcol: "z"})
    # Slice Grande Muraille
    m = (df["z"].between(0.04, 0.12)) & (df["ra"].between(130, 240)) & (df["dec"].between(-5, 60))
    df = df[m]
    # Cap pour maîtrise du calcul
    if len(df) > max_rows:
        df = df.sort_values("ra").head(max_rows)
    return df.reset_index(drop=True)


def ra_dec_z_to_cartesian(df, H0=70.0, c=300000.0):
    dist = (df["z"].to_numpy() * c) / H0  # Mpc/h approx
    ra = np.deg2rad(df["ra"].to_numpy())
    dec = np.deg2rad(df["dec"].to_numpy())
    x = dist * np.cos(dec) * np.cos(ra)
    y = dist * np.cos(dec) * np.sin(ra)
    zc = dist * np.sin(dec)
    coords = np.vstack([x, y, zc]).T
    # recentrer
    coords -= coords.mean(axis=0, keepdims=True)
    return coords


def build_strength_weights(coords, k=10):
    # k-NN distances, mode='distance'
    A = kneighbors_graph(coords, n_neighbors=k, mode="distance", include_self=False)
    # convertir en poids 1/d
    A.data = 1.0 / (A.data + 1e-9)
    # strength = somme des poids
    strength = np.array(A.sum(axis=1)).ravel()
    return strength


def compute_weighted_power_spectrum(points, weights, grid=64):
    # normalisation dans un cube [0, box]^3
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    span = maxs - mins
    span[span == 0] = 1.0
    box_size = span.max() * 1.05
    norm = (points - mins) / box_size  # [0,1]
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
    kfreq = np.fft.fftfreq(grid) * 2 * np.pi  # en unités de box^{-1}
    kx, ky, kz = np.meshgrid(kfreq, kfreq, kfreq, indexing="ij")
    kmag = np.sqrt(kx**2 + ky**2 + kz**2).ravel()
    power = power.ravel()
    # bin log
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


def main():
    df = load_sdss_slice()
    coords = ra_dec_z_to_cartesian(df)
    strength = build_strength_weights(coords, k=10)
    weights_matter = np.ones(len(coords))
    weights_topo = strength / strength.mean()

    k_m, p_m = compute_weighted_power_spectrum(coords, weights_matter, grid=64)
    k_t, p_t = compute_weighted_power_spectrum(coords, weights_topo, grid=64)

    slope_m = fit_slope(k_m, p_m)
    slope_t = fit_slope(k_t, p_t)

    # biais moyen P_topo / P_matter
    common = (k_m is not None) and (k_t is not None)
    bias = np.nan
    if common:
        # align by nearest k (same bins built, so lengths match)
        msk = (p_m > 0) & (p_t > 0) & np.isfinite(p_m) & np.isfinite(p_t)
        if msk.sum() > 0:
            bias = np.mean(np.sqrt(p_t[msk] / p_m[msk]))

    plt.figure(figsize=(8, 6))
    plt.loglog(k_m, p_m, "o-", label=f"Masse (slope {slope_m:.2f})", color="orange")
    plt.loglog(k_t, p_t, "x--", label=f"Centralité/Force (slope {slope_t:.2f})", color="blue")
    plt.xlabel("k")
    plt.ylabel("P(k)")
    plt.title(f"Substitution de masse : biais b ≈ {bias:.2f}")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3, ls="--")
    plt.tight_layout()
    OUT_IMG.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_IMG, dpi=200)
    print(f"Saved {OUT_IMG}")
    print(f"SDSS slice points: {len(coords)}")
    print(f"Slope mass   : {slope_m:.3f}")
    print(f"Slope topo   : {slope_t:.3f}")
    print(f"Bias (mean sqrt P_topo/P_mass): {bias:.3f}")


if __name__ == "__main__":
    main()
