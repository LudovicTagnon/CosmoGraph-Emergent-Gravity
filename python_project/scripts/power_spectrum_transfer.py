import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.neighbors import kneighbors_graph

DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "sdss_100k_galaxy_form_burst.csv"
OUT_IMG = Path(__file__).resolve().parent.parent / "outputs" / "transfer_function_spectrum.png"


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
    coords -= coords.mean(axis=0, keepdims=True)
    return coords


def build_strength_weights(coords, k=10):
    A = kneighbors_graph(coords, n_neighbors=k, mode="distance", include_self=False)
    A.data = 1.0 / (A.data + 1e-9)
    strength = np.array(A.sum(axis=1)).ravel()
    return strength


def compute_weighted_power_spectrum(points, weights, grid=64):
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


def main():
    # Données réelles
    df = load_sdss_slice()
    coords = ra_dec_z_to_cartesian(df)
    N = len(coords)
    # Positions aléatoires (même cube)
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    box_size = (maxs - mins).max() * 1.05
    rng = np.random.default_rng(42)
    rand_coords = rng.uniform(-box_size / 2, box_size / 2, size=(N, 3))

    # Poids masse / topo (strength k=10)
    strength_real = build_strength_weights(coords, k=10)
    strength_rand = build_strength_weights(rand_coords, k=10)
    w_mass_real = np.ones(N)
    w_topo_real = strength_real / strength_real.mean()
    w_mass_rand = np.ones(N)
    w_topo_rand = strength_rand / strength_rand.mean()

    # Spectres
    k_m_r, p_m_r = compute_weighted_power_spectrum(coords, w_mass_real, grid=64)
    k_t_r, p_t_r = compute_weighted_power_spectrum(coords, w_topo_real, grid=64)
    k_m_n, p_m_n = compute_weighted_power_spectrum(rand_coords, w_mass_rand, grid=64)
    k_t_n, p_t_n = compute_weighted_power_spectrum(rand_coords, w_topo_rand, grid=64)

    # Fonction de transfert (rapport réel / aléatoire, même k-bins)
    # Les binnings sont identiques par construction (mêmes dimensions), on suppose mêmes longueurs.
    T_m = p_m_r / p_m_n
    T_t = p_t_r / p_t_n

    slope_Tm = fit_slope(k_m_r, T_m)
    slope_Tt = fit_slope(k_t_r, T_t)

    plt.figure(figsize=(8, 6))
    plt.loglog(k_m_r, T_m, "o-", label=f"Matière nettoyée (slope {slope_Tm:.2f})", color="orange")
    plt.loglog(k_t_r, T_t, "x--", label=f"Topo nettoyée (slope {slope_Tt:.2f})", color="blue")
    plt.xlabel("k")
    plt.ylabel("T(k) = P_real / P_random")
    plt.title("Fonction de transfert (signal / bruit de grille)")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3, ls="--")
    plt.tight_layout()
    OUT_IMG.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_IMG, dpi=200)

    print(f"Saved {OUT_IMG}")
    print(f"N points: {N}")
    print(f"Slope T_mass : {slope_Tm:.3f}")
    print(f"Slope T_topo : {slope_Tt:.3f}")
    if np.isfinite(slope_Tm) and np.isfinite(slope_Tt):
        print(f"Delta slope (topo - mass) : {slope_Tt - slope_Tm:.3f}")


if __name__ == "__main__":
    main()
