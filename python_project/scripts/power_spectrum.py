import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "sdss_100k_galaxy_form_burst.csv"
OUT_IMG = Path(__file__).resolve().parent.parent / "outputs" / "power_spectrum_compare.png"


def load_sdss_slice(max_rows=20000):
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Missing SDSS file at {DATA_PATH}")
    df = pd.read_csv(DATA_PATH, comment="#")
    # column names may vary; try common ones
    zcol = "redshift" if "redshift" in df.columns else "z"
    if not {"ra", "dec", zcol}.issubset(df.columns):
        raise ValueError("SDSS file missing required columns ra, dec, z/redshift")
    df = df[["ra", "dec", zcol]].rename(columns={zcol: "z"})
    # slice Grande Muraille
    m = (df["z"].between(0.04, 0.12)) & (df["ra"].between(130, 240)) & (df["dec"].between(-5, 60))
    df = df[m]
    if len(df) > max_rows:
        df = df.sort_values("ra").head(max_rows)
    return df


def ra_dec_z_to_cartesian(df, H0=70.0, c=300000.0):
    dist = (df["z"].to_numpy() * c) / H0  # Mpc/h approximate
    ra = np.deg2rad(df["ra"].to_numpy())
    dec = np.deg2rad(df["dec"].to_numpy())
    x = dist * np.cos(dec) * np.cos(ra)
    y = dist * np.cos(dec) * np.sin(ra)
    zc = dist * np.sin(dec)
    coords = np.vstack([x, y, zc]).T
    # recenter
    coords -= coords.mean(axis=0, keepdims=True)
    return coords


def make_toy_points(n_total=3000, box=200.0, seed=42):
    rng = np.random.default_rng(seed)
    centers = np.array([[ -40,  10,  20],
                        [  50, -20, -30],
                        [ -10, -50,  40]])
    pts = []
    n_cluster = n_total // len(centers)
    for c in centers:
        pts.append(rng.normal(loc=c, scale=15.0, size=(n_cluster, 3)))
    pts = np.vstack(pts)
    # ensure within box
    pts = np.clip(pts, -box/2, box/2)
    return pts


def voxelize(points, grid=64):
    # normalize to [0,1]^3 then to indices
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    span = maxs - mins
    span[span == 0] = 1.0
    norm = (points - mins) / span
    idx = np.floor(norm * grid).astype(int)
    idx = np.clip(idx, 0, grid - 1)
    rho = np.zeros((grid, grid, grid), dtype=np.float64)
    for i, j, k in idx:
        rho[i, j, k] += 1.0
    return rho


def power_spectrum(rho):
    mean = rho.mean()
    if mean == 0:
        return None, None
    delta = (rho - mean) / mean
    delta_k = np.fft.fftn(delta)
    power = np.abs(delta_k) ** 2
    kfreq = np.fft.fftfreq(rho.shape[0])
    kx, ky, kz = np.meshgrid(kfreq, kfreq, kfreq, indexing="ij")
    kmag = np.sqrt(kx**2 + ky**2 + kz**2).ravel()
    power_flat = power.ravel()
    # binning
    mask = kmag > 0
    kmag = kmag[mask]
    power_flat = power_flat[mask]
    bins = np.linspace(0, kmag.max(), 50)
    digitized = np.digitize(kmag, bins)
    pk = np.array([power_flat[digitized == i].mean() for i in range(1, len(bins))])
    kbin = 0.5 * (bins[:-1] + bins[1:])
    # remove NaN
    m = np.isfinite(pk) & (pk > 0)
    return kbin[m], pk[m]


def main():
    # SDSS
    df_sdss = load_sdss_slice()
    coords_sdss = ra_dec_z_to_cartesian(df_sdss)
    rho_sdss = voxelize(coords_sdss, grid=64)
    k_sdss, p_sdss = power_spectrum(rho_sdss)

    # Toy
    pts_toy = make_toy_points()
    rho_toy = voxelize(pts_toy, grid=64)
    k_toy, p_toy = power_spectrum(rho_toy)

    plt.figure(figsize=(10, 5))
    ax1 = plt.subplot(1, 2, 1)
    if k_toy is not None:
        ax1.loglog(k_toy, p_toy, label="Toy Universe", color="cyan")
    ax1.set_xlabel("k")
    ax1.set_ylabel("P(k)")
    ax1.set_title("Toy Power Spectrum")
    ax1.grid(True, which="both", ls="--", alpha=0.3)
    ax1.legend()

    ax2 = plt.subplot(1, 2, 2)
    if k_sdss is not None:
        ax2.loglog(k_sdss, p_sdss, label="SDSS Slice", color="orange")
    ax2.set_xlabel("k")
    ax2.set_ylabel("P(k)")
    ax2.set_title("SDSS Power Spectrum (slice)")
    ax2.grid(True, which="both", ls="--", alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    OUT_IMG.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_IMG, dpi=200)
    print(f"Saved {OUT_IMG}")
    if k_sdss is not None:
        print(f"SDSS points: {len(df_sdss)}, k-range [{k_sdss.min():.4f}, {k_sdss.max():.4f}]")
    if k_toy is not None:
        print(f"Toy points: {len(pts_toy)}, k-range [{k_toy.min():.4f}, {k_toy.max():.4f}]")


if __name__ == "__main__":
    main()
