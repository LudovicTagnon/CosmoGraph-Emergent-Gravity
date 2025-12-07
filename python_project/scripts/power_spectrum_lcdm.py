import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.neighbors import kneighbors_graph

DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "sdss_100k_galaxy_form_burst.csv"
OUT_IMG = Path(__file__).resolve().parent.parent / "outputs" / "power_with_lcdm.png"


def load_sdss_slice(max_rows=5000):
    df = pd.read_csv(DATA_PATH, comment="#", low_memory=False)
    zcol = "redshift" if "redshift" in df.columns else "z"
    df = df[["ra", "dec", zcol]].rename(columns={zcol: "z"})
    m = (df["z"].between(0.04, 0.12)) & (df["ra"].between(130, 240)) & (df["dec"].between(-5, 60))
    df = df[m]
    if len(df) > max_rows:
        df = df.sort_values("ra").head(max_rows)
    return df.reset_index(drop=True)


def ra_dec_z_to_cartesian(df, H0=70.0, c=300000.0):
    dist = (df["z"].to_numpy() * c) / H0
    ra = np.deg2rad(df["ra"].to_numpy())
    dec = np.deg2rad(df["dec"].to_numpy())
    x = dist * np.cos(dec) * np.cos(ra)
    y = dist * np.cos(dec) * np.sin(ra)
    zc = dist * np.sin(dec)
    coords = np.vstack([x, y, zc]).T
    coords -= coords.mean(axis=0, keepdims=True)
    return coords


def build_strength(coords, k=10):
    A = kneighbors_graph(coords, n_neighbors=k, mode="distance", include_self=False)
    A.data = 1.0 / (A.data + 1e-9)
    return np.array(A.sum(axis=1)).ravel()


def compute_pk(points, weights, grid=64):
    mins = points.min(axis=0); maxs = points.max(axis=0)
    span = maxs - mins; span[span == 0] = 1.0
    box = span.max() * 1.05
    norm = (points - mins) / box
    idx = np.clip(np.floor(norm * grid).astype(int), 0, grid - 1)
    rho = np.zeros((grid, grid, grid), float)
    for (i, j, k), w in zip(idx, weights):
        rho[i, j, k] += w
    mean = rho.mean()
    delta = (rho - mean) / mean
    dk = np.fft.fftn(delta)
    P = np.abs(dk) ** 2
    kfreq = np.fft.fftfreq(grid) * 2 * np.pi
    kx, ky, kz = np.meshgrid(kfreq, kfreq, kfreq, indexing="ij")
    kmag = np.sqrt(kx**2 + ky**2 + kz**2).ravel()
    P = P.ravel()
    mask = kmag > 0
    kmag = kmag[mask]; P = P[mask]
    bins = np.logspace(np.log10(kmag.min()), np.log10(kmag.max()), 30)
    kc, pk = [], []
    for a, b in zip(bins[:-1], bins[1:]):
        m = (kmag >= a) & (kmag < b)
        if np.any(m):
            kc.append(np.mean(kmag[m])); pk.append(np.mean(P[m]))
    return np.array(kc), np.array(pk)


def fit_slope(k, p, kmin=0.05, kmax=0.3):
    m = (k >= kmin) & (k <= kmax) & (p > 0)
    if m.sum() < 5:
        return np.nan
    coef = np.polyfit(np.log10(k[m]), np.log10(p[m]), 1)
    return coef[0]


def lcdm_reference(k, ns=0.96):
    # Simple power-law reference P ~ k^{ns-4}; normalized to match mass spectrum at mid-k
    slope = ns - 4  # ~ -3.04
    return k ** slope


def main():
    df = load_sdss_slice()
    coords = ra_dec_z_to_cartesian(df)
    N = len(coords)
    rng = np.random.default_rng(42)
    box_size = (coords.max(axis=0) - coords.min(axis=0)).max() * 1.05
    rand_coords = rng.uniform(-box_size / 2, box_size / 2, size=(N, 3))

    strength = build_strength(coords, k=10)
    strength_rand = build_strength(rand_coords, k=10)
    w_m = np.ones(N)
    w_t = strength / strength.mean()
    w_m_r = np.ones(N)
    w_t_r = strength_rand / strength_rand.mean()

    k_m, Pm = compute_pk(coords, w_m, grid=64)
    k_t, Pt = compute_pk(coords, w_t, grid=64)
    k_mr, Pmr = compute_pk(rand_coords, w_m_r, grid=64)
    k_tr, Ptr = compute_pk(rand_coords, w_t_r, grid=64)

    # match LCDM amplitude at median k of mass spectrum
    k_ref = k_m
    lcdm = lcdm_reference(k_ref)
    # scale to match Pm at mid k
    mid = len(k_ref) // 2
    if lcdm[mid] != 0:
        scale = Pm[mid] / lcdm[mid]
        lcdm *= scale

    slope_m = fit_slope(k_m, Pm)
    slope_t = fit_slope(k_t, Pt)
    Tm = Pm / Pmr
    Tt = Pt / Ptr
    slope_Tm = fit_slope(k_m, Tm)
    slope_Tt = fit_slope(k_t, Tt)

    plt.figure(figsize=(10, 5))
    plt.loglog(k_m, Pm, 'o-', label=f'Mass (slope {slope_m:.2f})', color='orange')
    plt.loglog(k_t, Pt, 'x--', label=f'Topo (slope {slope_t:.2f})', color='blue')
    plt.loglog(k_ref, lcdm, ':', label='LCDM-like ref (ns~0.96)', color='gray')
    plt.xlabel('k'); plt.ylabel('P(k)'); plt.title('P(k) with LCDM-like reference')
    plt.legend(); plt.grid(True, which='both', alpha=0.3)
    plt.tight_layout(); plt.savefig(OUT_IMG, dpi=200)
    print(f"Saved {OUT_IMG}")

    # save arrays for reuse
    np.savez(OUT_IMG.with_suffix('.npz'), k_m=k_m, Pm=Pm, k_t=k_t, Pt=Pt, k_mr=k_mr, Pmr=Pmr, k_tr=k_tr, Ptr=Ptr, lcdm=k_ref)

if __name__ == "__main__":
    main()
