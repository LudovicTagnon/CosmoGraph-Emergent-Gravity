import numpy as np
from numpy.fft import fftn, ifftn

# Small reproducible sandbox to compare proxies (geometry vs entropie locale)
# on synthetic 1/f anisotropic fields with periodic Poisson solver.


def generate_1overf_field(size=32, alpha=1.8, anisotropy=(1.0, 0.8, 0.6), seed=0):
    rng = np.random.default_rng(seed)
    kx = np.fft.fftfreq(size)[:, None, None]
    ky = np.fft.fftfreq(size)[None, :, None]
    kz = np.fft.fftfreq(size)[None, None, :]
    # Apply anisotropy by scaling axes
    ax, ay, az = anisotropy
    k2 = (kx / ax) ** 2 + (ky / ay) ** 2 + (kz / az) ** 2
    k2[0, 0, 0] = 1.0  # avoid div by zero
    amp = 1.0 / np.power(np.sqrt(k2), alpha / 2.0)
    phases = rng.normal(size=(size, size, size)) + 1j * rng.normal(size=(size, size, size))
    rho = np.real(ifftn(amp * phases))
    rho -= rho.min()
    rho /= rho.max() + 1e-9
    return rho


def poisson_potential(rho):
    n = rho.shape[0]
    rho_k = fftn(rho)
    kx = np.fft.fftfreq(n)[:, None, None]
    ky = np.fft.fftfreq(n)[None, :, None]
    kz = np.fft.fftfreq(n)[None, None, :]
    k2 = kx * kx + ky * ky + kz * kz
    k2[0, 0, 0] = 1.0  # avoid divide by zero for DC
    phi_k = -4 * np.pi * (rho_k - rho_k[0, 0, 0]) / k2
    phi_k[0, 0, 0] = 0.0
    phi = np.real(ifftn(phi_k))
    return phi


def grad_norm(phi):
    gx, gy, gz = np.gradient(phi)
    return np.sqrt(gx * gx + gy * gy + gz * gz)


def laplacian(f):
    return (
        -6 * f
        + np.roll(f, 1, 0)
        + np.roll(f, -1, 0)
        + np.roll(f, 1, 1)
        + np.roll(f, -1, 1)
        + np.roll(f, 1, 2)
        + np.roll(f, -1, 2)
    )


def local_entropy(field, bins=16):
    from numpy.lib.stride_tricks import sliding_window_view

    window = sliding_window_view(field, (3, 3, 3))
    patches = window.reshape(-1, 27)
    # Quantize
    scaled = np.clip((field - field.min()) / (field.max() - field.min() + 1e-9), 0, 1)
    scaled_win = sliding_window_view(scaled, (3, 3, 3)).reshape(-1, 27)
    digitized = np.clip((scaled_win * bins).astype(int), 0, bins - 1)
    ent = np.empty(digitized.shape[0], dtype=np.float32)
    for i, row in enumerate(digitized):
        counts = np.bincount(row, minlength=bins)
        p = counts / counts.sum()
        mask = p > 0
        ent[i] = -np.sum(p[mask] * np.log(p[mask]))
    # center crop to match valid region
    n = field.shape[0]
    valid = (slice(1, n - 1), slice(1, n - 1), slice(1, n - 1))
    return ent.reshape(n - 2, n - 2, n - 2)


def make_dataset(rho):
    n = rho.shape[0]
    phi = poisson_potential(rho)
    y = grad_norm(phi)
    # center region to avoid boundary artifacts for entropy
    core = (slice(1, n - 1), slice(1, n - 1), slice(1, n - 1))
    m = rho[core]
    y = y[core]
    grad2 = sum(np.gradient(rho)[i][core] ** 2 for i in range(3))
    lap = laplacian(rho)[core]
    ent = local_entropy(rho)
    X_baseline = np.stack([m.ravel()], axis=1)
    X_geom = np.stack([m.ravel(), grad2.ravel(), lap.ravel()], axis=1)
    X_ent = np.stack([m.ravel(), ent.ravel(), lap.ravel()], axis=1)
    return X_baseline, X_geom, X_ent, y.ravel()


def train_test_split(X, y, rng, test_ratio=0.4):
    idx = np.arange(len(y))
    rng.shuffle(idx)
    t = int(len(y) * (1 - test_ratio))
    train, test = idx[:t], idx[t:]
    return X[train], X[test], y[train], y[test]


def fit_ols(X, y):
    Xb = np.c_[X, np.ones(len(X))]
    coef, *_ = np.linalg.lstsq(Xb, y, rcond=None)
    return coef


def predict(X, coef):
    Xb = np.c_[X, np.ones(len(X))]
    return Xb @ coef


def metrics(y_true, y_pred):
    yt = y_true
    yp = y_pred
    yt_center = yt - yt.mean()
    yp_center = yp - yp.mean()
    denom = np.sqrt((yt_center ** 2).sum() * (yp_center ** 2).sum()) + 1e-12
    pearson = float((yt_center @ yp_center) / denom)
    rmse = np.sqrt(np.mean((yt - yp) ** 2))
    relrmse = float(rmse / (np.mean(np.abs(yt)) + 1e-12))
    return pearson, relrmse


def run_seed(seed=0, size=32, alpha=1.8, anisotropy=(1.0, 0.8, 0.6)):
    rng = np.random.default_rng(seed)
    rho = generate_1overf_field(size=size, alpha=alpha, anisotropy=anisotropy, seed=seed)
    Xb, Xg, Xe, y = make_dataset(rho)
    out = {}
    for name, X in [
        ("baseline_m", Xb),
        ("geom_m_grad2_lap", Xg),
        ("ent_m_entropy_lap", Xe),
    ]:
        Xtr, Xte, ytr, yte = train_test_split(X, y, rng)
        coef = fit_ols(Xtr, ytr)
        pred = predict(Xte, coef)
        p, r = metrics(yte, pred)
        out[name] = {
            "pearson": p,
            "relrmse": r,
            "coef": coef.tolist(),
        }
    return out


def aggregate(seeds=(0, 1, 2), **kwargs):
    rows = []
    for s in seeds:
        res = run_seed(seed=s, **kwargs)
        for k, v in res.items():
            rows.append((k, s, v["pearson"], v["relrmse"], v["coef"]))
    return rows


def summarize(rows):
    import math
    summary = {}
    for name in set(r[0] for r in rows):
        vals_p = [r[2] for r in rows if r[0] == name]
        vals_r = [r[3] for r in rows if r[0] == name]
        summary[name] = {
            "pearson_mean": float(np.mean(vals_p)),
            "pearson_std": float(np.std(vals_p)),
            "relrmse_mean": float(np.mean(vals_r)),
            "relrmse_std": float(np.std(vals_r)),
            "n": len(vals_p),
        }
    return summary


def main():
    seeds = list(range(3))
    rows = aggregate(seeds=seeds, size=32, alpha=1.8, anisotropy=(1.0, 0.8, 0.6))
    summary = summarize(rows)
    print("# Entropie vs géométrie — 1/f anisotrope size=32, alpha=1.8, seeds=", len(seeds))
    for name, stats in summary.items():
        print(f"{name}: pearson={stats['pearson_mean']:.6f} ±{stats['pearson_std']:.6f} ; "
              f"relRMSE={stats['relrmse_mean']:.6f} ±{stats['relrmse_std']:.6f} (n={stats['n']})")

    # Write detailed rows for traceability
    import json
    out = {
        "setup": {
            "size": 32,
            "alpha": 1.8,
            "anisotropy": (1.0, 0.8, 0.6),
            "seeds": seeds,
        },
        "runs": rows,
        "summary": summary,
    }
    with open("IA5/entropy_results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
