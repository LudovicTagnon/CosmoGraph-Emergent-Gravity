import numpy as np


def idx(i: int, j: int, k: int, ny: int, nz: int) -> int:
    return (i * ny + j) * nz + k


def build_K3d(
    nx: int,
    ny: int,
    nz: int,
    k_coupling: float = 1.0,
    mass: float = 0.0,
    defect: tuple[int, int, int] | None = None,
    defect_mass: float | None = None,
    periodic: bool = True,
) -> np.ndarray:
    n = nx * ny * nz
    K = np.zeros((n, n), dtype=np.float64)

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                p = idx(i, j, k, ny, nz)
                diag = mass * mass
                neighbors = []
                if periodic:
                    neighbors = [
                        ((i + 1) % nx, j, k),
                        ((i - 1) % nx, j, k),
                        (i, (j + 1) % ny, k),
                        (i, (j - 1) % ny, k),
                        (i, j, (k + 1) % nz),
                        (i, j, (k - 1) % nz),
                    ]
                else:
                    if i > 0:
                        neighbors.append((i - 1, j, k))
                    if i < nx - 1:
                        neighbors.append((i + 1, j, k))
                    if j > 0:
                        neighbors.append((i, j - 1, k))
                    if j < ny - 1:
                        neighbors.append((i, j + 1, k))
                    if k > 0:
                        neighbors.append((i, j, k - 1))
                    if k < nz - 1:
                        neighbors.append((i, j, k + 1))
                for (ii, jj, kk) in neighbors:
                    q = idx(ii, jj, kk, ny, nz)
                    K[p, q] = -k_coupling
                    diag += k_coupling
                K[p, p] = diag

    if defect is not None and defect_mass is not None:
        di, dj, dk = defect
        p = idx(di, dj, dk, ny, nz)
        diag = defect_mass * defect_mass
        diag += 6 * k_coupling  # six neighbors
        K[p, p] = diag
    return K


def ground_state_covariance(K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    w, v = np.linalg.eigh(K)
    eps = 1e-8
    mask = w > eps
    w_nz = w[mask]
    v_nz = v[:, mask]
    w_sqrt = np.sqrt(w_nz)
    w_inv_sqrt = 1.0 / w_sqrt
    Cx = v_nz @ np.diag(0.5 * w_inv_sqrt) @ v_nz.T
    Cp = v_nz @ np.diag(0.5 * w_sqrt) @ v_nz.T
    return Cx, Cp


def mutual_information(i: int, j: int, Cx: np.ndarray, Cp: np.ndarray, scale: float) -> float:
    Vx_i = Cx[np.ix_([i], [i])]
    Vp_i = Cp[np.ix_([i], [i])]
    Vx_j = Cx[np.ix_([j], [j])]
    Vp_j = Cp[np.ix_([j], [j])]
    Vx_ij = Cx[np.ix_([i, j], [i, j])]
    Vp_ij = Cp[np.ix_([i, j], [i, j])]
    def symp_from_xp(Vx, Vp):
        vals = np.linalg.eigvals(Vx @ Vp).real
        vals = np.clip(vals, 0, None)
        return np.sqrt(np.sort(vals)) * scale
    def entropy(nu):
        eps = 1e-12
        S = 0.0
        for v in nu:
            v = max(v, 0.5 + eps)
            S += (v + 0.5) * np.log(v + 0.5) - (v - 0.5) * np.log(v - 0.5)
        return S
    S_i = entropy(symp_from_xp(Vx_i, Vp_i))
    S_j = entropy(symp_from_xp(Vx_j, Vp_j))
    S_ij = entropy(symp_from_xp(Vx_ij, Vp_ij))
    return max(0.0, S_i + S_j - S_ij)


def calibrate_scale(Cx: np.ndarray, Cp: np.ndarray, center: int) -> float:
    Vx = Cx[np.ix_([center], [center])]
    Vp = Cp[np.ix_([center], [center])]
    eig = np.linalg.eigvals(Vx @ Vp).real
    eig = np.clip(eig, 0, None)
    nu = np.sqrt(np.max(eig))
    if nu < 0.5:
        return 0.5 / nu
    return 1.0
