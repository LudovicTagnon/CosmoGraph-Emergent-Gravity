import numpy as np


def build_K2d(
    nx: int,
    ny: int,
    k_coupling: float = 1.0,
    mass: float = 1.0,
    defect: tuple[int, int] | None = None,
    defect_mass: float | None = None,
    periodic: bool = False,
) -> np.ndarray:
    """Stiffness matrix for a 2D grid of oscillators with optional defect. Supports PBC if periodic=True."""
    n = nx * ny
    K = np.zeros((n, n), dtype=np.float64)

    def idx(i, j):
        return i * ny + j

    for i in range(nx):
        for j in range(ny):
            p = idx(i, j)
            diag = mass * mass
            neighbors = []
            if periodic:
                neighbors = [
                    ((i + 1) % nx, j),
                    ((i - 1) % nx, j),
                    (i, (j + 1) % ny),
                    (i, (j - 1) % ny),
                ]
            else:
                if i > 0:
                    neighbors.append((i - 1, j))
                if i < nx - 1:
                    neighbors.append((i + 1, j))
                if j > 0:
                    neighbors.append((i, j - 1))
                if j < ny - 1:
                    neighbors.append((i, j + 1))
            for (ii, jj) in neighbors:
                q = idx(ii, jj)
                K[p, q] = -k_coupling
                diag += k_coupling
            K[p, p] = diag

    if defect is not None and defect_mass is not None:
        di, dj = defect
        p = idx(di, dj)
        diag = defect_mass * defect_mass
        if periodic or di > 0:
            diag += k_coupling
        if periodic or di < nx - 1:
            diag += k_coupling
        if periodic or dj > 0:
            diag += k_coupling
        if periodic or dj < ny - 1:
            diag += k_coupling
        K[p, p] = diag
    return K


def ground_state_covariance(K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    w, v = np.linalg.eigh(K)
    w_sqrt = np.sqrt(w)
    w_inv_sqrt = 1.0 / w_sqrt
    Cx = v @ np.diag(0.5 * w_inv_sqrt) @ v.T
    Cp = v @ np.diag(0.5 * w_sqrt) @ v.T
    return Cx, Cp


def subsystem_covariance(Cx: np.ndarray, Cp: np.ndarray, subset: list[tuple[int, int]], nx: int, ny: int) -> np.ndarray:
    idx = np.array([i * ny + j for (i, j) in subset], dtype=int)
    Cx_sub = Cx[np.ix_(idx, idx)]
    Cp_sub = Cp[np.ix_(idx, idx)]
    Sigma = np.block([[Cx_sub, np.zeros_like(Cx_sub)], [np.zeros_like(Cp_sub), Cp_sub]])
    return Sigma


def symplectic_eigenvalues(Sigma: np.ndarray) -> np.ndarray:
    n = Sigma.shape[0] // 2
    Omega = np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])
    eigvals = np.linalg.eigvals(1j * Omega @ Sigma)
    vals = np.sort(np.abs(eigvals))
    return vals[n:]


def ent_entropy_von_neumann(Sigma: np.ndarray) -> float:
    nu = symplectic_eigenvalues(Sigma)
    eps = 1e-12
    def f(x: float) -> float:
        xp = max(x, 0.5 + eps)
        term1 = xp + 0.5
        term2 = xp - 0.5
        return term1 * np.log(term1) - term2 * np.log(term2)
    return float(np.sum([f(v) for v in nu]))
