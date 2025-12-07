import numpy as np


def build_K(
    n: int,
    k_coupling: float = 1.0,
    onsite: float = 0.0,
    defect_index: int | None = None,
    defect_mass: float | None = None,
    coupling_boost: float | None = None,
) -> np.ndarray:
    """Stiffness matrix for a 1D chain (open b.c.), with onsite term (gap), optional local mass defect or coupling boost."""
    K = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        diag = onsite
        if i > 0:
            K[i, i - 1] = -k_coupling
            diag += k_coupling
        if i < n - 1:
            K[i, i + 1] = -k_coupling
            diag += k_coupling
        K[i, i] = diag
    if defect_index is not None:
        if defect_mass is not None:
            diag = defect_mass * defect_mass
            if defect_index > 0:
                diag += k_coupling
            if defect_index < n - 1:
                diag += k_coupling
            K[defect_index, defect_index] = diag
        if coupling_boost is not None:
            if defect_index > 0:
                K[defect_index, defect_index - 1] = -coupling_boost
                K[defect_index - 1, defect_index] = -coupling_boost
            if defect_index < n - 1:
                K[defect_index, defect_index + 1] = -coupling_boost
                K[defect_index + 1, defect_index] = -coupling_boost
            if defect_index > 0:
                K[defect_index, defect_index] += coupling_boost - k_coupling
            if defect_index < n - 1:
                K[defect_index, defect_index] += coupling_boost - k_coupling
    return K


def ground_state_covariance(K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return <x x^T> and <p p^T> for ground state of H=1/2(p^T p + x^T K x)."""
    w, v = np.linalg.eigh(K)
    eps = 1e-8
    # Project out zero modes explicitly
    mask = w > eps
    w_nz = w[mask]
    v_nz = v[:, mask]
    w_sqrt = np.sqrt(w_nz)
    w_inv_sqrt = 1.0 / w_sqrt
    # Covariance with ħ=1/2 convention: <x x^T> = Vx, <p p^T> = Vp
    Cx = v_nz @ np.diag(0.5 * w_inv_sqrt) @ v_nz.T
    Cp = v_nz @ np.diag(0.5 * w_sqrt) @ v_nz.T
    return Cx, Cp


def symplectic_eigenvalues(Sigma: np.ndarray) -> np.ndarray:
    """Symplectic eigenvalues via Williamson: positive eigvals of |i Ω Σ|."""
    n = Sigma.shape[0] // 2
    Omega = np.block([[np.zeros((n, n)), np.eye(n)], [-np.eye(n), np.zeros((n, n))]])
    M = Omega @ Sigma
    vals = np.linalg.eigvals(1j * M)
    vals = np.abs(vals)
    vals = np.sort(vals)
    return vals[n:]


def ent_entropy_von_neumann(Sigma: np.ndarray) -> float:
    """Von Neumann entropy for a Gaussian state from its covariance matrix block."""
    nu = symplectic_eigenvalues(Sigma)
    eps = 1e-8
    def f(x: float) -> float:
        xp = max(x, 0.5 + eps)
        return (xp + 0.5) * np.log(xp + 0.5) - (xp - 0.5) * np.log(xp - 0.5)
    S = np.sum([f(v) for v in nu])
    return float(max(S, 0.0))


def subsystem_covariance(Cx: np.ndarray, Cp: np.ndarray, subset: list[int]) -> np.ndarray:
    """Extract covariance matrix for canonical variables of a subset of sites."""
    idx = np.array(subset, dtype=int)
    Cx_sub = Cx[np.ix_(idx, idx)]
    Cp_sub = Cp[np.ix_(idx, idx)]
    Sigma = np.block([[Cx_sub, np.zeros_like(Cx_sub)], [np.zeros_like(Cp_sub), Cp_sub]])
    return Sigma
