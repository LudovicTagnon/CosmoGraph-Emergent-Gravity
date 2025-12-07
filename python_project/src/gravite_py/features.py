from typing import Dict, Tuple, Iterable, List

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .fields import grad_norm, laplacian, poisson_potential


def local_entropy(field: np.ndarray, patch_size: int = 3, bins: int = 16) -> np.ndarray:
    """Entropie de Shannon sur patch patch_size³ (recadrée au centre)."""
    if patch_size % 2 == 0 or patch_size < 3:
        raise ValueError("patch_size doit être un entier impair ≥ 3")
    scaled = np.clip((field - field.min()) / (field.max() - field.min() + 1e-9), 0, 1)
    window = sliding_window_view(scaled, (patch_size, patch_size, patch_size))
    patches = window.reshape(-1, patch_size**3)
    digitized = np.clip((patches * bins).astype(int), 0, bins - 1)
    ent = np.empty(digitized.shape[0], dtype=np.float32)
    for i, row in enumerate(digitized):
        counts = np.bincount(row, minlength=bins)
        p = counts / counts.sum()
        mask = p > 0
        ent[i] = -np.sum(p[mask] * np.log(p[mask]))
    n = field.shape[0]
    out_shape = n - patch_size + 1
    return ent.reshape(out_shape, out_shape, out_shape)


def make_dataset(rho: np.ndarray, patch_size: int = 3, bins: int = 16) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Construit les features baseline, géométriques et entropiques + cible y=||∇Φ||."""
    if patch_size % 2 == 0 or patch_size < 3:
        raise ValueError("patch_size doit être un entier impair ≥ 3")
    n = rho.shape[0]
    phi = poisson_potential(rho)
    y_full = grad_norm(phi)
    offset = patch_size // 2
    valid = (slice(offset, n - offset), slice(offset, n - offset), slice(offset, n - offset))
    m = rho[valid]
    y = y_full[valid]
    grad2 = sum(np.gradient(rho)[i][valid] ** 2 for i in range(3))
    lap = laplacian(rho)[valid]
    ent = local_entropy(rho, patch_size=patch_size, bins=bins)
    var = local_variance(rho, patch_size=patch_size)
    X_baseline = np.stack([m.ravel()], axis=1)
    X_geom = np.stack([m.ravel(), grad2.ravel(), lap.ravel()], axis=1)
    X_ent = np.stack([m.ravel(), ent.ravel(), lap.ravel()], axis=1)
    X_var = np.stack([m.ravel(), var.ravel(), lap.ravel()], axis=1)
    X_ent_var = np.stack([m.ravel(), ent.ravel(), var.ravel(), lap.ravel()], axis=1)
    return {"baseline_m": X_baseline, "geom_m_grad2_lap": X_geom, "ent_m_entropy_lap": X_ent, "var_m_var_lap": X_var, "entvar_m_entropy_var_lap": X_ent_var}, y.ravel()


def local_variance(field: np.ndarray, patch_size: int = 3) -> np.ndarray:
    """Variance locale sur patch patch_size³ (recadrée au centre)."""
    if patch_size % 2 == 0 or patch_size < 3:
        raise ValueError("patch_size doit être un entier impair ≥ 3")
    window = sliding_window_view(field, (patch_size, patch_size, patch_size))
    patches = window.reshape(-1, patch_size**3)
    var = np.var(patches, axis=1)
    n = field.shape[0]
    out_shape = n - patch_size + 1
    return var.reshape(out_shape, out_shape, out_shape)


def build_feature_sets(rho: np.ndarray, patch_size: int = 3, bins: int = 16) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    """Retourne un dict de matrices X par modèle et la cible y aplatie."""
    features, y = make_dataset(rho, patch_size=patch_size, bins=bins)
    return features, y


def count_components(binary_patch: np.ndarray) -> int:
    """Compte les composantes connexes 6-voisins dans un patch binaire 3D."""
    visited = np.zeros(binary_patch.shape, dtype=bool)
    comps = 0
    dims = binary_patch.shape
    neighbors = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    for x, y, z in np.argwhere(binary_patch):
        if visited[x, y, z]:
            continue
        comps += 1
        stack: List[Tuple[int, int, int]] = [(x, y, z)]
        visited[x, y, z] = True
        while stack:
            cx, cy, cz = stack.pop()
            for dx, dy, dz in neighbors:
                nx, ny, nz = cx + dx, cy + dy, cz + dz
                if 0 <= nx < dims[0] and 0 <= ny < dims[1] and 0 <= nz < dims[2]:
                    if binary_patch[nx, ny, nz] and not visited[nx, ny, nz]:
                        visited[nx, ny, nz] = True
                        stack.append((nx, ny, nz))
    return comps


def local_b0(rho: np.ndarray, patch_size: int = 7, thresholds: Iterable[float] = (0.35, 0.5, 0.65)) -> np.ndarray:
    """Betti0 (nb de composantes) par patch et par seuil."""
    if patch_size % 2 == 0 or patch_size < 3:
        raise ValueError("patch_size doit être un entier impair ≥ 3")
    windows = sliding_window_view(rho, (patch_size, patch_size, patch_size))
    patches = windows.reshape(-1, patch_size, patch_size, patch_size)
    thresholds = tuple(thresholds)
    b0_vals = []
    for th in thresholds:
        vals = np.empty(len(patches), dtype=np.float32)
        for i, patch in enumerate(patches):
            vals[i] = count_components(patch > th)
        b0_vals.append(vals)
    b0_stack = np.stack(b0_vals, axis=1)  # [num_patches, num_thresholds]
    out_shape = rho.shape[0] - patch_size + 1
    return b0_stack.reshape(out_shape, out_shape, out_shape, len(thresholds))


def build_feature_sets_with_topology(
    rho: np.ndarray, patch_size: int = 3, bins: int = 16, topo_patch: int | None = None, topo_thresholds: Iterable[float] = (0.35, 0.5, 0.65)
) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    """Retourne un dict de matrices X incluant b0 multi-seuils (patch topo distinct possible)."""
    features, y = make_dataset(rho, patch_size=patch_size, bins=bins)
    tp = topo_patch or max(patch_size, 7)
    b0 = local_b0(rho, patch_size=tp, thresholds=topo_thresholds)
    # recadrage pour aligner m/lap et b0
    offset = (rho.shape[0] - b0.shape[0]) // 2
    slice_b0 = (slice(offset, offset + b0.shape[0]),) * 3
    m_center = rho[slice_b0].ravel()
    lap_center = laplacian(rho)[slice_b0].ravel()
    b0_flat = b0.reshape(-1, b0.shape[-1])
    X_topo = np.concatenate([m_center[:, None], b0_flat, lap_center[:, None]], axis=1)
    features["topo_m_b0multi_lap"] = X_topo
    return features, y
