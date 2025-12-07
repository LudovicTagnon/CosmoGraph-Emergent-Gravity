from typing import Tuple

import numpy as np
from numpy.fft import fftn, ifftn


def generate_1overf_field(
    size: int = 32, alpha: float = 1.8, anisotropy: Tuple[float, float, float] = (1.0, 0.8, 0.6), seed: int = 0
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    kx = np.fft.fftfreq(size)[:, None, None]
    ky = np.fft.fftfreq(size)[None, :, None]
    kz = np.fft.fftfreq(size)[None, None, :]
    ax, ay, az = anisotropy
    k2 = (kx / ax) ** 2 + (ky / ay) ** 2 + (kz / az) ** 2
    k2[0, 0, 0] = 1.0  # éviter la division par zéro
    amp = 1.0 / np.power(np.sqrt(k2), alpha / 2.0)
    phases = rng.normal(size=(size, size, size)) + 1j * rng.normal(size=(size, size, size))
    rho = np.real(ifftn(amp * phases))
    rho -= rho.min()
    rho /= rho.max() + 1e-9
    return rho


def poisson_potential(rho: np.ndarray) -> np.ndarray:
    n = rho.shape[0]
    rho_k = fftn(rho)
    kx = np.fft.fftfreq(n)[:, None, None]
    ky = np.fft.fftfreq(n)[None, :, None]
    kz = np.fft.fftfreq(n)[None, None, :]
    k2 = kx * kx + ky * ky + kz * kz
    k2[0, 0, 0] = 1.0
    phi_k = -4 * np.pi * (rho_k - rho_k[0, 0, 0]) / k2
    phi_k[0, 0, 0] = 0.0
    return np.real(ifftn(phi_k))


def grad_norm(phi: np.ndarray) -> np.ndarray:
    gx, gy, gz = np.gradient(phi)
    return np.sqrt(gx * gx + gy * gy + gz * gz)


def laplacian(f: np.ndarray) -> np.ndarray:
    return (
        -6 * f
        + np.roll(f, 1, 0)
        + np.roll(f, -1, 0)
        + np.roll(f, 1, 1)
        + np.roll(f, -1, 1)
        + np.roll(f, 1, 2)
        + np.roll(f, -1, 2)
    )
