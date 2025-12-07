import os
import numpy as np
import matplotlib.pyplot as plt

def fake_lcdm(k: np.ndarray, slope=-2.9, norm=1e5):
    # Simple power law placeholder; replace with CAMB/CLASS if available
    return norm * k ** slope

def load_spectrum(path):
    data = np.load(path, allow_pickle=True).item() if path.endswith('.npy') else None
    raise FileNotFoundError

# For simplicity we reuse existing power_spectrum_combined plot by overlaying a reference curve

def main():
    # Retrieve existing Pk from png is not feasible; so we regenerate from saved arrays if present.
    # Here we assume power_spectrum_combined.py saved arrays for mass/topo; if not, skip.
    arrays_path = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'power_arrays.npz')
    if not os.path.exists(arrays_path):
        print("No saved arrays (power_arrays.npz). Skipping LCDM overlay.")
        return
    data = np.load(arrays_path)
    k_m = data['k_m']; Pm = data['Pm']; k_t = data['k_t']; Pt = data['Pt']

    # LCDM reference
    k_ref = k_m
    P_lcdm = fake_lcdm(k_ref)

    plt.figure(figsize=(8,5))
    plt.loglog(k_m, Pm, 'o-', label='Mass weighting', color='orange')
    plt.loglog(k_t, Pt, 'x--', label='Topo weighting', color='blue')
    plt.loglog(k_ref, P_lcdm, ':', label='LCDM ref (power-law)', color='gray')
    plt.xlabel('k [a.u.]'); plt.ylabel('P(k)')
    plt.title('Power spectra with LCDM reference (placeholder)')
    plt.legend(); plt.grid(True, which='both', alpha=0.3)
    out = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'power_with_lcdm.png')
    plt.tight_layout(); plt.savefig(out, dpi=200)
    print(f"Saved {out}")

if __name__ == '__main__':
    main()
