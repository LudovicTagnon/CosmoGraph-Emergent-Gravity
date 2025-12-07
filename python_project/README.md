# 1D Quantum Gravity Toy Model via Scalar Fields

## Final Results
- Emergent entropic gravity in critical scalar lattices :
  - Static MI decay: dimension-dependent power laws (α≈0.64 in 1D ; α≈1.45–1.52 in 2D/3D). Long-range interaction distinct from Newton (α=2) but geometric.
  - Geometry: MI dip around a massive defect → metric bump g\_xx ∝ 1/MI (see `outputs/gravity_well_heatmap.png`).
  - Dynamics: wavepacket deflection (~0.6 shift) and strong bending/capture attempts (see `dynamics_lensing.png`, `orbit.png`).
- Rapport final : `reports/FINAL_QUANTUM_GRAVITY_REPORT.md`.

## Résumé
Petit univers jouet 1D : chaîne d’oscillateurs couplés en vide critique (gap quasi nul). Un défaut massif modifie l’intrication du vide : l’entropie chute près du cut, la force entropique suit une loi de puissance longue portée, et l’information mutuelle entre voisins s’effondre autour du défaut, signalant une “courbure” émergente (métrique effective g\_xx ∝ 1/MI).

## Méthodologie
- Hamiltonien discret bosonique : H = 1/2 (p^T p + x^T K x), bords ouverts, couplage k, terme onsite faible (près du critique). Défaut = masse locale ou boost de couplage.
- États gaussiens : covariance au sol via diagonalisation de K (projet zéro modes).
- Entropie/forces : entropie de Von Neumann par valeurs propres symplectiques ; force ∝ dS/dx.
- Géométrie : Information Mutuelle I(i,i+1) = S\_i + S\_{i+1} – S\_{ij} (spectre de Vx·Vp, calibration Heisenberg). Déformation métrique : g\_xx(i) ∝ 1/I(i,i+1).

## Résultats clés
- Phase 2 (Force) : vide critique (onsite→0), N=128 → ΔS saturant, force dS/dx décroît en puissance avec exposant α ≈ 0.64 (portée longue, α<1, cohérent 1D critique). Fit exponentiel donne ξ≈8 (longue portée relative).
- Phase 3 (Géométrie) : profil MI (N=64, onsite=1e-6, défaut m=20) : vide ~constante (≈1.6–1.7), avec défaut MI chute jusqu’à ≈0–0.8 au centre. Interprétation : la masse dilate l’espace (g\_xx↑ là où MI↓).
- Phases 4–5 (2D/3D) : MI entre deux défauts décroît en puissance avec α ≈ 1.45 (2D) et α ≈ 1.52 (3D, PBC). Interaction longue portée, dimension-dépendante.
- Phases 6–7 (Dynamique) : lentille (déviation ~0.6) et trajectoire fortement courbée/capture (voir `dynamics_lensing.png`, `orbit.png`). Carte du puits entropique : `gravity_well_heatmap.png`.

### Phase 2 : Gravité Vectorielle (Succès)
- Modèle d’élasticité 2D (champ vectoriel ux, uy) : avec un défaut de masse modéré (m=3) la MI décroît en loi de puissance avec α ≈ 2.04, compatible avec 1/r².
- Commandes :
  - Scaling : `python3 scripts/run_vector_scaling.py --L 16 --m 1e-6 --defect-mass 3`
  - Dynamique (déviation) : `python3 scripts/run_vector_dynamics.py`
  - Orbite (capture) : `python3 scripts/run_vector_orbit.py`
- Artefacts : `outputs/vector_scaling_defect_L16_m3.json`, `outputs/vector_scaling_defect_L16_m3.png`, `outputs/vector_dynamics_orbit.png`, `outputs/vector_orbit_M30_P0.50.png`.

## Reproduire
1) Prérequis : Python 3, venv conseillé.
```bash
cd python_project
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```
2) Force (Phase 2) :
```bash
python3 scripts/run_qit_scaling.py --n 128 --onsite 0.0 --defect-masses 2 4 8 16 32 --output qit_scaling_critical_N128.json
```
Résultats : `outputs/qit_scaling_critical_N128.json` (ΔS, forces, fits via post-process).
3) Géométrie (Phase 3) :
```bash
python3 scripts/run_qit_geometry.py --n 64 --onsite 1e-6 --defect-mass 20 --output qit_geometry.json
```
Profil MI et métrique effective : `outputs/qit_geometry.json`, `outputs/qit_geometry.png`.

## Arborescence (simplifiée)
- `src/gravite_py/qit.py` : chaîne 1D, covariance, entropie, forces.
- `scripts/run_qit_scaling.py` : Phase 2 (force vs distance et masse, vide critique).
- `scripts/run_qit_geometry.py` : Phase 3 (information mutuelle, métrique).
- `scripts/run_2d_gravity.py`, `scripts/run_3d_gravity.py` : Phases 4–5 (MI entre défauts en 2D/3D).
- `scripts/run_dynamics_lensing.py`, `scripts/run_orbit.py` : Phases 6–7 (dynamique, lentille, capture).
- `scripts/map_gravity_well.py` : Heatmap du puits (MI vs centre).
- `outputs/` : JSON/PNG des runs (scaling, géométrie, dynamique).
- `reports/FINAL_QUANTUM_GRAVITY_REPORT.md` : synthèse finale.
