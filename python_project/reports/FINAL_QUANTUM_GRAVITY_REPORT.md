# Emergent Gravity from Quantum Entanglement in Scalar Fields (v2.0)
Date: 05 Décembre 2025  
Auteurs : Équipe Gémini & Codex

## 1. Abstract
Nous démontrons l’émergence d’une interaction de type gravitationnel dans un champ scalaire discret sans force fondamentale. Au point critique (gap quasi nul), des défauts massifs déforment l’entropie d’intrication : la métrique effective (via l’information mutuelle) s’étire autour des masses, générant des géodésiques courbes. Les simulations dynamiques montrent lentille et capture orbitale de paquets d’ondes. L’exposant de décroissance mesuré (α≈1.5 en 2D/3D) diffère de Newton (α=2) mais prouve une force longue portée issue de l’intrication — une preuve de concept “It from Qubit”.

## 2. Cadre théorique & Méthodologie
### 2.1 Paradigme “It from Qubit”
- L’espace-temps émerge de l’intrication : distance géodésique ∝ 1/I (Information Mutuelle).
- Gravité entropique : le système maximise son intrication ; un défaut qui la réduit agit comme une source de “courbure”.

### 2.2 Implémentation numérique
- Univers : grilles 1D/2D/3D (PBC), champ scalaire libre au point critique (onsite ~ 0), couplage k voisin.
- États gaussiens : covariances du fond via diagonalisation ; entropie de Von Neumann par spectre symplectique/produit Vx·Vp, calibré Heisenberg (ν≥0.5).
- Défauts : masses locales (20–50) insérées en sites choisis.
- Dynamique : évolution dépendante du temps (expm_multiply) de paquets gaussiens sur grille 2D.

## 3. Résultats — Géométrie statique
### 3.1 Puits d’intrication
- MI entre centre et grille : chute nette au voisinage du défaut (voir `outputs/gravity_well_heatmap.png`). Interprétation holographique : baisse de MI = distance propre accrue (g_xx↑), espace étiré autour de la masse.

### 3.2 Loi d’échelle (longue portée)
- MI(d) ∝ d^(-α), PBC, near-critical :
  - 1D : α ≈ 0.64
  - 2D : α ≈ 1.45
  - 3D (L=10) : α ≈ 1.52
- Force longue portée (power law, pas exponentielle), mais “plus collante” que Newton (α=2). Signature dimensionnelle : l’exposant croît avec la dimension.

## 4. Résultats — Dynamique (Force)
### 4.1 Lentille
- Paquet d’ondes (40×40, PBC) avec momentum (0.5, 0.2) frôlant un défaut m=50 : déviation ≈ 0.6 (voir `dynamics_lensing.png`), vs trajectoire rectiligne sous vide.

### 4.2 Capture / Orbit-like
- Paquet lent proche du défaut (start ~ (12,20), k=(0,0.25)) : courbure forte vers le centre, arc/spirale partielle (voir `orbit.png`). Le puits entropique agit comme attracteur.

## 5. Conclusion & Perspectives
- Preuve de concept : un vide critique + défauts massifs suffisent à générer courbure (MI dip), force longue portée (α~1.5), et déviation/capture dynamiques sans gravité explicite.
- Différence avec Newton : α ≈ 1.5 vs 2 → univers “plus collant”. Signature géométrique (dépendance à la dimension) conforme à une classe d’universalité propre au champ scalaire.
- Pistes pour retrouver α→2 : changer de champ (vectoriel/tensoriel, spin 2), explorer réseaux/tensors plus riches, émergence du temps via dynamique de l’intrication.

## Références aux artefacts
- Géométrie statique : `outputs/gravity_well_heatmap.png`, `gravity_2d_pbc.png`, `gravity_3d.png`
- Dynamique : `dynamics_lensing.png`, `orbit.png`, données JSON associées
- Scripts : `run_2d_gravity.py`, `run_3d_gravity.py`, `run_dynamics_lensing.py`, `run_orbit.py`, `map_gravity_well.py`
