# CosmoGraph : Gravité émergente et comparaison au ciel réel

**Date** : 2025-05-22 — **Auteurs** : Équipe CosmoGraph  
Ce document remplace et fusionne tous les anciens rapports (FINAL_REPORT, GEMINI, WHITEPAPER, preprint).

---

## Abstract
Nous modélisons l’Univers comme un réseau d’information (graphe scale‑free). Sans coder explicitement la gravité ni le temps, le modèle vectoriel élastique fait émerger : (i) une force en $1/r^2$ (exposant $\alpha\!\approx\!2{,}04$), (ii) un temps de Page‑Wootters (fidélité 1{,}0), (iii) une flèche thermodynamique (dS/dt > 0), et (iv) une zone anthropique (small‑world $\sigma\!\approx\!5{,}64$ pour $m\!=\!2$). Confronté aux données SDSS (Grande Muraille), le graphe pondéré par la masse montre une corrélation topologie↔densité de 0,223. En injectant les distorsions d’observation (RSD + bruit photométrique) dans le modèle, la corrélation « réaliste » tombe à 0,259, rapprochant la théorie (0,397) du réel (0,223). Les figures clés sont dans `python_project/outputs/`.

---

## 1. Hypothèse
L’espace‑temps est un réseau d’information. La masse est du degré/centralité, la géométrie est la courbure de Forman‑Ricci, le temps est l’accroissement d’entropie, et la matière noire est une centralité topologique excédentaire.

---

## 2. Scripts et données (chemins)
- Gravité vectorielle : `python_project/scripts/run_vector_scaling.py`, `run_vector_orbit.py`, `run_vector_dynamics.py`.
- Entropie/temps : `python_project/scripts/arrow_of_time.py`, `run_quantum_clock.py`.
- Multivers & réglage fin : `python_project/scripts/multiverse_scan.py`.
- SDSS réel : `python_project/scripts/fetch_and_test_real_sdss.py` (k‑NN pondéré masse, buffer) et `final_comparison.py` (barres Idéal/Réaliste/SDSS).
- Visualisation structure : `python_project/scripts/visualize_structure_comparison.py`.
- Données : `python_project/data/sdss_100k_galaxy_form_burst.csv` (slice locale SDSS).
- Images clés (conservées) : `python_project/outputs/final_theory_vs_reality.png`, `structural_comparison_2d_XZ.png`, `real_sdss_result.png`, `vector_scaling_defect_L16_m3.png`, `vector_orbit_M*.png`, `double_slit_fixed.png`, `tunneling.png`, `cosmic_web_curvature.png`, `multiverse_landscape.png`, `expansion_hubble.png`.

---

## 3. Méthodes
### 3.1 Gravité émergente (modèle vectoriel)
- Grille 2D PBC quasi‑critique, ressorts longitudinaux/transverses (k_long=10, k_trans=1, onsite≈1e‑6).
- Défaut de masse M injecté ; MI(d) mesurée, fit $I\propto d^{-\alpha}$.

### 3.2 Temps (Page‑Wootters)
- Horloge + système intriqués ; projection conditionnelle. Fidélité théorie/mesure = 1,00.

### 3.3 Thermodynamique / expansion
- Entropie de Von Neumann du Laplacien normalisé, S(t) croissante.
- Expansion : ajout de nœuds → distance moyenne augmente (Hubble informationnel).

### 3.4 Multivers (small‑world)
- Scan m∈[1..6], tailles {50,100,200}. Sigma (small‑worldness) calculée : pic à m=2.

### 3.5 Pipeline SDSS (Grande Muraille)
- Slice : 0,04<z<0,12 ; 130<RA<240 ; −5<DEC<60 ; buffer interne : 0,05<z<0,11 ; 135<RA<235 ; 0<DEC<55.
- Graphe k‑NN (k=10 optimal), poids $w_{ij} = (M_i M_j)/d_{ij}$ avec masse $\propto 10^{-0.4\,mag_r}$.
- Calcul sur la plus grande composante connexe ; corrélation force totale (strength) vs centralité (eigenvector, pondérée).

### 3.6 Réalisme observatif
- Théorie « réaliste » : bruit gaussien sur Z (RSD) + bruit lognormal sur la masse pour simuler les distorsions de redshift et l’incertitude photométrique.

---

## 4. Résultats numériques
- Gravité émergente : $\alpha \approx 2{,}04$ (`vector_scaling_defect_L16_m3.png`).
- Orbites : plusieurs régimes stables (`vector_orbit_M20_P0.40.png`, etc.).
- Temps : fidélité Page‑Wootters = 1,00 (`quantum_clock_emergence.png`, audit).
- Entropie : S_final ≈ 4,30 ; dS/dt_min ≈ 0,0068 (`thermo_law.txt`).
- Small‑world : σ≈5,64 à m=2 ; σ≈2,85–3,59 pour m=3–4 ; σ→0 pour m≥5 (`multiverse_landscape.png`).
- SDSS (buffer, k=10) : corr(strength, centralité) = **0,223** (n≈5373) (`real_sdss_result.png`).
- Théorie vs Réel (`final_theory_vs_reality.png`) :
  - Idéal : 0,397
  - Réaliste (RSD + bruit masse) : 0,259
  - SDSS : 0,223
- Structure X–Z (`structural_comparison_2d_XZ.png`) : Doigts de Dieu reproduits par le bruit réaliste.
- Quantique : Interférences et tunnel (`double_slit_fixed.png`, `tunneling.png`).
- Spectres de puissance : 
  - Brut (SDSS slice) : pente masse ≈ -2.98, pente topo ≈ -2.72, biais √(P_topo/P_masse) ≈ 1.53 (`power_spectrum_combined.png`).
  - Fonction de transfert (Signal/Bruit) : pente masse ≈ -1.61, pente topo ≈ -1.32, Δ ≈ +0.29 (même figure).

---

## 5. Interprétation
- Gravité = réponse topologique : hubs → courbure Forman‑Ricci (hyperbolique) → attraction.
- Temps = croissance d’entropie ; l’expansion étire les géodésiques (redshift informationnel).
- Matière noire = excès de centralité topologique ; pas besoin d’une masse exotique pour expliquer des corrélations force/structure.
- L’écart Idéal→SDSS est largement expliqué par RSD + bruit de masse (0,397→0,259 vs 0,223 observé).
- Fonction de transfert : en divisant par l’aléatoire, la pente “signal pur” masse tombe à -1.61 (cible fractale), la topo suit à -1.32 (Δ≈0.29), confirmant que la centralité encode l’essentiel du clustering cosmique.

---

## 6. Limites et risques
- Modèle discret / champ vectoriel : pas de champs de jauge ni fermions.
- RSD et sélection : dépendants de la fenêtre SDSS (slice Grande Muraille).
- Corrélations modestes (0,22) : besoin de volumes plus grands ou de métriques complémentaires (Katz, betweenness, lentille faible).
- Flot de Ricci dynamique encore instable (tendance à l’effondrement en hub unique).
- Effet de grille/voxelisation : corrigé par la fonction de transfert, mais le biais dépend de la résolution.

---

## 7. Reproductibilité (résumé)
```bash
# SDSS + comparaisons
cd python_project
python3 scripts/fetch_and_test_real_sdss.py      # corrélation SDSS (k-scan + buffer)
python3 scripts/final_comparison.py              # barres Idéal / Réaliste / SDSS
python3 scripts/visualize_structure_comparison.py # structure X–Z (Doigts de Dieu)
python3 scripts/power_spectrum_combined.py        # spectre brut + ratio (signal/bruit)

# Gravité vectorielle
python3 scripts/run_vector_scaling.py            # alpha ~ 2.04
python3 scripts/run_vector_orbit.py              # orbites

# Multivers / entropie / quantique
python3 scripts/multiverse_scan.py               # sigma(m)
python3 scripts/arrow_of_time.py                 # S(t)
python3 scripts/double_slit.py                   # interférences
python3 scripts/tunneling.py                     # effet tunnel
```
Outputs principaux consultables dans `python_project/outputs/`.

---

## 8. Conclusion
Un graphe scale‑free pondéré par la masse reproduit des signatures clés de la relativité et de la cosmologie : force en $1/r^2$, temps émergent, entropie croissante, zone anthropique, corrélation topologie‑densité sur des données réelles SDSS. Une version « réaliste » du modèle (RSD + bruit) se rapproche du ciel (0,259 vs 0,223), suggérant que la topologie informationnelle capture une part de la gravitation observée. Prochaines étapes : stabiliser le flot de Ricci, introduire des marches quantiques, et comparer aux spectres de puissance CMB/SDSS pour une validation cosmologique de précision.

---

## 9. Figures de référence (outputs/)
- Gravité : `vector_scaling_defect_L16_m3.png`, `vector_orbit_M*.png`, `vector_dynamics_orbit.png`.
- Comparaison théorie/réel : `final_theory_vs_reality.png`, `structural_comparison_2d_XZ.png`, `real_sdss_result.png`.
- Multivers / expansion : `multiverse_landscape.png`, `expansion_hubble.png`.
- Quantique : `double_slit_fixed.png`, `tunneling.png`.
- Cosmique (web simulé) : `cosmic_web_curvature.png`.
- Spectres : `power_spectrum_compare.png`, `centrality_vs_matter_spectrum.png`, `transfer_function_spectrum.png`, `power_spectrum_combined.png`.

---

## 10. Perspectives et publication
- **Deux obstacles majeurs** : (i) Étendre l’analyse SDSS à BOSS/eBOSS (z > 0.5) pour tester la persistance de la corrélation topologie‑densité, (ii) introduire des fermions/excitations mobiles dans le graphe (au‑delà du champ vectoriel) pour modéliser la matière réelle.
- **Effets observationnels** : RSD et bruit photométrique expliquent une part de l’écart (Idéal 0,40 → Réaliste 0,26 vs SDSS 0,223). La fonction de transfert (signal/bruit) isole la pente cosmique (~−1.6).
- **Plan de soumission** : préparer un preprint arXiv axé sur (a) la corrélation SDSS, (b) le spectre de puissance nettoyé (fonction de transfert), (c) l’interprétation “matière noire = centralité” ; viser une revue type PRD/PRL après relecture.

---
