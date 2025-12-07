# IA5 – Proxies entropiques vs géométriques (sandbox rapide)

Ce dossier contient une première expérience autonome pour comparer des proxys locaux (géométrie vs entropie) sur un champ 1/f anisotrope, en reprenant l’esprit des tests IA1–IA4 mais avec un script léger exécutable.

## Fichiers
- `entropy_experiments.py` : génère un champ 1/f anisotrope, résout Poisson périodique, calcule des features locales (m, |∇m|², Δm, entropie 3×3×3), entraîne des régressions linéaires et mesure Pearson/relRMSE sur un split train/test.
- `entropy_results.json` : sorties chiffrées (runs par seed et statistiques agrégées).

## Paramètres du run courant
- Grille : 32³ ; α=1.8 ; anisotropie (1.0, 0.8, 0.6) ; seeds = {0,1,2}.
- Modèles : baseline (m), géométrie (m + |∇m|² + Δm), entropie (m + entropie locale 3×3×3 + Δm).
- Split : 60/40 aléatoire sur voxels intérieurs (patch 3×3×3 valide).

## Résultats (moyenne ± écart-type, n=3)
- baseline m : Pearson ≈ 0.018 ± 0.015 ; relRMSE ≈ 1.732 ± 0.008
- géométrie m+|∇m|²+Δm : Pearson ≈ 0.031 ± 0.018 ; relRMSE ≈ 1.730 ± 0.049
- entropie m+H₃×₃×₃+Δm : Pearson ≈ 0.034 ± 0.011 ; relRMSE ≈ 1.676 ± 0.107

## Lecture rapide
- Sur ce setup minimal, aucun proxy n’explique bien y=‖∇Φ‖ (relRMSE ≫1), mais l’entropie locale est légèrement meilleure que la géométrie pure, et toutes dépassent le baseline m.
- Le script offre un cadre reproductible pour itérer (champs fractals/anisotropes, diffusions, autres proxys). Ajuster `alpha`, `anisotropy`, ou ajouter des features/patchs plus larges pour rapprocher les performances des runs IA1–IA4.
