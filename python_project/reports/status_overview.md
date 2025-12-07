# État des lieux — Proxies et tests (IA5)

## Synthèse rapide
- Champs 1/f anisotropes (synthétiques) : entropie/variance locales sur patchs larges (≥7) apportent des gains nets sur ‖∇Φ‖ (relRMSE ~0.81–0.96 selon α/aniso), supérieurs à la géométrie seule ; topo b₀ faible.
- Diffusion (iso/guidée) : l’avantage entropie/variance persiste à t=0→2.
- Toy halos (N-body light maison) : gain significatif (relRMSE ~0.72 vs baseline ~0.85).
- Toy disque baryonique : gains nuls/faibles (relRMSE ≈ baseline ±1–3%).
- Snapshot Enzo (64³) extrait : gains quasi nuls (~0.5%), corr ~0.

Conclusion : le médiateur entropie/variance n’est pas robuste sur des cas physiques simples (disque, Enzo). Topo b₀ reste coûteux et peu performant. La partie “entropie locale explique la gravité manquante” est au mieux contextuelle (1/f, halos jouet), non validée en général.

## Fichiers et scripts clés
- Scripts : `run_entropy.py`, `run_sweep.py`, `run_diffusion.py`, `run_disk.py`, `run_nbody_grid.py`, `run_batch.py`.
- Sorties : `python_project/outputs/` (notamment `nbody_enzo_results.json`, `nbody_light_results.json`, `disk_results.json`, `sweep_results.json`, etc.).
- Données Enzo : archive `Enzo_64.tar.gz` extraite dans `/tmp/enzo_tmp` (utilisée pour `nbody_enzo_rho.npz`).

## Résultats chiffrés (exemples)
- 1/f, patch 9, α=2.0, aniso 1/0.6/0.4 : entropie+variance relRMSE ~0.81 vs baseline ~0.82.
- Diffusion guidée, patch 7, α=1.8 : entropie/variance relRMSE ~0.96 vs baseline ~0.98 à t=0→2.
- N-body light (halos jouet, grille 64³) : entropie/variance relRMSE ~0.72 vs baseline ~0.85.
- Enzo 64³ (snapshot) : entropie/variance relRMSE ~0.0605 vs baseline ~0.0607 (gain négligeable ; corr ~0).
- Disque baryonique : gains nuls ou négatifs (relRMSE ~0.23–0.25 selon patch/aniso).

## Prochaines décisions
- Acter l’échec sur cas physiques (disque, Enzo) : partie “entropie locale universelle” non validée.
- Si persistance : tester d’autres snapshots (ou niveaux Enzo) pour confirmer le no-signal.
- Si pivot : redéfinir les hypothèses (ou explorer la piste Gémini “It from Qubit” avec de nouveaux observables), documenter les no-go.
