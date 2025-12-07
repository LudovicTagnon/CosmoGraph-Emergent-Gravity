# No-Go — Proxies entropie/variance comme médiateur universel

## Constat
- Sur champs 1/f synthétiques (patch ≥7, α≈2), entropie/variance locales améliorent nettement ‖∇Φ‖.
- Sur halos jouet (N-body light), l’effet persiste (relRMSE ~0.72 vs baseline ~0.85).
- Sur cas physiques simples (toy disque baryonique, snapshot Enzo 64³), le signal est nul ou marginal :
  - Disque : gains ±1–3% ou dégradation ; corr ≈ 0.
  - Enzo : corr ≈ 0 ; relRMSE quasi identique (differences ~0.5%).
- Topologie b₀ : signal faible et coûteux ; pas de gain significatif.

## Conclusion
- L’hypothèse “entropie/variance locale explique la gravité manquante de façon générale” est falsifiée dans notre cadre : elle dépend fortement de la structure du champ et ne se transfère pas aux datasets physiques testés.
- Les proxys locaux (entropie, variance, b₀, géométrie simple) ne fournissent pas de médiateur robuste pour la gravité émergente sur des snapshots cosmologiques (Enzo) ni sur un disque baryonique.

### État des jouets “It from Qubit” (entanglement)
- Chaîne 1D d’oscillateurs : déplacement d’un défaut de masse fait varier l’entropie (cloche inversée), confirmant le lien position ↔ intrication. Mais :
  - ΔS vs masse sature (~0.112 max), pas ∝ m.
  - dS/dx décroît vite (power ~ x^-0.8 ou exp ~ e^-x/2.7), portée courte, pas 1/r².
  - Coupling boost augmente ΔS (~0.22) mais même décroissance courte.
- Grille 2D (8×8) : ΔS faible (~0.05 pour défaut mass=4), pas de signal long-portée observé.
→ Force entropique locale observée, mais courte portée et saturante : analogie gravitationnelle non retrouvée dans ces jouets.

## Éléments à conserver / rejeter
- À conserver (ciblé) : usage éventuel d’entropie/variance sur sandbox 1/f ou halos jouet (application locale limitée, pas une théorie générale).
- À rejeter (dans le cadre “universel”) : entropie/variance locale comme signature globale de gravité émergente ; topo b₀ comme médiateur clé.

## Trace des fichiers
- Résultats Enzo : `python_project/outputs/nbody_enzo_results.json`, `nbody_enzo_patch5.json`, `nbody_enzo_patch9.json`.
- Disque baryonique : `python_project/outputs/disk_results.json`.
- Synthèse globale : `python_project/reports/status_overview.md`.
