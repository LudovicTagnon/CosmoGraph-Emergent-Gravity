# Méga-rapport (IA1 → IA5) — Gravité émergente

## Panorama chronologique
- IA1 (FMM, champs gaussiens) : hypothèse contenu → géométrie (K̃ = −Δ log c, K_Regge) → cinématique. Gains nets sur temps d’eikonale T : ΔrelRMSE ≈ −31 % pour K̃ (patch 10³, L=18), ≈ −6–7 % pour Regge (Lc=8, marge ≥2). Gradient |∇m| rivalise avec K̃ ; fusion m+K̃+Regge > m seul. Split damier train/test, seeds multiples, robustesse t=0→2 (back-reaction rapide).
- IA2 (Annexes) : confirmations sur variantes et ablations (gradient vs K̃, marge Regge, bootstrap blocs), ΔC-index > 0, relRMSE ↓. K_Regge seul modeste ; K̃ et fusion dominent. Placebos (phase scramble) annulent le gain. Marge intérieure Regge ≥2 stabilise.
- IA3 (médiation statistique) : t=0 validé avec instruments forts (F>1900, ΔC-index +0.01–0.02), mais échec à t>0 avec instruments stricts (p_J rejette). Signal sensible au médiateur et à la dynamique ; médiateur dynamique ou ΔK/ΔT non conclusifs.
- IA4 (topologie b₀ sur champs 1/f anisotropes L=32) : gains modestes mais significatifs (α=2.0 patch 7³ : Δρ≈+0.0817, ΔrelRMSE≈−0.0065 ; α=1.5 patch 9³ similaire). Preuves formelles : monotonie E[‖∇Φ‖|b₀=k], irréductibilité vs descripteurs lisses. Sandbox N-body : b₀ quasi constant/inconsistant → piste “b₀ révolutionnaire” falsifiée ; forme/concentration explique mieux (compatible Poisson).
- IA5 (refonte Python, proxys entropie/variance, diffusion, N-body light, Enzo) : pipeline modulable, tests exhaustifs, no-go formalisés.

## Synthèse par proxy / contexte
### Géométrie (K̃, Regge, grad²)
- Champs gaussiens (IA1/IA2) : K̃ fort (ΔrelRMSE ≈ −31 %), Regge modéré (≈ −6–7 %), gradient proche de K̃. Fusion > m seul. Gains persistent t=0→2 sous diffusion rapide.
- Champs 1/f (IA5) : géométrie seule faible (relRMSE ~0.97–1.73 selon patch/aniso), inférieure à entropie/variance.
- Enzo, disque : géométrie n’apporte rien ou marginal.

### Topologie b₀
- IA4 : gains modestes sur 1/f (Δρ ~ +0.08 ; ΔrelRMSE ~ −0.0065). Preuves mathématiques (monotonie, irréductibilité) dans le cadre 1/f.
- N-body sandbox : b₀ quasi constant → signal nul. IA5 : b₀ coûteux, relRMSE ~0.97 (patch topo 7) < entropie/variance.
- Conclusion : non robuste, falsifié comme médiateur clé.

### Entropie / variance locales
- 1/f anisotrope (IA5) : patch ≥7, α≈2 → relRMSE ~0.81–0.96 (entropie/variance) vs baseline ~0.82–0.98. Diffusion (iso/guidée) : avantage persistant t=0→2.
- N-body light (halos jouet sur grille 64³) : relRMSE ~0.72 (ent/entvar) vs baseline ~0.85 ; Pearson ↑ (~0.63 vs 0.42).
- Disque baryonique toy : gains nuls/faibles (±1–3 %, corr ≈0) ; parfois dégradation.
- Snapshot Enzo 64³ : corr ~0 ; relRMSE quasi identique (gains ~0.5 % ou nuls, patch 5/7/9). Confirmé sur plusieurs patchs.
- Conclusion : médiateur efficace dans des sandboxes 1/f/halos jouet, falsifié sur cas physiques (disque, Enzo).

### Diffusion / temporalité
- IA1/IA2 : back-reaction rapide t=0→2 conserve l’avantage K̃/fusion.
- IA5 diffusion iso/guidée : entropie/variance gardent un léger avantage (relRMSE ~0.96 vs ~0.98 baseline) sur 1/f.
- Médiation IA3 : signal disparaît t>0 sous instruments stricts.

## Falsifications actées
- Entropie/variance locale comme médiateur universel de gravité émergente : non validé (échec sur Enzo et disque).
- Topologie b₀ comme clé révolutionnaire : falsifiée (N-body, tests IA5).
- Proxys géométriques simples seuls : insuffisants hors champs lisses/gaussiens.

## Ce qui reste exploitable (local/limité)
- K̃/Regge/grad sur champs gaussiens FMM (preuve de principe locale).
- Entropie/variance sur sandboxes 1/f ou halos jouet (application contextuelle, pas générale).

## Fichiers et repères
- IA1/IA2 : `IA1/IA1 SAVE.txt`, `Documents drive/Annexes Rapport IA2 v2/*.csv`, résultats K̃/Regge/grad, bootstrap.
- IA3 : `IA3/Synthèse IA3.txt`, annexe PRC, instruments et médiation.
- IA4 : `IA4/rapport_complet_gravite_emergente.txt`, `IA4/rapport_exhaustif_gravite_emergente.txt`, `GEMINI.md`, `IA4/Conversation complète.txt`.
- IA5 code : `python_project/` (scripts, src, outputs). Rapports : `python_project/reports/status_overview.md`, `python_project/reports/no_signal_entropy.md`.
- Résultats IA5 : `python_project/outputs/` (sweeps, diffusion, disk, nbody_light, Enzo).
- Données Enzo extraites : `python_project/outputs/nbody_enzo_rho.npz`; résultats : `nbody_enzo_results.json`, `nbody_enzo_patch5.json`, `nbody_enzo_patch9.json`.
- Méga-rapport (ce fichier) à partager pour onboarding.

## Enseignements / recherche par élimination
- Les proxies locaux (entropie/variance/b₀/géométrie simple) ne donnent pas de signature robuste sur un snapshot cosmologique ni sur un disque baryonique. L’hypothèse “observable local universel” est à rejeter.
- Le signal dépend de la structure statistique : présent sur 1/f/halos jouet, absent sur Enzo → pas une nouvelle loi de gravité, plutôt un ajustement contextuel.
- Les approches topologiques naïves (b₀) ne survivent pas aux tests réalistes.

## Pistes prochaines (au-delà du no-go)
1) Documenter le no-signal (déjà dans `reports/no_signal_entropy.md`) et clore la piste entropie/variance universelle.
2) Explorer un cadre alternatif : “It from Qubit” (Théorie de Gémini) → définir des observables d’intrication/cohérence/renormalisation plutôt que stats locales ; tester sur données synthétiques et snapshots.
3) Si besoin de retests : autres snapshots (ou niveaux Enzo) pour confirmer l’absence de signal ; sinon pivot complet.
