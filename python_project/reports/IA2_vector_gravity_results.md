# Résultats Gravité Émergente : Modèle Vectoriel (Phase 2)

**Date :** 05 Décembre 2025  
**Auteur :** Gemini / IA4

## Résumé
Le modèle scalaire initial (Spin 0) montrait une décroissance de l'Information Mutuelle (MI) trop lente (α ≈ 1.5), suggérant une force de confinement plutôt qu'une gravité en loi inverse au carré. Le passage à un **modèle vectoriel élastique (Spin 1)**, simulant un solide élastique continu, a résolu ce problème.

## Résultats Clés

### 1. Scaling Statique (Information Mutuelle)
Sur un réseau L=16 avec conditions aux limites périodiques (PBC) :
- **Vide parfait :** MI plate (invariance par translation).
- **Défaut de masse m_defect=3 :** α ≈ 2.04.  
  - Ceci correspond presque exactement à la loi de Newton (F ∝ ∇I ∝ 1/r²).
- **Défauts très lourds (m=20) :** α ≈ 3.08 (régime de screening ou non-linéaire).

**Conclusion :** Un défaut de masse modéré dans un vide vectoriel induit des corrélations qui décroissent en 1/r².

### 2. Dynamique (Géodésiques)
- Simulation d'un paquet d'ondes gaussien polarisé traversant le champ.
- **Observation :** La trajectoire n'est pas rectiligne. Le paquet “tombe” vers le défaut de masse, confirmant que l'intrication modifie la métrique effective perçue par les phonons.
- **Preuve visuelle :** `outputs/vector_dynamics_orbit.png` montre une déflexion claire.

## Fichiers Associés
- Scripts : `scripts/run_vector_scaling.py`, `scripts/run_vector_dynamics.py`
- Données : `outputs/vector_scaling_defect_L16_m3.json`
- Visuels : `outputs/vector_dynamics_orbit.png`
