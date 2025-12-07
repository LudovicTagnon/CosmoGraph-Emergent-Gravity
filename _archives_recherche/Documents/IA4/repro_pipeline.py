# repro_pipeline.py — Pipeline minimal (lecture humaine)
# Voir rapport_complet_gravite_emergente.txt pour les détails.

# Étapes (pseudocode):
# 1) Générer champ 1/f anisotrope (L=32), normaliser [0,1].
# 2) Poisson périodique (Φ_k=-4πG(ρ_k-ρ̄δ_0)/||k||^2 ; Φ_{k=0}=0), y=||∇Φ||.
# 3) Features baseline: [ρ, FD local, moy, var, skew, kurt] sur patch 7^3/9^3.
# 4) Topologie: b0(θ) pour θ∈{0.35,0.50,0.65} sur le patch (6-connexe).
# 5) Ridge λ=0.8, standardisation par train, split 60/40.
# 6) Mesures: Δρ (gain Pearson), ΔrelRMSE (diff. d'erreur relative), IC95% (CLT seeds).
