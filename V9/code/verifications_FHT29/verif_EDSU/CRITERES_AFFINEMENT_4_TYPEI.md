# Affinement 4 — admixture de type I et θ23, δ : critères écrits AVANT tout calcul (16 nov. 2026, ~21 h 15 UTC)

## Attente analytique, écrite avant le calcul
Sous (M3′), tous les blocs neutrinos (Y_ν, M_R, Y_Δ) sont réels dans la base géométrique, donc m_ν = m_I + m_II est réelle.
Le Thm « Reality gives μτ-reflection » d'A2 vaut pour **toute** matrice de masse légère réelle, avec ou sans la transposition résiduelle.
On attend donc :
- une admixture de type I sans (H3) déplace sin²θ12 (colonne TM1, table P6:tab:rg) ;
- elle laisse **exactement** s23 = 0.5 et δ = ±90° ;
- la course avec un Y_ν réel non universel (termes (Y_νᵀ Y_ν) κ + κ (Y_νᵀ Y_ν), réels dans la base géométrique) conserve la réalité, donc la μτ-réflexion ;
- seul Y_e, diagonal dans la base propre de C3 et donc complexe dans la base géométrique, la brise, et son effet a été chiffré à l'affinement 3.

## Vérification numérique
- Textures d'h4_H3_tm1.py (FHT_26), recopiées sans changer la physique :
  - Y_ν C3-invariant réel (αP + β(1 − P) + γJ) ;
  - M_R et Y_Δ symétriques réels invariants par (0 1).
- R = ‖m_II‖/‖m_I‖ ∈ {0.3, 1, 3, 10, 30, 100} ; 300 tirages par R, pour chacun des trois spectres de M_R.
- Variante avec course : κ → (1 + ε S) κ (1 + ε S), avec S = Y_νᵀ Y_ν / ‖Y_νᵀ Y_ν‖ et ε ∈ {0.05, 0.2} (bien au-delà de y_ν² ln(Λ/M_R)/16π² ≈ 0.01).
- Angles calculés dans la base des leptons chargés Q = (v1, v_ω, v_ω̄), avec e = v1 et l'orientation d'A2.
- Critère : max |s23 − 0.5| < 10⁻⁸ et max ||δ| − 90°| < 10⁻⁶ sur tous les tirages, sans course et avec course. On rapporte aussi la médiane de |Δ sin²θ12|.
- Si le critère échoue, l'attente analytique est fausse, et on le dit.

## Conséquence pour A2 (si le critère passe)
Sous (M3′), aucune admixture de type I ni aucune course du secteur neutrino ne peut abaisser θ23. Le coût Δχ² ≈ 4.8 de (M3′) est irréductible.
Seule la brisure de C3 dans le secteur chargé (Georgi–Jarlskog) déplace θ23, et seulement jusqu'à 0.490.
