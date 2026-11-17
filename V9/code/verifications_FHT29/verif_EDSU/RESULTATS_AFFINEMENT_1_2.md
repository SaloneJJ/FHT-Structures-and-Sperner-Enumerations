# Résultats des affinements 1 et 2 (critères : CRITERES_AFFINEMENT_1_2.md ; script affinement_1_2.py ; journal affinement_1_2.log)
16 nov. 2026, ~20 h 05 UTC.
- Contrôle de non-régression GJ (s13 = 0.0220) : OK, bornes d'A2 retrouvées (0.2688–0.3740 ; 0.4933–0.5097 ; 68.5°–111.6°).
- Diagnostic : l'image CP-conjuguée est déjà atteinte par le balayage en α (distance 0).
- Changement de code déclaré, sans effet sur la physique : scipy est absent de la VM, la p-valeur du χ² (1–3 ddl) et l'inversion gaussienne ont été codées à la main (formules fermées, bissection).

## Le contrôle B (T12/DCP 2D + T23 1D) diffère de plus de 1 de la construction A dans deux cas : signalé, les deux valeurs sont données
Les deux constructions sont des approximations du Δχ² à trois paramètres, que NuFIT ne publie pas. On retient l'intervalle [B, A].

## TByes, ordre normal (variante principale)
| modèle | χ² (A ; B) | ddl | p | équiv. | verdict |
|---|---|---|---|---|---|
| TM1 exact (s23 = 0.5, δ = −90°) | 8.96 ; 6.67 | 3 | 0.030 ; 0.083 | 2.2σ ; 1.7σ | en tension / au bord |
| TM1 seul (sans (M3′)) | 4.20 ; 4.30 | 2 | 0.12 | 1.5σ | compatible |
| GJ, vrai ajustement | 5.58 ; 4.06 | 2 | 0.061 ; 0.13 | 1.9σ ; 1.5σ | compatible (au bord) |
- Coût de (M3′) : 4.8 (1 ddl, p = 0.03). Le rapport χ²(TM1 seul)/χ²(TM1 exact) vaut 0.47 < 0.8.
  **Selon le critère écrit d'avance, la tension n'est PAS attribuée à TM1 : environ la moitié vient de (M3′)**, c'est-à-dire de θ23 = 45° exact.
  Mon attente avant le calcul (« la tension vient de TM1 ») est donc **réfutée** ; je l'avais dite dans le chat.
- Prédiction de TM1 seul : s23 ∈ [0.453, 0.477] et δ ∈ [−101.8°, −95.7°] à Δχ² ≤ min + 1 ; s23 ∈ [0.443, 0.495] et δ ∈ [−104.5°, −91.2°] à + 4. Seul δ < 0 est retenu.
- GJ : meilleur ajustement α = 42.5°, orientation échangée, en (s12, s23, δ) = (0.310, 0.490, −89°). Le minimum est **au bord** de la courbe (s23 minimal) : les données tirent θ23 plus bas que ce que GJ permet.
  - Prédiction à + 1 : s12 ∈ [0.303, 0.316], s23 ∈ [0.490, 0.494], δ ∈ [−97.8°, −86.1°].
  - À + 4 : s12 ∈ [0.296, 0.322], s23 ∈ [0.490, 0.510], δ ∈ [−100.5°, −81.7°].
- Sensibilité à s13 (±1σ) : moins de 0.2 sur tous les χ².

## Autres variantes
- TBoff NO : TM1 exact 10.2 ; 7.7 (2.4σ ; 1.9σ). TM1 seul 5.1 (1.8σ). GJ 6.8 ; 5.0 (2.1σ ; 1.7σ). Coût de (M3′) : 5.0.
- IO (relatif au minimum IO, TByes) : TM1 exact 4.8 (1.3σ), TM1 seul 2.0, GJ 2.5. Mais l'IO est lui-même à Δχ² = 5.9.
  En IO, TM1 seul préfère s23 = 0.55 et δ = −77°.

## Lecture
1. TM1 avec (P0) et le marquage n'est pas en difficulté (p = 0.12) ; sa prédiction propre est δ ≈ −99° ± 3° avec θ23 ≈ 42.9°–43.6°, soit la relation de somme TM1 connue (R).
2. La réalité des blocs neutrinos (M3′) impose θ23 = 45° (μτ-réflexion), ce qui coûte Δχ² ≈ 4.8. GJ réduit ce coût en écartant θ23 jusqu'à 0.490, sans pouvoir aller plus bas.
3. Test décisif : DUNE et Hyper-K sur θ23. Si sin²θ23 < 0.485 est établi, (M3′) est réfuté dans les deux scénarios, TM1 seul pouvant survivre ; si θ23 ≈ 45° est confirmé, l'énoncé fort de (M3′) l'emporte.

## Correction après la relecture indépendante bis
Le coût de (M3′) vaut 4.77 avec la construction A et 2.37–2.43 avec la construction B. Dans A2, on écrit donc « entre un tiers et la moitié (Δχ² = 2.4–4.8) ».
GJ sur la plage complète des deux orientations, s23 ∈ [0.490, 0.510] : 5.07 (TByes), 6.26 (TBoff), au lieu de 5.56 et 6.75.
