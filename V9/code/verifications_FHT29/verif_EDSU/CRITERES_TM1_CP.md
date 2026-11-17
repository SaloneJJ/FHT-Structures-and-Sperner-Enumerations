# Énoncé « TM1 ⇒ violation de CP forte » : critères écrits AVANT le calcul (16 nov. 2026, ~21 h 40 UTC)
- Relation de somme TM1 (Thm TM1 d'A2) : cos δ = −(1 − 5 s²) cot 2θ23 / (2√2 s √(1 − 3 s²)), avec s = sin θ13.
- Entrées : sin²θ13 dans [0.02077, 0.02419] (NuFIT 6.1 NO ±3σ, avec 0.00057 par σ) ; sin²θ23 dans son intervalle à 3σ, NO [0.435, 0.585] (IO [0.440, 0.584] en contrôle).
- Sorties :
  - l'intervalle de |δ| sur ce domaine ;
  - les valeurs de sin²θ23 qui donneraient δ = 0 et δ = 180°, au s13 central et aux bornes.
- Énoncé à mettre dans A2 si : |δ| ∈ [a, b] avec a > 60° et b < 120° sur tout le domaine, et sin²θ23(δ = 0 ou 180°) hors de [0.435, 0.585] pour tout s13 du domaine.
  Sinon, on écrit les nombres tels quels, sans l'énoncé « CP conservée exclue ».
- Contrôle croisé : on recalcule |δ| directement sur des matrices TM1 explicites (U_TBM R23(θ, φ)), par l'invariant de Jarlskog, sur une grille de (θ, φ).
