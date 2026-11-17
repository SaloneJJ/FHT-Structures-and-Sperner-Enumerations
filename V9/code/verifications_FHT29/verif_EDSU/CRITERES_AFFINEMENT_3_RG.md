# Affinement 3 — transport par le groupe de renormalisation : critères écrits AVANT tout calcul (16 nov. 2026, ~20 h 55 UTC)

## Cadre et justification
- Les relations TM1 exact et μτ-réflexion (s23 = 0.5, δ = ±90°, phases de Majorana 0 ou π) sont posées à l'échelle Λ = 10¹⁶ GeV (celle du triplet d'A2). Les données sont à M_Z.
- Une boucle. Sous les seuils, opérateur effectif κ (Weinberg) du modèle standard :
  16π² dκ/dt = C[(Y_e†Y_e)ᵀ κ + κ (Y_e†Y_e)] + α κ, avec C = −3/2 (Antusch et al. 2001).
  - Le terme α κ est universel : il change les masses, pas les angles. Il est omis pour les angles.
  - Dans la base où Y_e est diagonal, et en gardant y_τ constant, la solution est exacte : κ(M_Z) = D κ(Λ) D, avec D = diag(1, 1 + ε_μ, 1 + ε_τ) et ε_ℓ = C y_ℓ² ln(M_Z/Λ)/(16π²).
- Au-dessus des seuils (triplet M_Δ, neutrinos droits M_R) :
  - les termes en Y_Δ seul conservent les vecteurs propres de Y_Δ ;
  - sous (H3), Y_ν ∝ 1 et les termes en Y_ν†Y_ν sont universels ;
  - les seuls termes qui changent les angles sont en Y_e, avec un coefficient d'ordre 1 différent de C.
  On les couvre par une marge : on refait le calcul avec |C| multiplié par 2 et par 10, et dans les deux signes.
  Le cas sans (H3) (Y_ν non universel, admixture de type I) relève de l'affinement 4 et n'est pas traité ici.
- Valeurs : y_τ = 0.0102, y_μ = 6.1×10⁻⁴ (valeurs à M_Z, prises constantes ; elles décroissent vers le haut, donc c'est conservatif).
  Δm² et s13 de NuFIT 6.1 (NO), imposés à basse énergie : les valeurs à Λ sont ajustées par itération.

## Balayage
- Ordre normal. m1 ∈ [0, 4.25] meV (borne DESI DR2 en ΛCDM), 41 valeurs.
- Phases de Majorana : les quatre choix de signes (η2, η3), avec η1 = +1.
- δ(Λ) = −90° (le signe choisi par les données).
- Sorties : Δsin²θ12, Δsin²θ23, Δδ entre Λ et M_Z, et leurs maxima en valeur absolue.

## Verdict
- **Rigide** si, pour tout le balayage et avec la marge ×2 : |Δsin²θ23| < 0.002, |Δδ| < 1° et |Δsin²θ12| < 0.001.
  On l'écrit alors dans A2 : les prédictions sont transportées sans changement mesurable, et la course ne peut pas réconcilier θ23 = 45° avec sin²θ23 ≈ 0.465.
- **Non rigide** sinon : on donne les décalages et leur dépendance en m1.
- Grandeur complémentaire : le m1 (avec des masses quasi dégénérées, sans la borne DESI) pour lequel |Δsin²θ23| atteint 0.035, l'écart entre 0.5 et le meilleur ajustement de TM1 seul. On la compare à la borne cosmologique.
