# Résultats de l'affinement 3 (critères : CRITERES_AFFINEMENT_3_RG.md ; script affinement_3_rg.py ; journal affinement_3_rg.log)
16 nov. 2026, ~21 h UTC.
- Point de départ à Λ = 10¹⁶ GeV : TM1 exact (s12 = 0.31800, s23 = 0.5, δ = −90°), phases de Majorana 0 ou π. Ordre normal, m1 ∈ [0, 4.25] meV, quatre choix de signes.
- Décalages maximaux de Λ à M_Z (une boucle, modèle standard, base de Y_e diagonale) :
  - facteur 1 (ε_τ = 3.2×10⁻⁵) : |Δs12| ≤ 1.8×10⁻⁵, |Δs23| ≤ 2.2×10⁻⁵, |Δδ| ≤ 0.002° ;
  - marge ×2, deux signes : 3.5×10⁻⁵, 4.4×10⁻⁵, 0.004° ;
  - marge ×10 : 1.8×10⁻⁴, 2.2×10⁻⁴, 0.022°.
- **Verdict : RIGIDE**, avec un facteur d'environ 50 sous les seuils fixés (0.001 ; 0.002 ; 1°).
- Grandeur complémentaire : en masses quasi dégénérées, |Δs23| vaut 0.026 à m1 = 1 eV et reste sous 0.035. Il faudrait donc m1 > 1 eV, exclu par KATRIN et par la cosmologie.
  La course du modèle standard ne peut pas réconcilier θ23 = 45° avec le sin²θ23 ≈ 0.465 préféré par les données.
- Limites déclarées :
  - une boucle ;
  - y_τ constant (conservatif) ;
  - au-dessus des seuils, Y_ν universel (H3) ;
  - le cas Y_ν non universel (admixture de type I) est l'affinement 4.
  - Petite erreur d'affichage dans le journal (Σmν ≈ 3 m1 imprimé aussi pour m1 = 0.01 eV, où ce n'est pas valable) : sans effet sur le verdict.

## Correction après la relecture indépendante bis (16/11, ~23 h UTC)
Le verdict « rigide » ne vaut que **sous M_Δ** (modèle standard). Entre Λ et M_Δ, la boucle du triplet renormalise L par un terme en Y_Δ†Y_Δ (Schmidt, PRD 76, 073010, 2007), avec Y_Δ = √2 g m̂ d'ordre 1.
Ce terme entre dans la course de Y_e, fait tourner la base des leptons chargés hors de la base de C3 et brise la μτ-réflexion.
Estimation du relecteur : |Δ sin²θ23| ~ 10⁻², de signe non fixé. Elle n'est pas encore calculée ici.
Les critères de l'affinement 3 avaient omis ce terme : erreur de cadrage, reconnue. A2 est corrigé : « irréductible au niveau de l'arbre et sous M_Δ ».
