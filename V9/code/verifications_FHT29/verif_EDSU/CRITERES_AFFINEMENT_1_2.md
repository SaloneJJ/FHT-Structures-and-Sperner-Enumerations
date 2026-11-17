# Affinements 1 et 2 — critères écrits AVANT tout calcul (16 nov. 2026, ~19 h 55 UTC)
Données : tables NuFIT 6.1 (nufit61/, identifiées et contrôlées par sha256). Variante principale : TByes (IC24 + SK-atm) ; TBoff en contrôle. Ordre principal : NO ; IO rapporté.

## Construction du χ² (approximation déclarée)
- χ²(s12, s23, δ) ≈ Δχ²₂D[T23/DCP](s23, δ) + Δχ²₁D[T12](s12) − m, où m est le minimum de l'ordre (pour que la somme reparte de 0 au minimum de l'ordre). On rapporte aussi la valeur relative au minimum global.
- Hypothèse : θ12 est presque décorrélé de (θ23, δ). Contrôle : refaire le calcul avec T12/DCP (2D) + T23 (1D). Si les deux évaluations diffèrent de plus de 1 en Δχ² au point retenu, on le signale et on donne les deux.
- s13 fixé à la valeur de NuFIT 6.1 (0.02248 NO, 0.02262 IO). Sensibilité : ±1σ (0.00057), rapportée sans être minimisée.
- Interpolation bilinéaire sur les grilles, sans lissage.

## Affinement 1 — coût de TM1 et coût de (M3′)
- Modèle « TM1 seul » ((P0) + (M1)–(M2), sans (M3′)) :
  - sin²θ12 = 1 − 2/(3(1 − s13²)) ;
  - cos δ = −(1 − 5 s13²) cot 2θ23 / (2√2 s13 √(1 − 3 s13²)) (Thm TM1 d'A2), les deux signes de δ étant permis ;
  - paramètre libre : s23 ∈ [0.35, 0.65], pas de 0.0005, sur la partie de la courbe où |cos δ| ≤ 1.
- Modèle « TM1 exact » (+ (M3′)) : s23 = 0.5, δ = ±90°, sans paramètre libre.
- Résultats demandés :
  - χ²_min de TM1 seul, à 2 ddl (3 observables, 1 paramètre), avec la p-valeur et l'équivalent gaussien bilatéral ;
  - χ² de TM1 exact, à 3 ddl ;
  - coût de (M3′) = χ²(TM1 exact) − χ²_min(TM1 seul), à 1 ddl.
- Prédiction de TM1 seul : intervalles de s23 et de δ sur la courbe pour Δχ² ≤ χ²_min + 1 et + 4.
- Lecture prévue : si χ²_min(TM1 seul) ≥ 0.8 × χ²(TM1 exact), la tension est attribuée à TM1, donc à (P0) avec le marquage, et non à (M3′).

## Affinement 2 — vrai ajustement de Georgi–Jarlskog
- Courbe : script V358_B_bande (texture κ = 1, r = −3, masses V356 inchangées), recopiée sans changer la physique. On y remplace s13 = 0.0220 par 0.02248 (changement d'entrée déclaré).
  α ∈ [0°, 360°) au pas de 0.5° ; toutes les solutions en θ (une ou deux branches).
- Orientation (lequel de v_ω, v_ω̄ est le muon) : choix binaire libre (A2). On l'implémente exactement, en échangeant les lignes μ et τ de U, puis on recalcule s23 et δ. Les deux images sont incluses.
  Diagnostic : on vérifie si l'image CP-conjuguée (δ → −δ) est déjà atteinte par le balayage en α.
- Paramètre libre : α (et la branche). χ²_min à 2 ddl, avec la p-valeur et l'équivalent gaussien.
- Prédiction : pour les points de la courbe avec χ² ≤ χ²_min + 1 et + 4 (1 paramètre), les intervalles de s12, s23 et δ, à transmettre à DUNE et Hyper-K.
- Contrôle de non-régression : avec s13 = 0.0220, le balayage doit redonner les bornes publiées dans A2 (s12 ∈ [0.269, 0.374], s23 ∈ [0.493, 0.510], |δ| ∈ [68.5°, 111.6°]) à 0.002 et 0.5° près. Sinon, on s'arrête.

## Verdicts (identiques pour les trois modèles)
p > 0.05 : compatible ; 0.0027 < p ≤ 0.05 : en tension ; p ≤ 0.0027 : exclu.
A2 est mis à jour avec les nombres obtenus, quel qu'en soit le sens.
