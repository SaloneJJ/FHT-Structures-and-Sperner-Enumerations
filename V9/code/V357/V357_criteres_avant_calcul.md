# V357 — critères du relecteur écrits AVANT ses calculs (13/11, heure du cloud = AST, Guadeloupe)

R341 (V341, borne ε > 1 − 1/√3) — analyse papier du relecteur : dans le secteur γ_F = s, F est conservé ; la masse de profondeur (coefficient de Γ5)
vaut m_K − μ pour le survivant (s f = −1) et m_K + μ pour le jumeau (s f = +1), m_K = √3(ε − 1). Chaîne de Wilson–Kaplan (r = 1) :
une copie porte un mode de bord ssi sa masse est dans (−2r, 0). Survivant : ε > 1 − 1/√3. Jumeau : m_K ∈ (−3, −1), soit ε ∈ (1 − √3, 1 − 1/√3).
Donc pour ε < 0,423 (et > −0,73) il y a encore UN Weyl, porté par l'autre copie (goût opposé), de même chiralité (V331 §4), de décroissance
λ_j = |2 − √3(1 − ε)|. La borne de V341 est une transition de goût du survivant, pas une borne d'existence du Weyl.
Prédictions (arbre N = 3, q = 1, 64 pas, mode unif, μ = r = t_c = 1, V336_FHT_flot.py tel quel) :
  * ε = 0,2 (λ_j = 0,614, comparable à ε = 0,6 de V341 : λ′ = 0,69 → −0,55) : F(3) ∈ [−0,85 ; −0,30], même signe que V341 (négatif).
  * ε = 0,1 (λ_j = 0,441, comparable à ε ≈ 0,75) : F(3) ∈ [−0,90 ; −0,40], négatif.
  Si |F(3)| < 0,15 pour les deux : mon analyse est fausse sur l'arbre (doubleurs latéraux ou autre) et la borne de V341 tient.
R340 (V340) — analyse papier : pour un motif η uniforme par niveau, la chiralité nette des modes nuls d'une région [k, N] contenant le bord
vaut 1 (somme télescopique sur les deux goûts) pour tout η ; le Weyl est déplacé (vers la dernière paroi), pas détruit ; les parois
intérieures portent des paires de chiralités opposées (goûts opposés). Contrôle : chaîne réduite (Kaplan en profondeur), N = 3 et N = 10.
R345 (V345 K1) : premier ordre exact : m(h) = ε_D |ψ(h)|² = ε_D (1 − λ′²) λ′^{2h} ; chaîne indépendante de 40 couches : écart < 1 % pour ε_D = 1e-3.
R342 (V342 M1, M2) : algèbre recalculée avec mes propres matrices ; M3 : modèle « tranche » (dalle de profondeur uniforme, Bloch latéral) :
  gap des survivants = ε_D (AY) et 0 au premier ordre (Γ5).
R352 (V352) : contrôles indépendants : (i) n_eff → 2 et d_s → 3 sur la dalle ; (ii) Z_ins/Z_0 = e^{−tε_D²} à 1e-3 près pour t ≥ 5 ;
  (iii) MS à une boucle recalculé ; (iv) relance de V352_noyau.py à N = 2 (m = 8) pour n_eff(t = 20) ∈ [1,9 ; 2,3] (après retrait de la contrepartie).
