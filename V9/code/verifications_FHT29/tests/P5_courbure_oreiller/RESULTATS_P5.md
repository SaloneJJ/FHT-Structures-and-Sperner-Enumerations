# P5 : lecture provisoire de N = 4 (16 nov. 2026, ~21 h 15 UTC), en attendant le contrôle de référence (Wilson L = 32)
Critères : CRITERES_P5.md avec son amendement de 8 h 35. Journal : P5_v2.log ; données : P5_N4.json.

## État du calcul
- `mesure 4 96` est terminé (16 h 18 AST).
- `wilson 5 240` tourne depuis 16 h 18 AST sans rien imprimer. Le message « inversion spectrale » est imprimé après la factorisation LU, donc le calcul est encore dans `splu`, sur une matrice de dimension 4·2·8⁵ = 262 144 (réseau 3D, remplissage important).
  À N = 4 (dimension 32 768), factorisation et diagonalisation prenaient 66 s. À dimension ×8 en 3D, il faut compter plusieurs dizaines de minutes à plusieurs heures, et une mémoire qui peut dépasser la limite de 14 Go.

## Lecture littérale (procédure écrite d'avance : niveaux regroupés à 10⁻⁴ près, comparaison par rang)
- ρ_n Wilson : 1, 1.5423, 1.8448, 1.9928, 2.0573, 2.428.
- ρ_n natif : 1, 1.531, 1.532, 1.822, 1.9649, 2.0243.
- Max |écart| = 17 %.
Le critère amendé exige aussi le contrôle de stabilité de Wilson entre L = 16 et L = 32 ; sans lui, le verdict est « non résolu ».

## Défaut de la procédure, constaté après lecture et déclaré comme tel
- Le regroupement à 10⁻⁴ coupe en deux des niveaux natifs presque dégénérés : 0.33009 et 0.33031 (écart relatif 6.5×10⁻⁴), 0.43646 et 0.43677 (7×10⁻⁴). Chacun de ces couples correspond à un seul niveau de Wilson de multiplicité 16.
- La comparaison par rang aligne donc des niveaux différents. Le rapport des multiplicités le signale : 0.25 au lieu de 0.5 attendu.
- Relecture post hoc (couples fusionnés, comparaison niveau à niveau, multiplicités natif/Wilson toutes égales à 0.5) :

| niveau | Wilson | natif | écart |
|---|---|---|---|
| 1 | 1.5423 | 1.5315 | −0.7 % |
| 2 | 1.8448 | 1.8221 | −1.2 % |
| 3 | 1.9928 | 1.9649 | −1.4 % |
| 4 | 2.0573 | 2.0250 | −1.6 % |
| 5 | 2.4280 | 2.3752 | −2.2 % |

  Max 2.2 % ≤ 5 %, ce qui serait **« standard »**. Ce n'est **pas** un résultat au sens des critères : c'est une relecture après coup, à confirmer par une procédure réécrite avant un nouveau calcul (appariement par multiplicité).
- La trace de la chaleur, qui ne dépend d'aucun appariement, donne un rapport natif/Wilson stable vers 0.518–0.536 pour t de 5 à 60, soit environ 1/2 : une espèce de Weyl contre deux. C'est cohérent avec la lecture « standard ».
- L'ajustement v = 0.81, ξ_eff = 0.83 du journal repose sur l'appariement erroné : il est **non valide**.

## Suite proposée
1. Décider du sort de `wilson 5 240` : le laisser finir si la mémoire tient, sinon l'arrêter.
2. Réécrire avant tout nouveau calcul la procédure d'appariement : fusion des niveaux à 2×10⁻³ près, appariement par multiplicité, ξ_eff recalculé. L'appliquer à N = 4 (données déjà là) et au contrôle de Wilson.

## Application de l'amendement 2 (21 h 30 UTC ; script appariement_v2.py ; journal appariement_v2.log)
- Grappes, seuil 2×10⁻³ : 10 niveaux dans la fenêtre de chaque côté. Les **10 paires sont valides**, toutes avec un rapport de multiplicités natif/Wilson exactement 1/2.
- Écarts relatifs des ρ_n : −0.7, −1.2, −1.4, −1.6, −2.2, −2.5, −2.6, −2.6, −2.7 %. Max 2.7 % ≤ 5 % : **STANDARD, sous réserve du contrôle de Wilson à L = 32.**
- Ajustement sur les 10 paires : v = 0.961, ξ_eff = +0.13.
- Sensibilité au seuil de fusion (10⁻³ et 5×10⁻³) : résultats identiques.
- Constat : l'écart croît régulièrement avec le niveau, toujours dans le même sens. Cela ressemble à une différence de relation de dispersion entre les deux discrétisations plutôt qu'à un effet de courbure, mais seul le contrôle L = 32 peut le dire.
- La trace de la chaleur (rapport stable vers 0.52) va dans le même sens.
- Lecture provisoire, statut [E] partiel : le Dirac natif se couple à la courbure de l'oreiller comme un fermion de Weyl ordinaire, à 3 % près sur les dix premiers niveaux. Le facteur 1.46 du tore n'est pas reproduit sur la racine fermée.

## Contrôle de la référence (amendement 3) : PASSÉ (22 h UTC ; P5_W5_iter.log, P5_W5_iter.json)
- Essai bloquant à N = 4 : 112 valeurs retrouvées à 1.0×10⁻¹⁴ près, VALIDÉ.
- Wilson L = 32 (dimension 262 144, 13 itérations, 427 s, résidu < 10⁻⁷) : mêmes multiplicités que L = 16 sur les 10 premiers niveaux (8, 16, 8, 8, 16, 8, 8, 16, 16, 8).
  Les énergies sont divisées par deux, comme attendu (L doublé).
- Variation des ρ_n entre L = 16 et L = 32 : au plus 0.76 % (seuil 3 %). **Contrôle passé.**
- Écart natif–Wilson réestimé contre la référence à L = 32 : −0.8, −1.0, −1.1, −1.4, −1.4, −2.8, −2.6, −2.4, −2.4 % ; max 2.8 % ≤ 5 %.
- **Verdict P5 : STANDARD, confirmé.** Le Dirac natif se couple à la courbure de l'oreiller comme un fermion de Weyl ordinaire, à 3 % près sur les dix premiers niveaux. Le facteur 1.46 du tore n'est pas reproduit sur la racine fermée.
- Réserves maintenues :
  - la procédure d'appariement a été réécrite après la première lecture de N = 4 (amendement 2, déclaré) ;
  - l'écart résiduel de 1 à 3 %, régulier, est plus grand que la dérive de Wilson. Il peut venir des artefacts propres à la discrétisation native, qui n'existe qu'à N = 4 dans la fenêtre (à N = 3, aucun niveau sous la bande). On ne peut pas le vérifier sans N = 5 natif, hors de portée de la mémoire du Mac.

## Correction après la relecture indépendante bis
Puissance du test : un terme de courbure 1.46 fois celui de Weyl déplacerait les ρ_n d'environ 4 % (estimation sur S³ rond par le relecteur), sous le seuil de 5 %. ξ_eff = +0.13 est du même signe.
La formulation « facteur 1.46 non reproduit » est donc remplacée par « ni confirmé ni exclu » ; c'est d'ailleurs la formulation fixée d'avance dans CRITERES_P5 (« n'est pas confirmé »).
A3 donne maintenant le résultat littéral (17 %) et l'historique des seuils (3 % puis 5 %).
Test à faire : ajouter un terme de courbure à la référence, ou comparer le coefficient a₁ de la chaleur sur l'oreiller.
