# T1 — Résultats (Fable, 16/11/2026)

Critères : `CRITERES_T1_T2.md`. Script : `T1_isospin_octet.py` (numpy seul, quelques secondes ; sortie dans `T1_isospin_octet.log`). Conventions de `V335_FHT_C3prime_flot.py` (`Xg`, `Yg`, `G5`, `Am`, `Bm`, `G4`, `FG`).

## Partie algébrique : tous les critères (A1)–(A5) sont satisfaits

| critère | résultat |
|---|---|
| (A1) | le commutant de `{Γ_5, α_d}` est de dimension complexe 4 ; `{1, iΓ_4A, iΓ_4B, AB}` sont quatre matrices réelles, de carrés `−1`, anticommutantes : une algèbre de quaternions `H`, forme réelle fixée par la conjugaison des fonctions sur les huit sites |
| (A2) | `γ = iY_1Y_2Y_3` (égal à `−iΓ_4Γ_5`) est de carré 1, commute aux `α_d` et à `H`, anticommute à `Γ_5` et `Γ_4` |
| (A3) | `P_L H ⊕ P_R H` est de dimension 8 ; `P_L H ⊕ P_R {1, AB}` de dimension 6 ; sur `P_R`, `AB` est une structure complexe et `F = iAB` vaut `+1` et `−1` avec multiplicité 2 chacun |
| (A4) | à `k = K`, `t_f = t_c` : `|E| = m + μ` et `|m − μ|`, quatre fois chacun (vérifié pour trois couples) |
| (A5) | les trois commutateurs valent 1,2·10⁻³, 4,8·10⁻³, 1,9·10⁻², 7,6·10⁻² pour `q` ≈ 0,05, 0,1, 0,2, 0,4 : ordre `q²` |

Lecture : au point `K`, l'octet est un spineur de Dirac tensorisé par un doublet `E`. L'algèbre gauche–droite `H_L ⊕ H_R` y est native (chiralité ⊗ quaternions de saveur), et la partie isospin de `A_F`, `C ⊕ H`, s'y plonge : `H` sur les gauchers, `C = {1, AB}` sur les droitiers, où `F` distingue le haut du bas. Le terme `(m + μF)Γ_5` de A3 (Lemme 9.3(iv)) est une masse `m + μ` pour l'un et `m − μ` pour l'autre : `m` est la masse moyenne du doublet, `μ` la demi-différence haut/bas.

## Lecture d'un journal existant (post hoc)

`FHT_26/03_code/CalculCodes/V335_logs.txt`, mode `bbhG`, `N = 3`, `γ_F = +1`, 64 pas :

| `m_tot` | `μ` | flot de la couche fine |
|---|---|---|
| −1,0 | 0,0 | **−1,844** |
| 0,0 | 1,0 | −0,943 |

Rapport 1,96. Sans masse de saveur, le secteur porte deux branches de Weyl de même chiralité ; c'est ce que la Prop 9.2 de A3 appelle les jumeaux. Une seule exécution à `μ = 0` figure dans les journaux que j'ai lus ; le secteur `γ_F = −1` et `N = 4` manquent.

## Compte

* Lecture actuelle de A3 : 96 états de `F` × 2 (`E`) = 192 branches ; Majorana : 96 ; masse de saveur posée : 48.
* Lecture testée : 2 chiralités × (1 + 3) × 3 générations = 24 secteurs ; chacun porte un doublet : **48 fermions de Weyl**, sans masse de saveur, sans les Choix 9.x qui servent à lever le doublement.

## Ce qui reste à faire, et les commandes

(M1), `N = 3`, environ 5 minutes par exécution sur le Mac. Le module `rvA_common.py` est dans `V331_relecture/` ; adaptez le `PYTHONPATH` si vous lanciez le script autrement.

```
cd ~/Desktop/2026FHT/FHT_26/03_code/CalculCodes
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=V331_relecture ../../../venv/bin/python V335_FHT_C3prime_flot.py 3 -1.0 0.0 -1 64 1 bbhG | tee "../../../FHT_29 Final/00_RETOURS_FABLE/tests/T1_M1_N3_mu0_sm1.log"
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=V331_relecture ../../../venv/bin/python V335_FHT_C3prime_flot.py 3 -1.0 0.0 1 64 1 bbhG | tee "../../../FHT_29 Final/00_RETOURS_FABLE/tests/T1_M1_N3_mu0_sp1.log"
```

(M2), `N = 4` (dimension 37 448 ; durée longue, à lancer quand le Mac est libre) : mêmes commandes avec `4` à la place de `3`, journaux `T1_M2_N4_mu0_sm1.log` et `T1_M2_N4_mu0_sp1.log`.

Je n'ai pas pu exécuter ces commandes (pas de scipy dans mon environnement) ; je ne garantis pas la ligne de commande exacte, que je déduis de l'en-tête du script et des journaux.

## Réserves

1. L'isospin gauche `P_L H` n'est une symétrie qu'à l'ordre dominant près du cône ; les termes de Wilson anisotropes le brisent à l'ordre `q²` (A5). C'est la brisure de goût des fermions décalés ; jauger une symétrie ainsi brisée demande un argument de limite continue.
2. La lecture change ce que A3 tensorise (`F'` de dimension 24 au lieu de `H_F`) et donc le Dirac fini : les Yukawa haut et bas deviennent `(m ± μ)`, c'est-à-dire que `μ` devrait dépendre du secteur lepton/couleur et de la génération. Cela n'est pas testé.
3. Le groupe ponctuel agit sur `E` : l'isospin se transforme sous les rotations de la tuile. C'est cohérent avec le Thm 4.3 de A3 (le groupe ponctuel est une jauge), pas avec une lecture où il serait une rotation d'espace.

## (M1) — exécutions du 16/11 sur le Mac : critère satisfait

Journaux `T1_M1_N3_mu0_sm1.log` et `T1_M1_N3_mu0_sp1.log` (mode `bbhG`, `N = 3`, `m_tot = −1`, `μ = 0`, 64 pas, environ 140 s chacun).

| secteur | flot net par couche | couche fine | critère |
|---|---|---|---|
| `γ_F = −1` | −0,049 ; −1,091 ; −0,544 ; **+1,844** | +1,844 | dans [+1,65 ; +2,05] : **satisfait** |
| `γ_F = +1` | +0,049 ; +1,091 ; +0,544 ; **−1,844** | −1,844 | reproduit le journal `V335_logs.txt` au chiffre près |

Lecture. Sans masse de saveur, chaque secteur de chiralité porte deux branches de Weyl dans la couche fine (1,844 contre 0,943 avec `μ = 1` : rapport 1,96), de signes opposés entre les deux secteurs.

Deux limites à dire.
1. Les deux exécutions sont exactement opposées, à tous les chiffres : `γ_F ↦ −γ_F` est une symétrie exacte du code. Le secteur `−1` n'apporte donc pas d'information indépendante ; (M1) confirme la reproductibilité et le signe, pas davantage.
2. Aucune exécution ne peut, à elle seule, départager « goût » et « isospin » : c'est le même opérateur. Ce que les nombres établissent est le compte (deux branches par secteur à `μ = 0`). Ce qui porte la lecture est la partie algébrique : `H_L ⊕ H_R` et `C ⊕ H` natifs, et `(m + μF)Γ_5` comme masse haut/bas.

Reste (M2), `N = 4`, pour vérifier que le facteur 2 ne dépend pas de la taille.

## Prédiction supplémentaire, écrite avant tout calcul (M3)

Dans la lecture testée, les deux membres du doublet ont les masses de paroi `m_tot + μ` et `m_tot − μ` ; chacun ne porte une branche que si sa masse est du côté topologique. À `m_tot = −1` fixé, le flot de la couche fine doit donc former un palier proche de 1,84 pour `|μ| < 1`, puis tomber vers 0,9 quand `|μ|` dépasse 1. Commandes (`N = 3`, environ 2 minutes chacune) : les mêmes que (M1) avec `0.5` puis `1.5` à la place de `0.0` (troisième argument), journaux `T1_M3_N3_mu05.log` et `T1_M3_N3_mu15.log`.
Critère : `|flot| ∈ [1,65 ; 2,05]` à `μ = 0,5` ; `|flot| ∈ [0,7 ; 1,2]` à `μ = 1,5`. Échec utile : un flot qui varie continûment avec `μ` contredirait l'image « deux membres, deux masses ».

## Correction du 16/11 : (M1) existait déjà — recherche dans FHT_26, FHT_27, FHT_28

JJ avait raison : l'exécution `N = 3`, `m_tot = −1`, `μ = 0` avait déjà été faite, trois fois, dans `FHT_26/03_code/CalculCodes` :
* `V335_logs.txt` et `logs_c3/log_m-1_mu0_s1.txt` (même résultat, −1,844) ;
* `V339_relecture/log_N3_64_alt_1_m-1_mu0_s1_cos.txt`, relecture indépendante avec quatre fenêtres.
(M1) n'a donc rien apporté de neuf ; je l'ai demandé à tort, faute d'avoir cherché au-delà de `V335_logs.txt`.

La relecture `V339` apporte en revanche une réserve que je n'avais pas : le facteur 2 dépend de la fenêtre.

| fenêtre | `μ = 0`, `m_tot = −1` | `μ = 1`, `m_tot = 0` | rapport |
|---|---|---|---|
| box | −1,844 | −0,944 | 1,95 |
| box/2 | −1,507 | −0,908 | 1,66 |
| box×2 | −1,846 | −0,942 | 1,96 |
| gauss | −1,695 | −0,875 | 1,94 |

Le rapport va de 1,66 à 1,96. « Deux branches par secteur » tient pour trois fenêtres sur quatre ; la demi-fenêtre donne 1,66. C'est à dire tel quel.

Inventaire des exécutions de flot trouvées dans FHT_26 à FHT_28 (aucun journal de flot dans FHT_27 ni FHT_28) :
* `μ = 0` : seulement `N = 3`, `m_tot = −1`, mode `bbhG` ;
* `N = 4` : seulement `m_tot = 0`, `μ = 1` (sept exécutions) et `μ = 0,8` (une), mode `unif` ;
* `m_tot = −1` avec `μ ≠ 0` : aucune.
Donc (M2) (`N = 4`, `μ = 0`) et (M3) (`m_tot = −1`, `μ = 0,5` et `1,5`) n'existent pas dans ces trois dossiers, sauf sous un autre nom ou dans un format que ma recherche n'a pas reconnu.

## Commandes données le 16/11 pour (M3) et (M2)

* (M3), `N = 3`, `bbhG`, `m_tot = −1`, `γ_F = +1`, `μ = 0,5` puis `1,5` : script `V335_FHT_C3prime_flot.py` (solveur dense), environ 2 min chacune.
* (M2), `N = 4`, `bbhG`, `m_tot = −1`, `μ = 0`, `γ_F = +1` : script `V338_FHT_flot.py` (même Hamiltonien, solveur lobpcg ; c'est celui des exécutions `N = 4` existantes), environ 2 h.
* Référence pour (M2) : l'exécution existante `N = 4`, `μ = 1` en mode `unif` (couche fine −0,904). À `N = 3`, les modes `unif` et `bbhG` donnent le même flot à `μ = 1`, couche par couche ; une exécution `bbhG` à `N = 4`, `μ = 1` est facultative.
Critères inchangés (écrits avant) : (M3) `|flot| ∈ [1,65 ; 2,05]` à `μ = 0,5`, `[0,7 ; 1,2]` à `μ = 1,5` ; (M2) `|flot(μ = 0)| / 0,904 ∈ [1,8 ; 2,2]`.

## Lecture de (M3) — Claude, 16/11/2026, ≈ 7 h 50 (Guadeloupe), selon le critère écrit avant par Fable
Journaux `T1_M3_N3_mu05.log` et `T1_M3_N3_mu15.log` (exécutés par JJ, `bbhG`, `N = 3`, `m_tot = −1`, `γ_F = +1`, 64 pas) :

| `μ` | flot net par couche | couche fine | critère | verdict |
|---|---|---|---|---|
| 0 (M1) | +0,049 ; +1,091 ; +0,544 ; −1,844 | 1,844 | — | référence |
| 0,5 | +0,047 ; +0,714 ; −0,015 ; −0,637 | 0,637 | ∈ [1,65 ; 2,05] | **échec** |
| 1,5 | −0,703 ; +0,228 ; +0,262 ; +0,055 | 0,055 | ∈ [0,7 ; 1,2] | **échec** |

* **(M3) échoue aux deux points.** La couche fine passe de 1,844 à 0,637 puis 0,055 : pas de palier pour |μ| < 1, une décroissance continue. C'est l'« échec utile » prévu
  dans le critère : l'image « deux membres du doublet, masses m_tot ± μ, chacun porte une branche s'il est du côté topologique » est **contredite à N = 3**.
* Constat : quand μ croît, le flot **quitte la couche fine** et se reporte vers les couches grossières (μ = 1,5 : −0,703 sur la couche 0). La masse de goût déplace le survivant
  en profondeur plutôt qu'elle ne sépare un doublet. (Lecture, non testée.)
* Conséquence pour A3 V9 : la partie algébrique de T1 (A1)–(A5) reste acquise ([T]) ; la lecture « saveur = isospin » ne reçoit **pas** d'appui numérique de (M3) et reste [C] ;
  (M2) (N = 4, μ = 0, en cours sur le Mac) peut confirmer le facteur 2 à μ = 0, mais pas réparer (M3).

## (M2) — exécuté par JJ le 16/11 (≈ 20 min, solveur lobpcg) ; lecture Claude, 8 h
`T1_M2_N4_mu0.log` : N = 4, m_tot = −1, μ = 0, γ_F = +1, flot net par couche [+0,056 ; +1,195 ; +0,588 ; +0,043 ; **−1,811**].
Critère (écrit avant) : |flot(μ = 0)| / 0,904 ∈ [1,8 ; 2,2] ⇒ 1,811 / 0,904 = **2,00 : satisfait**. Le facteur 2 (deux branches par secteur sans masse de goût) ne dépend
pas de la taille entre N = 3 et N = 4. Il ne répare pas (M3).

## Correction de ma lecture de (M3) — Claude, 16/11/2026, 8 h 30 (Guadeloupe)
J'avais appliqué le critère de Fable sans le vérifier. Or **le critère était mal posé**, et (M3) n'est pas négatif : il est **non concluant**.

* **La règle utilisée.** C'est celle de la Prop `A3:prop:branches`, tirée de la chaîne réduite et démontrée [T] : une masse de Wilson (r = 1) porte un mode de bord si et seulement si elle est dans (−2, 0), avec la décroissance λ = |1 + m|.
* **Ce que prédit l'image « doublet »** à `m_tot = −1` : les deux membres ont les masses −1 ± μ.
  - Pour |μ| < 1, il y a **deux branches**, avec λ = |μ|.
  - Pour |μ| > 1, il n'y a **aucune branche**. Le critère de Fable attendait 0,7–1,2 à μ = 1,5 : c'était une confusion avec le point uniforme `m_tot = 0`, où μ = 1 donne une branche.
* **μ = 1,5** : on mesure 0,055 et l'image prédit 0. **Accord.**
* **μ = 0,5** : on mesure 0,637 et l'image prédit 2 branches de λ = 0,5. On sait déjà que l'arbre à ces profondeurs ne porte une branche que pour λ ≲ 0,3 :
  - la branche jumelle passe de −0,763 à ε = 0 (λ = 0,27) à −0,084 à ε = 0,2 (λ = 0,61) ;
  - à `m_tot = 0`, μ = 0,8 donne 0,583 à N = 3 puis 0,809 à N = 4.

  Une valeur basse est donc attendue par taille finie. **Non concluant.**

L'article est corrigé dans ce sens (remarque `A3:rem:isospin`, ligne 13c, résumé, conclusion).

## (M3′) — critères écrits AVANT tout calcul (16/11/2026, 8 h 30)
Même script, même mode : `bbhG`, `m_tot = −1`, `γ_F = +1`, 64 pas. Grandeur lue : |flot net de la couche fine|.

* **(M3′a)** `N = 3`, μ = 0,25. Masses −1,25 et −0,75, λ = 0,25 : on est dans la plage où l'arbre porte les branches.
  - Image doublet : |flot| ≥ 1,5 → **réussite**.
  - |flot| ≤ 1,0 → **échec** : l'image est contredite.
  - Entre les deux → non concluant.
* **(M3′b)** `N = 4`, μ = 0,5. On compare à 0,637 (N = 3).
  - |flot| ≥ 0,85 (croissance d'au moins un tiers, comme pour μ = 0,8 à `m_tot = 0`) → **réussite** : la valeur basse de N = 3 était de la taille finie.
  - |flot| ≤ 0,70 → **échec**.
  - Entre les deux → non concluant.
* Les deux réussites ensemble feraient passer la lecture « saveur = isospin » de [C] à [C] appuyée par trois calculs, sans en faire un théorème. Un échec de (M3′a) la contredirait vraiment.

Commandes :
```
cd ~/Desktop/2026FHT/FHT_26/03_code/CalculCodes
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=V331_relecture ../../../venv/bin/python V335_FHT_C3prime_flot.py 3 -1.0 0.25 1 64 1 bbhG | tee "../../../FHT_29 Final/00_RETOURS_FABLE/tests/T1_M3p_N3_mu025.log"
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=V331_relecture ../../../venv/bin/python V335_FHT_C3prime_flot.py 4 -1.0 0.5 1 64 1 bbhG | tee "../../../FHT_29 Final/00_RETOURS_FABLE/tests/T1_M3p_N4_mu05.log"
```
(environ 3 min puis environ 20 min, d'après (M3) et (M2))

## Résultats (M2) et (M3) — journaux du 16/11 (lus à 17 h 37 UTC)

| journal | flot de la couche fine | critère écrit avant | verdict |
|---|---|---|---|
| `T1_M2_N4_mu0.log` (`N = 4`, `m_tot = −1`, `μ = 0`) | −1,811 ; résidu 5,6·10⁻⁸ ; saturées 0 | `|flot| / 0,904 ∈ [1,8 ; 2,2]` | 2,00 : **réussi** |
| `T1_M3_N3_mu05.log` (`m_tot = −1`, `μ = 0,5`) | −0,637 | `|flot| ∈ [1,65 ; 2,05]` | **échec** |
| `T1_M3_N3_mu15.log` (`m_tot = −1`, `μ = 1,5`) | +0,055 (couche 0 : −0,703) | `|flot| ∈ [0,7 ; 1,2]` | **échec** |
| `T1_M3p_N3_mu025.log` (`μ = 0,25`, lancé par JJ) | en-tête seul, aucune ligne de résultat cinq heures après le lancement | — | à vérifier sur le Mac (une exécution `N = 3` dure 2 à 3 minutes) |

Lecture.
1. **Le facteur 2 à `μ = 0` tient à `N = 4`** (rapport 2,00) comme à `N = 3` (1,96) : sans masse de saveur, chaque secteur porte deux branches de Weyl. C'est robuste.
2. **L'image simple « deux membres de masses de paroi `m_tot ± μ`, chacun avec une branche si sa masse est dans la fenêtre topologique » échoue sur le critère que j'avais écrit**, qui supposait que `−0,5` et `−1,5` sont tous deux dans la fenêtre. À `μ = 0,5` on obtient 0,64 et non environ 1,8.
3. Lecture après coup, déclarée comme telle : les six points disponibles en mode `bbhG` à `N = 3` seraient compatibles avec cette image si la fenêtre topologique est plus étroite, à peu près `0,75 < −m < 2` (à `m_tot = 0` : `μ = 0,5` donne 0, `μ = 1` donne 0,94 ; à `m_tot = −1` : `μ = 1,5` donne 0,05). Ce n'est pas une réussite : c'est une hypothèse, que le test (M4) ci-dessous rend réfutable.

Conséquence pour les articles : « saveur = isospin » reste hors de A3 V9. Le compte de 48 Weyl par la lecture isospin n'est pas établi ; seul le facteur 2 à masse de saveur nulle l'est.

## (M4) — test d'additivité, critère écrit le 16/11 à 17 h 40 UTC, avant toute exécution

Hypothèse : le flot de la couche fine est la moyenne des flots de deux membres de masses de paroi `m_tot + μ` et `m_tot − μ`, chacun pris sans masse de saveur :
`flot(m, μ) = ½ [F₀(m + μ) + F₀(m − μ)]`, où `F₀(m)` est le flot à `μ = 0` (deux branches).

Mesure de `F₀` (mode `bbhG`, `N = 3`, `γ_F = +1`, 64 pas) aux sept masses `m ∈ {−2,5 ; −1,5 ; −0,5 ; +0,5 ; +1,0 ; +1,5 ; +2,5}` ; `F₀(−1) = −1,844` est connu.

Critère : avec ces `F₀`, la formule doit reproduire à 0,15 près les six points existants, tous en mode `bbhG`, `N = 3`, 64 pas :

| `m_tot` | `μ` | flot mesuré | source |
|---|---|---|---|
| −1 | 0,5 | −0,637 | `T1_M3_N3_mu05.log` |
| −1 | 1,5 | +0,055 | `T1_M3_N3_mu15.log` |
| 0 | 0,5 | +0,001 | `V335_logs.txt` |
| 0 | 1,0 | −0,943 | `V335_logs.txt` |
| 0 | 1,5 | −0,653 | `V335_logs.txt` |
| 0 | 2,5 | +0,371 | `V335_logs.txt` |

Réussite : les six écarts sont inférieurs à 0,15. Échec : deux écarts ou plus au-dessus de 0,3. Entre les deux : non concluant. Un échec réfute la lecture « deux membres d'un doublet de masses `m ± μ` », et donc l'interprétation isospin du doublement sous sa forme simple.

### Diagnostic du calcul `μ = 0,25` et changement de solveur (16/11, avant les exécutions de (M4))
Le processus `V335_FHT_C3prime_flot.py 3 -1.0 0.25 1 64 1 bbhG` tournait depuis cinq heures (≈ 500 % CPU). Cause probable : pour une dimension supérieure à 1 500, `V335` utilise `eigsh` en inversion décalée autour de `σ = 10⁻⁴` ; si une valeur propre tombe très près de `σ`, la factorisation est presque singulière et ARPACK peut ne pas converger. `V338_FHT_flot.py` (même Hamiltonien) utilise la diagonalisation dense jusqu'à la dimension 6 000, donc pour `N = 3` (dimension 4 680).
Décision : les exécutions de (M4) se font avec `V338`. Contrôle ajouté, critère écrit avant : `m_tot = −1`, `μ = 0` avec `V338` doit redonner −1,844 à 0,02 près ; sinon les points `V335` et `V338` ne sont pas comparables et (M4) est suspendu.
Mise à jour (16/11, 18 h UTC, avant tout résultat de (M4)) : la diagonalisation dense de `V338` à `N = 3` ne tourne que sur un cœur sur le Mac (≈ 100 % CPU, plus de 10 min pour le seul contrôle). On force le solveur lobpcg de `V338` (variable `SOLV=lobpcg`), celui de l'exécution `N = 4` réussie ; le contrôle (−1,844 à 0,02 près) et le critère de (M4) sont inchangés, et le résidu maximal imprimé doit rester inférieur à 10⁻⁶.

### (M4) — lecture (16/11, ~20 h 10 UTC), par le critère de 17 h 40, sans changement
Contrôle : m = −1, μ = 0 avec V338/lobpcg → −1.844, résidu 1.9e-8. Passé, points comparables.
F₀ de la couche fine (N = 3, bbhG, γ_F = +1) :

| m | −2.5 | −1.5 | −1 | −0.5 | +0.5 | +1.0 | +1.5 | +2.5 |
|---|---|---|---|---|---|---|---|---|
| F₀ | 0.000 | −0.855 | −1.844 | −0.015 | 0.000 | 0.000 | 0.000 | 0.000 |

Tous les résidus sont inférieurs ou égaux à 2.1e-8.

| m_tot | μ | prédit ½[F₀(m+μ)+F₀(m−μ)] | mesuré | écart |
|---|---|---|---|---|
| −1 | 0.5 | −0.435 | −0.637 | 0.20 |
| −1 | 1.5 | 0.000 | +0.055 | 0.06 |
| 0 | 0.5 | −0.008 | +0.001 | 0.01 |
| 0 | 1.0 | −0.922 | −0.943 | 0.02 |
| 0 | 1.5 | −0.428 | −0.653 | 0.23 |
| 0 | 2.5 | 0.000 | +0.371 | 0.37 |

**Verdict selon le critère : NON CONCLUANT.** Trois écarts dépassent 0.15, donc ce n'est pas une réussite ; un seul dépasse 0.3, donc ce n'est pas un échec.
Constats, sans interprétation :
- l'additivité reproduit trois points à 0.06 près, dont (0, 1.0) ;
- elle sous-estime (−1, 0.5) et (0, 1.5) d'environ 0.2 ;
- elle ne peut pas produire le flot positif mesuré en (0, 2.5), puisque tous les F₀ sont ≤ 0 ;
- F₀(−0.5) ≈ 0 alors que F₀(−1.5) = −0.855, à la même distance |1 + m| = 0.5 du centre de la fenêtre de Wilson : l'image de la chaîne réduite (fenêtre (−2, 0) symétrique autour de −1) ne s'applique pas telle quelle à l'arbre à N = 3.
La lecture isospin reste une conjecture non confirmée. Le test suivant prévu est N = 4, μ = 0.5 (M3′), lancé par JJ vers 19 h UTC.

### (M3′b) — lecture (16/11, ~20 h 30 UTC), par le critère de 8 h 30, sans changement
`T1_M3p_N4_mu05.log` : N = 4, m_tot = −1, μ = 0.5, bbhG, γ_F = +1, 64 pas ; dim 37 448 ; solveur lobpcg ; résidu 3.3e-8 ; durée 2 681 s.
Flot net par couche : [+0.036 +0.489 +0.091 −0.128 **−0.637**].
**|flot| = 0.637 ≤ 0.70 : ÉCHEC selon le critère.**
Le flot de la couche fine est identique à 10⁻³ près à celui de N = 3 (−0.637, `T1_M3_N3_mu05.log`). Ce sont deux exécutions distinctes : dimensions 4 680 et 37 448, profils par couche différents.
La valeur basse n'est donc pas un effet de taille finie : l'hypothèse post hoc de la ligne 176 (« λ = 0.5 hors de la plage où l'arbre porte une branche, effet de taille ») est contredite pour ce point.
Bilan de la lecture isospin :
- algèbre [T] ;
- deux branches à μ = 0 [E] ;
- (M3) échec ;
- (M4) non concluant ;
- (M3′b) échec.
Seul (M3′a) (N = 3, μ = 0.25) reste ouvert. La lecture est rétrogradée de « conjecture » à « conjecture contredite par deux tests numériques sur trois conclusifs ».

### (M3′a) — exécution (16/11, ~20 h 40 UTC, avant tout résultat)
Changement de solveur déclaré, le même que pour (M4) : V338 avec SOLV=lobpcg, car V335 (eigsh shift-invert) a figé 5 h sur ce point.
Critère de 8 h 30 inchangé : |flot fin| ≥ 1.5 réussite ; ≤ 1.0 échec ; entre les deux non concluant. Le résidu maximal doit rester < 10⁻⁶.
Journal : T1_M3p_N3_mu025.log (l'ancien, qui ne contient que l'en-tête, est écrasé).

### (M3′a) — lecture (16/11, ~22 h 30 UTC), par le critère de 8 h 30, sans changement
`T1_M3p_N3_mu025.log` : N = 3, m_tot = −1, μ = 0.25, bbhG, γ_F = +1, 64 pas ; V338, SOLV=lobpcg ; durée 231 s.
Le journal imprime un résidu max de 1.4×10⁻⁸ (< 10⁻⁶).
lobpcg signale que deux vecteurs n'atteignent pas sa tolérance interne de 10⁻⁷ (6×10⁻⁶ et 3×10⁻⁶). Déclaré ; c'est sous le seuil de 10⁻⁶ du critère pour les modes de la fenêtre.
Flot net par couche : [+0.049 +1.033 +0.424 **−1.406**].
**|flot| = 1.406, entre 1.0 et 1.5 : NON CONCLUANT selon le critère.**
Constat, sans interprétation : à N = 3 et m_tot = −1, le flot de la couche fine décroît de façon monotone avec μ :
1.844 (μ = 0), 1.406 (0.25), 0.637 (0.5), 0.055 (1.5).
L'image « doublet de masses −1 ± μ », pour laquelle deux branches persistent tant que |μ| < 1, prévoyait un plateau, et non cette décroissance régulière.

## Bilan final de la lecture « saveur = isospin » (16/11, 22 h 30 UTC)
- Algèbre : [T].
- Deux branches à μ = 0 : [E], N = 3 et 4.
- Tests en masse de saveur :
  - (M3) échec ;
  - (M4) non concluant ;
  - (M3′b) échec, la valeur basse n'étant pas un effet de taille ;
  - (M3′a) non concluant.
- Aucun test ne réussit, deux échouent. La lecture reste [C], contredite sous sa forme simple (doublet de masses m ± μ). La masse de saveur reste un choix posé dans A3.
