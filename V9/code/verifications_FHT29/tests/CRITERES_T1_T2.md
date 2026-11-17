# Critères des tests T1 et T2 (Fable, 16/11/2026)

Règle du programme : critères écrits avant le calcul ; ce qui est ajouté après est déclaré. Rien n'est modifié hors de `00_RETOURS_FABLE/` ; les scripts de `FHT_26/03_code/CalculCodes` sont lus, jamais réécrits.

## T1 — L'indice de saveur `E` de l'octet est-il l'indice d'isospin ?

Hypothèse testée (LOG_A3, A3-D11) : le « doublement de saveur » de A3 est le doublet d'isospin de A2 compté deux fois. Lecture testée : l'octet est tensorisé par `F' = (deux chiralités) ⊗ (lepton + trois couleurs) ⊗ (trois générations)`, de dimension 24, et non par `H_F` entier ; le doublet est fourni par `E`.

Critère écrit avant tout calcul (LOG_A3, partie VI, rédigé avant l'ouverture du dépôt) : le compte donne 48 fermions de Weyl sans masse de saveur ajoutée.

Critères ajoutés le 16/11 **après** lecture du dépôt (déclarés post hoc pour (L) ; écrits avant exécution pour (A1)–(A5) et (M1)–(M2)) :
* (A1) le commutant de `{Γ_5, α_d}` contient une algèbre de quaternions réelle pour la structure réelle naturelle de l'octet ;
* (A2) il existe un opérateur de chiralité `γ` qui commute aux termes cinétiques et à cette algèbre, et anticommute à `Γ_5`, `Γ_4` ;
* (A3) `⟨γ⟩ ⊗ H` est de dimension réelle 8 (`H ⊕ H`), et `P_L H ⊕ P_R {1, AB}` de dimension 6 (`H ⊕ C`) ;
* (A4) à `k = K`, `t_f = t_c`, le spectre de `H + (m + μF)Γ_5` est `{m + μ, |m − μ|}`, quatre fois chacun ;
* (A5) les commutateurs de `H(k)` avec `F`, `P_L H`, `P_R AB` sont d'ordre `q²` près du cône.
* (L) lecture d'un journal existant (`V335_logs.txt`, post hoc) : le flot dans la couche fine à `μ = 0` vaut environ deux fois celui à `μ = 1`.
* (M1) calcul à lancer sur le Mac, secteur `γ_F = −1`, `μ = 0`, `N = 3`, `m_tot = −1`, mode `bbhG` : réussite si le flot de la couche fine est dans `[+1,65 ; +2,05]` (opposé de la valeur connue `−1,844`, à 10 % près).
* (M2) idem à `N = 4`, les deux secteurs : réussite si `|flot(μ = 0)| / |flot(μ = 1)|` est dans `[1,8 ; 2,2]`. Échec : un rapport proche de 1 (une seule branche) réfute la lecture.

Ce qui réfuterait la lecture, indépendamment des nombres : que les deux copies n'aient pas la même chiralité (elles l'ont, Prop 9.2 de A3) ; que `E` ne soit pas la représentation `2` de `S4` (il l'est, `v6_dictionary.py`) ; que `H ⊕ C` ne se plonge pas dans `⟨γ⟩ ⊗ H`.

## T2 — Trois générations sur un bloc à deux niveaux

État : **le calcul n'est pas encore bien posé**, et je ne dépose pas de script. Dans la lecture bimodule (LOG_A3, A3-D8), l'espace de Hilbert fini est un espace d'opérateurs (`End(S)`), et `D_F` y agit par multiplication à gauche et à droite ; « comprimer le Dirac natif d'un bloc `4³` » demande de dire d'abord sur quel espace d'états. Deux définitions possibles, à choisir avant d'écrire le code :
1. états = modes de basse énergie du Dirac natif sur le bloc `4³` au point `K` du niveau fin, classés par les nombres quantiques du niveau parent (trois directions) ;
2. états = opérateurs `Cl(A, B) · span{α_d, Γ_4} · Λ^{1,2}(c')`, avec `D_F` = commutateur gradué avec le Dirac natif.
Critère commun, écrit maintenant : les blocs non nuls ont la forme `S ⊗ 1` ou `1 ⊗ S` (Thm 5.2 de A2), le bloc entre `Λ¹(c')` et `Λ²(c')` est symétrique, et le groupe du niveau parent agit par `3 ⊕ 3'`. Échec utile : le niveau parent n'apporte qu'une multiplicité 8 non structurée, ou des blocs interdits apparaissent.

## T2 — Décision du 16/11 (prise par Fable à la demande de JJ) : définition 1

Je retiens la définition 1 : les états sont les modes du Dirac natif sur un bloc à deux niveaux, classés par les nombres quantiques du niveau parent. Raisons : c'est un calcul sur l'opérateur de l'article et non une construction ; il peut échouer ; la définition 2 réussirait par construction pour les caractères.

Objet : le réseau BBH uniforme de `V327`/`V335` (trois chaînes SSH, flux π), sur un bloc périodique `4³` (huit octets, un parent). Les nombres quantiques du niveau parent sont les moments `k ∈ {0, π}³` du réseau des octets.

Critère, écrit avant l'exécution numérique mais après le calcul à la main du spectre (déclaré) :
* **réussite** : le secteur le plus léger contient trois moments du niveau parent (un triplet), apparié par l'échange `k ↦ k + (π, π, π)` à un autre triplet de même `|E|` ;
* **échec** : le secteur le plus léger ne contient qu'un moment ; ou les deux triplets n'ont pas la même énergie ; ou le triplet est à l'échelle du réseau.
Contrôle : les multiplicités par niveau, à `t_f = t_c` et pour deux dimérisations.
