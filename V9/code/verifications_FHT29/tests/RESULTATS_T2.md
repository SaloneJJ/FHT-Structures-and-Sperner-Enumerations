# T2 — Résultats (Fable, 16/11/2026)

Critères : `CRITERES_T1_T2.md` (définition 1). Script : `T2_bloc_deux_niveaux.py` (numpy seul ; sortie dans `T2_bloc_deux_niveaux.log`).

## Résultat : échec selon le critère écrit

Sur le bloc périodique `4³` du réseau BBH uniforme, les moments du niveau parent sont `k ∈ {0, π}³`. Avec `w` le nombre de directions où `k_d = 0` :

`m_w = √( w (t_f + t_c)² + (3 − w)(t_f − t_c)² )`, avec les multiplicités (1, 3, 3, 1) en moments, huit états par moment.

| `t_f`, `t_c` | `w = 0` | `w = 1` | `w = 2` | `w = 3` |
|---|---|---|---|---|
| 1, 1 | 0 | 2 | 2,828 | 3,464 |
| 0,7, 1 | 0,520 | 1,752 | 2,423 | 2,944 |
| 0,4, 1 | 1,039 | 1,637 | 2,069 | 2,425 |

Le calcul numérique et la formule coïncident (64 valeurs propres dans chaque cas).

* Le secteur le plus léger ne contient qu'un moment (`w = 0`, le point `K`).
* Le triplet `w = 1` est à l'échelle du réseau (`2t` au point uniforme).
* Les deux triplets n'ont pas la même énergie ; l'échange `k ↦ k + (π, π, π)` envoie `H(t_f, t_c)` sur `H(t_f, −t_c)` et ne commute pas à `H`.

Conclusion : **les générations ne sont pas les moments du niveau parent du réseau uniforme.** Le candidat de LOG_A3 (A3-D12) est écarté sous sa forme la plus simple.

## Ce que le calcul apporte quand même

L'espace des étiquettes du niveau parent est bien `Λ•(R³)`, de dimensions (1, 3, 3, 1) : la structure cherchée existe, mais elle porte une hiérarchie de masses `m_w` et non une multiplicité. Deux remarques pour la suite :
1. l'échange `w ↔ 3 − w` est `t_c ↦ −t_c`, c'est-à-dire la jauge alternée (`bbhG`) de A3 : la dualité de Hodge du niveau parent relie les deux conventions de l'article, pas deux secteurs d'un même opérateur ;
2. sur l'arbre (et non sur le réseau périodique), le couplage vertical et l'échelle des sauts par niveau peuvent changer ce tableau ; je n'ai pas pu le calculer (le code de l'arbre demande scipy).

## Ce qui reste du pont `n = d + 1`

Deux mécanismes testés, deux échecs : les grandes diagonales d'un octet (A3-D9), les moments du niveau parent (ce test). Le pont reste une coïncidence de cardinaux (`|G0| = |X4| = d + 1`), sans mécanisme. Pour les articles, la décision est donc la sortie 1 de LOG_A2 (partie V) : `N_gen = C(n−1, n/2 − 1)` sous (Deg1), `n = 4` lu sur le gabarit, sans postulat (M), et sans affirmer d'origine spatiale des générations.
