# Vérification indépendante de A3 (point 1) — 16/11/2026

Code : `verif_independante_A3.py` (+ `verif_c_modes_zero.py`), numpy seul, écrit à partir des seules définitions de A3 V9. Le code V352/Fable n'a été ni importé ni lu. Les générateurs de Cl(6) sont pris dans une représentation de Jordan–Wigner (jauge de flux π), différente de celle de V352. La comparaison est légitime parce que Cl(6) n'a qu'une seule représentation irréductible de dimension 8, à équivalence près. Sorties : `verif_a.json`, `verif_b.json`, `verif_c.json`, `verif_c_modes_zero.json`.

## (a) Thm P7:thm:bbh — **confirmé**
Le réseau octupolaire est construit en espace réel : (2L)³ sites (L=4), sauts t_f intra-cellule et t_c inter-cellule, signes x:+1, y:(−1)^x, z:(−1)^{x+y}. Le spectre complet coïncide avec E² = Σ_d(t_f²+t_c²+2t_f t_c cos k_d), avec ±E chacun quatre fois.

| t_f, t_c | écart max | \|E(K)\| | √3\|t_f−t_c\| |
|---|---|---|---|
| (1, 1) | 2.3e-14 | 0 | 0 |
| (1.3, 0.7) | 2.7e-14 | 1.039230 | 1.039230 |
| (0.6, 1.1) | 2.3e-14 | 0.866025 | 0.866025 |

## (b) Prop A3:prop:heat, niveau plat seul — **confirmé**
Formule analytique : z = 8[e^{−2t} I₀(2t)]³ par cellule. La même valeur sort d'une somme de Bloch sur L=64.

| t | n_eff | d_s |
|---|---|---|
| 10 | 8.156 | **3.040** (A3 : 3.04) |
| 20 | 8.076 | 3.019 |
| 40 | **8.038** (A3 : 8.04) | 3.010 |

## (c) Prop A3:prop:heat, arbre N=2, m=8 — **confirmé, avec une précision à ajouter**
D_s est construit littéralement à partir de la formule de la sous-section *Setting* :
- T_d sur chaque couche ;
- V parent→8 enfants, de poids 1/√8 ;
- jauge alternée t_c → −t_c ;
- t_f = t_c = r = μ = 1 ;
- moyenne de Bloch exacte sur k₀ ∈ (2πℤ/8)³.

Les deux secteurs s = ±1 donnent les mêmes valeurs.

**Modes zéro.** Il y en a 4 par secteur, tous à k₀ = 0 :
- 2 vivent sur la racine (poids 0.97–0.999 sur la couche 0) : c'est la « contrepartie » ;
- 2 vivent sur la couche fine (poids 0.94–1.0) : c'est l'échantillon q = 0 du fermion de Weyl.

Le tableau de A3 ne retire que les 2 modes de la racine. Avec la même prescription :

| t | 10 | 15 | 20 | 30 | 40 | 60 |
|---|---|---|---|---|---|---|
| n_eff (ici) | **2.158** | 2.096 | **2.072** | **2.052** | **2.060** | 2.205 |
| n_eff (A3, N=2) | 2.158 | — | 2.072 | 2.052 | 2.060 | — |
| d_s (ici) | 3.22 | 3.10 | 3.07 | 3.02 | 2.90 | 2.33 |

La ligne N=2 de A3 est donc reproduite à 10⁻³ près.

**Précision à porter dans le texte.** À N=2, la couche fine a un côté de 32. Les 2 modes zéro de la couche fine, gardés à juste titre, pèsent 2/(8³·8²)·(4πt)^{3/2} dans n_eff, ce qui fait 0.09 à t=10 mais **0.69 à t=40**. Le plateau à N=2 au-delà de t ≈ 25 est donc tenu en partie par ce seul échantillon q=0 : on y lit le volume fini et non plus le continuum. Le décollage à t=60 (2.205) le confirme, et c'est pourquoi A3 n'en donne pas la valeur.

Ordre de grandeur du volume fini : (L/2π)² ≈ 26 à N=2 (L=32) et ≈ 58 à N=3 (L=48). La fenêtre fiable est donc t ≲ 25 à N=2 et t ≲ 50 à N=3. Le tableau de A3 reste correct, mais la phrase « window bounded by the finite volume » mérite ces deux chiffres.

## Non fait ici
- La ligne N=3 : 216 diagonalisations de taille 4680, environ 30 min, ce qui dépasse la VM. Le même script peut la faire sur le Mac en remplaçant `test_c()` par `test_c(N=3, m=6)` (à lancer seulement si on le juge utile).
