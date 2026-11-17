# Relecture indépendante de V152, V154 et V149 (§3.2 Ec4, §6) — piste B, structure réelle, trilemme

Date : 2026-11-11. Relecteur indépendant (je n'ai pas écrit les calculs relus). Notes non modifiées.
Scripts de vérification : `/home/claude/fht/relecture/` (`rel_*.py`, journaux `rel_*.log`).
Étiquettes : **[E]** recalculé par moi avec mon propre code ; **[T]** vérifié par lecture du code ou par preuve ; **[C]** lecture ou interprétation.

## Résumé des verdicts

| point | objet | verdict |
|---|---|---|
| 1 | B1 : hamiltonien, secteurs, χ, comptages N = 3, 4 (et N = 5 partiel) | **confirmé avec réserve** (écart de formulation V150/V151 ; à N = 4, 12 partenaires d'hélicité − sont sous la fenêtre) |
| 2 | Contrôle plat 66 + 66 | **confirmé** ; contrôle équitable |
| 3 | Partenaires « à l'infini comme 2^N » | **confirmé** (énergie, jusqu'à N = 9 par réduction exacte) ; **localisation « vers la profondeur 3 » non confirmée** |
| 4 | R4 : J² = +1, JH = −HJ, ordre un, masses εμ | **confirmé** ; la loi des masses est une **identité exacte triviale** (pas un test du bord) ; J_∂² = −1 non calculé par les auteurs (vrai) |
| 5 | B2 : « bord relativiste seulement à m = −r ; z < 1 ailleurs » | **erreur d'interprétation (effet de taille finie)** pour m = −0,5 et −1,5 : v → 1,50, z → 1, E(1,1,0)/E(1,0,0) → √2. Anomalie réelle seulement si \|1 + m/r\| > 1/√2 |
| 6a | Ec4 = 12 sur 6 états neutres | **confirmé** (preuve à la main et calcul avec base complète) ; l'explication du « 42 » donnée en V149 est inexacte |
| 6b | V149 §6 « trilemme » (théorème) | **erreur partielle** : s(k*) = 0 est juste ; le pas (i) est faux, contre-exemple explicite ; « 8 nœuds » n'est pas un théorème ; seule la conclusion « spectre vectoriel pour D local » tient (Nielsen–Ninomiya) |
| 7 | Divers | multiplicité « 20 » à \|n\|² = 13 (artefact d'affichage) ; levées de dégénérescence réelles (pas du bruit) ; paire k = 0 non découplée à \|1+m/r\| = 1/2 ; B0 et pochoir « roi » non vérifiés |

---

## 1. B1 (V152 §0–1)

### 1.1 Hamiltonien [T]
Lecture de `V151_FHT_B1_paroi.py` :
* H₊ = Σ_j [(T_j − T_j†)/(2i) ⊗ Γ_j − (r/2)(T_j + T_j†) ⊗ Γ₅] + (S − S†)/(2i) ⊗ Γ₄ − (r/2)(S + S†) ⊗ Γ₅ + (m + 4r) Γ₅, avec (T_jψ)(x) = ψ(x + e_j) (symbole sin k_j), Wilson r(1 − cos k) par direction, terme sur place m + 4r (3 latérales + verticale). Signes corrects ; M(k = 0) = m.
* S : (Sψ)(c) = ψ(parent c)/√8 ; S†S = 1 (isométrie, moyenne sur les 8 frères). Conforme à V151 §2.1.
* À la racine (côté 1) et à la profondeur 1 (côté 2), T = T† : la partie cinétique latérale s'annule et le Wilson donne −r Γ₅ par direction. C'est cohérent (et identique dans mon code).
* **Réserve de formulation [T]** : V150 §2 écrit H_X = Σ(Γ ⊗ sauts + Wilson) + m γ_F ⊗ Γ₅, avec le **Wilson non multiplié par γ_F**. Le code (et V151 §2.1, qui le dit explicitement) multiplie **masse et Wilson** par γ_F. Ce n'est pas anodin : avec la formule littérale de V150, le secteur γ_F = −1 aurait une masse effective −m = +1 avec Wilson +r, donc serait dans la phase triviale et n'aurait **aucun** mode de bord. Le calcul suit la bonne version (V151), mais V150 §2 est à corriger.

### 1.2 Réduction par secteurs [E + T]
* Les demi-translations commutent avec la relation parent–enfant (parent(c + 2^{n−1}) = parent(c) + 2^{n−2}) [T] ; commutation exacte recalculée avec mes propres opérateurs : ‖[τ_d, H]‖ = 0 (N = 3, 4) [E].
* La base de secteurs des auteurs est orthonormée et complète (somme des dimensions = 4 N_c dans leurs journaux) [T].
* Surtout : j'ai diagonalisé **l'espace complet, sans secteurs** (N = 3 dense ; N = 4 shift-invert, 220 valeurs propres jusqu'à |E| = 0,97, résidu 2·10⁻¹²), avec mon propre constructeur et une **autre représentation des Γ** (Γ_j = σ_j ⊗ τ_z, Γ₄ = 1 ⊗ τ_x, Γ₅ = 1 ⊗ τ_y ; même classe de Cl(5) : Γ₁…Γ₅ = −1). Les comptages et énergies sont identiques à ceux des auteurs ⇒ la réduction n'a rien perdu [E].

### 1.3 Hélicité [T + E]
χ = −iΓ₄Γ₅ commute avec Γ₁, Γ₂, Γ₃ et anticommute avec Γ₄, Γ₅. Un mode de bout de la chaîne verticale (Γ₄, Γ₅) est vecteur propre de χ ; sur χ = ±1, Γ_j agit comme ±σ_j, donc H_bord = ±σ·k : χ est bien l'hélicité des états d'énergie positive [T]. Vérifié : Γ₅H₊Γ₅ = −H₋ exactement [E] ; mes valeurs χ = ±1,00 sur les modes nuls après diagonalisation de χ dans les sous-espaces dégénérés [E].

### 1.4 Comptages [E]
| N | fenêtre | bord, χ + | χ − hors racine | racine, χ − | source |
|---|---|---|---|---|---|
| 3 | 1,100 | 14 | 12 (0,81518, n̄ = 2,49) | 2 | mon code, espace complet ; identique au journal |
| 4 | 0,778 | 66 | **12** (0,74232, n̄ = 2,85) | 2 | mon code, espace complet |
| 5 (secteurs 100 et 000 seulement) | 0,421 | 4 et 14 | 0 (plus bas : 0,73917 et 0,74232) | 2 | mon code, mes secteurs |

* Énergies du bord N = 4 : 0,38438 (×6), 0,55024–0,55025 (×12), 0,68520 (×8), 0,73950 (×6) par signe, plus 4 états exactement nuls (2 au bord χ = +1, 2 à la racine χ = −1).
* **Réserve sur l'énoncé « 66 sous la fenêtre, 0 miroir hors racine, 2 à la racine »** : ce n'est vrai qu'à **N = 5**. À N = 4, 12 états d'hélicité − (les partenaires, 0,74232) sont **sous** la fenêtre (0,778), juste au-dessus de la couche |n|² = 4 (0,73950). Le tableau de V152 §1.2 le dit correctement ; la phrase du bilan V152 §0.1 (« sous la fenêtre, il ne reste que 2 états à la racine ») est ambiguë et ne vaut qu'à N = 5. Ces 12 états sont d'ailleurs classés « bord » par le critère n̄ ≥ N − 1,5.
* **La fenêtre ne cache pas de miroir** : à N = 4, dans l'espace complet jusqu'à |E| = 0,97, aucun état χ − hors racine sous 0,742 [E] ; au-dessus, une bande dense de χ − (0,815, 0,838, 0,904, 0,963…) à n̄ ≈ 3,5–3,8. À N = 5, le plus bas χ − hors racine est 0,739 dans les secteurs que j'ai recalculés, et dans la réduction exacte (§3). La fenêtre (1,1 × couche |n|² = 4) est donc adaptée ; mais l'« absence de miroir » est une affirmation **sous une énergie fixe ≈ 0,74 en unités de l'arbre**, pas à toute énergie.
* N = 5 complet (8 secteurs) **non recalculé** (coût) ; secteurs 010 et 001 équivalents à 100 par symétrie cubique ; 110, 101, 011, 111 non recalculés.

### 1.5 Isotropie (V152 §0.2, §1.2) [E partiel]
* Écarts arbre/plat confirmés pour |n|² = 1 (N = 4 : +0,44 % ; N = 5 : +0,04 %).
* **Les dégénérescences ne sont pas « exactes à 10⁻⁵ près » au sens d'un bruit numérique** : il y a de vraies levées, petites. À N = 4 (résidu 2·10⁻¹²), la couche |n|² = 2 se scinde en 0,55024 ×6 et 0,55025 ×6, et la couche suivante en 0,85909 ×12 et 0,86022 ×12 (0,13 %). C'est attendu : la symétrie exacte est (ℤ/2)³ × groupe cubique, dont les irréductibles sont de dimension ≤ 4 ; rien n'impose la dégénérescence 12 ou 24. À reformuler : « levées réelles, de 10⁻⁵ relatif pour les premières couches à N = 5 ».
* **Erreur mineure** : V152 §1.2 donne les multiplicités « 6, 12, 8, 6, 24, 24, 12, 30, 24, 24, 8, **20** » et dit que ce sont « exactement les nombres de vecteurs entiers ». Or r₃(13) = 24. Le 20 vient de la troncature de l'affichage (NPRINT = 400 lignes ; la dernière ligne imprimée est au milieu de la couche 13) [T]. Ce n'est pas un résultat ; il faut retirer le 20 ou recompter.

## 2. Contrôle plat (V152 §1.1) [E]
Recalcul par impulsion (mon code) : N = 4 et N = 5 : 66 états χ + sur la couche N et 66 états χ − sur la couche 0 (à N = 4, deux états nuls dégénérés sortent mélangés de la diagonalisation, 65 + 65 + 2) ; E(|n|² = 1) = 0,38268 et 0,19509. **Confirmé.**
Équité [T] : même pochoir, même Wilson, même m ; la seule différence est la géométrie (couches de 2^N cellules, S = identité). Pour les modes lisses, S†S = 1 dans les deux cas : la chaîne verticale vue par un mode lisse est la même. Le contrôle est équitable.

## 3. Partenaires (V152 §1.3, V154 §2.2)
**Outil [T + E]** : le sous-espace des fonctions indépendantes de (y, z) est **invariant** par H (les sauts en y, z et S, S† préservent les constantes). Sur lui, le modèle 3D se réduit **exactement** à un arbre binaire (poids 1/√2), avec Wilson–Dirac à 1 direction latérale. La couche n = (1,0,0) est dans le secteur impair. Vérifié : les énergies réduites figurent dans le spectre 3D complet à 4·10⁻¹⁶ près (N = 3) ; elles redonnent **tous** les chiffres des auteurs (0,73950 ; 0,38438 ; 0,19516 ; partenaires 0,81518 ; 0,74232 ; 0,73917). On peut ainsi aller à N = 9 (dense) ou N = 11 (shift-invert).

| N | m = −1 : partenaire (n̄) | rapport au mode physique | m = −0,5 : partenaire | rapport | m = −1,5 : partenaire | rapport |
|---|---|---|---|---|---|---|
| 5 | 0,73917 (2,91) | 3,79 | 0,59411 | 2,21 | 0,74796 | 3,12 |
| 6 | 0,73632 (3,55) | 7,51 | 0,52159 | 3,74 | 0,68922 | 5,29 |
| 7 | 0,73332 (4,10) | 14,9 | 0,48366 | 6,77 | 0,62210 | 9,06 |
| 8 | 0,73169 (4,55) | 29,8 | 0,46397 | 12,8 | 0,57788 | 16,3 |
| 9 | 0,73060 (5,02) | 59,5 | 0,45392 | 24,9 | 0,55041 | 30,5 |

* **Énergie [E]** : en unités de l'arbre, l'énergie du partenaire reste d'ordre 1 (≈ 0,73 à m = −1 ; ≈ 0,45–0,55 à m = −0,5 et −1,5, encore en lente descente) ; le rapport au mode physique double à chaque niveau. « Les partenaires partent à l'infini comme 2^N en unités physiques » est **justifié** pour cette couche jusqu'à N = 9, aux trois masses. Pour |n|² = 1 il n'y a, jusqu'à N = 9 et sous 0,7, **aucun** état χ − dans ce sous-espace à m = −1.
* **Localisation [E]** : la lecture « les partenaires se placent vers la profondeur 3 » (V152 §1.3) **ne tient pas au-delà de N = 5**. La profondeur moyenne croît avec N (2,9 → 5,0 de N = 5 à 9 à m = −1) et le profil s'étale : w(n) = [0 0,01 0,06 0,16 0,16 0,20 0,20 0,13 0,06 0,02] à N = 9. Le mécanisme proposé (« la chaîne devient triviale à la profondeur où k atteint l'échelle du réseau », profondeur fixe) n'est donc pas celui que montrent les données [C].
* Convention [C] : « unités physiques » suppose le saut latéral uniforme à toutes les profondeurs (métrique de type hyperbolique), et la maille du bord 2^{−N}. C'est la convention du modèle ; c'est elle qui fait diverger le rapport. À mentionner.
* Limite : la réduction ne couvre que les impulsions (n_x, 0, 0) et (n_x, n_y, 0) ; les autres directions ne sont vues qu'aux N ≤ 5 des auteurs.

## 4. R4 (V154 §1)
### 4.1 Modèle F [T + E]
* Le F de `V148_FHT_e_raccord.py` (repris par R4) est **identique** à celui de `V146_FHT_produit_bord.py` : mêmes définitions, même graine 146, même ordre des tirages jusqu'à D_F [T].
* Il est conforme à V146 §0 : I = A ⊗ B ⊗ 3, γ_A = (+, −, −, +), γ_F(II) = −γ_F(I) ; π agit sur A dans I par (λ, q, λ̄) et sur B dans II par (λ, m) ; J_F échange I et II avec conjugaison ; Yukawa s_a ↔ d₁ et s_b ↔ d₂ ; Majorana symétrique sur (a, b) = (0, 0), qui est bien l'état neutre ν_R [T].
* J'ai **reconstruit F à partir de la description de V146**, sans copier le code, avec une autre graine : J_F² = 1, J_Fγ_F = −γ_F J_F, J_F D_F = D_F J_F, D_F impair et hermitien ; ordre zéro et ordre un **exacts** sur la base réelle complète de A_F (24 × 24 paires), avec ‖[D_F, a]‖ = 3,2 en témoin. Ce sont les signes de KO-dimension 6 [E].

### 4.2 Structure réelle totale [T + E]
* Preuve : avec U_X = Γ₂Γ₄ (représentation des auteurs) on a U_X Γ̄_a U_X⁻¹ = Γ_a pour a = 1…5 ; U_X est unique à une phase près (1 solution, recalculé dans ma représentation, où U_X = Γ₂Γ₅) ; J_X² = −1 ; (J_Fγ_F)² = −1 ⇒ J² = +1. Ensuite J K J⁻¹ = −K (T réel, i ↦ −i), J (W ⊗ γ_F) J⁻¹ = −W ⊗ γ_F, J (Γ₅ ⊗ D_F) J⁻¹ = −Γ₅ ⊗ D_F ⇒ JH = −HJ [T].
* Recalcul (N = 3, dimension 224 640, mon F, mes Γ) : ‖UŪ − 1‖ = 0, ‖U H̄ U† + H‖ = 0 (‖U H̄ U† − H‖ = 6,0 en témoin), ordre un pour A_F = 0 (8 × 6 éléments de base) [E].
* R4-P2 : J envoie les vecteurs propres de H₊ (40 plus petites |E|, N = 3) sur ceux de H₋ avec E ↦ −E, à 1,3·10⁻¹⁴, et χ ↦ −χ exactement [E]. La partie « J_∂² = −1 » de R4-P2 n'a **pas** été calculée par les auteurs (absente du journal) ; elle est vraie : avec l'identification H₋ = −Γ₅H₊Γ₅, J_∂ = Γ₅J_X vérifie J_∂² = U_X Ū_X = −1 [T, E].
* Remarque : la mesure d'hélicité « [−0,999 ; +0,921] » du journal mélange les modes de bord et les partenaires ; elle est compatible avec « hélicité retournée état par état », que j'ai vérifiée directement.

### 4.3 Masses de bord [T + E] : confirmé, mais trivial
Avec D_F impair, les termes anticommutent : {K ⊗ 1, Γ₅ ⊗ D_F} = 0 et {W ⊗ γ_F, Γ₅ ⊗ D_F} = 0 (W ∝ Γ₅, γ_F D_F = −D_F γ_F). D'où l'**identité opératorielle exacte**
H² = H₀² + ε² (1 ⊗ D_F²), avec H₀ = K ⊗ 1 + W ⊗ γ_F,
vérifiée à 3·10⁻¹⁵ (N = 3, F complet) [E]. Comme H₀(γ_F = −1) = −Γ₅ H₀(γ_F = +1) Γ₅, on a E² = E₀² + (εμ_i)² **pour tous les états** (bord, racine, partenaires, bulk), à tout ε, sans correction O(ε²). Les accords « à 4 ou 5 chiffres » de V154 sont donc ceux de l'affichage. Conséquence de lecture : « le bord reproduit D = D_F ⊗ 1 + γ_F ⊗ σ·p, masses comprises » est vrai, mais c'est une propriété algébrique du plongement Γ₅ ⊗ D_F, pas un test dynamique du bord. Elle donne aussi la masse εμ aux 2 états de racine. La prédiction R4-P3 (« corrections O(ε²) ») était plus prudente que nécessaire.

## 5. B2 (V154 §2–3) : vitesse et exposant z
* **Définitions [T]** : v_N = E₁/sin(2π/2^N) et z_eff = log₂(E₁(N−1)/E₁(N)), comparé au réseau plat. Les définitions sont sensées, mais à N = 3–5 elles mesurent surtout un **régime de crossover**.
* **Recalcul exact à grand N [E]** (réduction du §3, mêmes nombres que les auteurs à N = 3–5) :

| m | λ = \|1 + m/r\| | v_N (N = 5 ; 7 ; 9 ; 11 ; 14) | z_eff (4→5 ; 6→7 ; 8→9 ; 13→14) |
|---|---|---|---|
| −1 | 0 | 1,0003 ; 1,0000 ; 1,0000 ; 1,0000 ; 1,0000 | 0,978 ; 0,998 ; 1,000 ; 1,000 |
| −0,5 | 0,5 | 1,375 ; 1,456 ; 1,489 ; 1,497 ; 1,500 | 0,929 ; 0,964 ; 0,994 ; 1,000 |
| −1,5 | 0,5 | 1,228 ; 1,400 ; 1,471 ; 1,492 ; 1,499 | 0,828 ; 0,922 ; 0,986 ; 0,999 |
| −0,25 | 0,75 | 1,998 ; 2,968 ; 4,502 ; 6,617 ; 10,97 | 0,72 ; 0,70 ; 0,72 ; 0,77 |
| −1,75 | 0,75 | 1,575 ; 2,489 ; 3,846 ; 5,707 ; 9,581 | 0,66 ; 0,68 ; 0,71 ; 0,76 |

* Isotropie : E(1,1,0)/E(1,0,0) (arbre quaternaire exact) = 1,3990 ; 1,3976 ; 1,4024 à N = 5, 6, 7 pour m = −0,5, et 1,3586 ; 1,3712 ; 1,3853 pour m = −1,5 : retour vers √2 = 1,4142 [E].
* **Conclusion** : pour m = −0,5 et −1,5, le bord **devient linéaire et isotrope**, avec une vitesse renormalisée v∞ = 1,50 (z → 1). L'écart observé par les auteurs est un **effet de taille finie**. L'estimation au premier ordre v∞ = (1 − λ²)/(1 − 2λ²) (poids λ^{2j} des niveaux, impulsion de réseau 2^j p) donne 1,50 pour λ = 0,5 et 1,110 pour λ = 0,3 (mesuré : 1,110). Une vitesse v ≠ 1 commune à toutes les espèces (H₋ = −Γ₅H₊Γ₅, et le Yukawa ne change pas v) n'empêche pas un bord relativiste. Il y a bien un **régime anormal**, mais pour λ > 1/√2, soit m ∈ (−0,29, 0) ∪ (−2, −1,71) : la somme Σ(2λ²)^j diverge, v_N croît sans borne et z_eff ≈ 0,7–0,77 jusqu'à N = 14 [E pour les valeurs, C pour le critère].
* Donc « un bord relativiste n'est établi qu'au point spécial m = −r » (V154 §0.3, §3) est **trop pessimiste et faux pour |1+m/r| ≤ 0,5** ; B2a-P4 (« décroissant non établi ») est tranché : décroissant. La « question ouverte » de V154 §3 est résolue pour les masses testées.
* **Point nouveau, paire k = 0 [E]** : l'écart Δ de la paire bord–racine varie comme λ^N (Δ = 0,0235 ; 0,0117 ; 0,0059 ; 0,0029 ; 0,0015 ; 0,0007 pour N = 4…9 à m = −0,5), donc Δ/E₁ ∝ (2λ)^N. À λ = 0,5, Δ/E₁ reste fixe (≈ 0,04) : en unités physiques, le mode k = 0 du bord **ne découple pas de la racine**, il garde un écart fixe de 4 % de la première couche. Pour λ > 0,5 (m = −0,35 : 0,12 → 0,25 entre N = 4 et 9), le rapport croît. Le mode k = 0 n'est exactement nul en unités physiques que pour λ < 1/2. V154 présente cette paire comme anodine ; ce n'est pas un miroir propagatif (2 états), mais c'est un couplage bord–racine à l'échelle physique, à mentionner.
* Non recalculé : le pochoir « roi » (B2b) et ses écarts de +1,1 % à +12 %.

## 6. V149
### 6a. Ec4 (§3.2) [T + E] : confirmé
* **Preuve à la main** : comme bimodule (π, J_FπJ_F⁻¹), H_F se décompose en classes (irréductible gauche, irréductible droite). Sur I : (1_λ, 1_λ̄), (1_λ, 3̄), (2, 1_λ̄), (2, 3̄), (1_λ̄, 1_λ̄), (1_λ̄, 3̄). Sur II : (1_λ, 1_λ̄), (1_λ, 2), (1_λ, 1_λ), (3, 1_λ̄), (3, 2), (3, 1_λ). Chaque classe a la multiplicité 3 (générations). La seule classe commune à I et II est (1_λ, 1_λ̄) = {ν_R, ν̄_R}, de multiplicité 6. γ_F est constant sur chaque autre classe, donc un B impair du commutant vit sur ce bloc : B = [[0, C†], [C, 0]], et l'invariance par J_F impose C = Cᵀ. Résultat : 3 × 3 complexe symétrique, soit **12 dimensions réelles sur 6 états** (le terme de Majorana). Témoins : commutant de dimension réelle 2(10 × 9 + 36) = 252 ; pairs hermitiens J-invariants : 45 + 21 = **66** (valeur des auteurs).
* **Recalcul** (`rel_Ec4.py`, mon F, base réelle complète des 2 × 24 générateurs, autre méthode de réduction) : commutant 252, Ec4 **12**, témoin 66, sans J_F 18, support = les 6 états (a, b) = (0, 0) de I et II, C symétrique [E].
* **Ensemble générateur** : les 16 générateurs des auteurs (4 éléments aléatoires, parties hermitienne **et** antihermitienne, et leurs images par J) engendrent bien l'algèbre ; un seul élément aléatoire avec ses deux parties suffit déjà (12) [E].
* **Correction de l'explication** : le premier calcul (42) n'est pas faux parce que « 2 éléments n'engendrent pas l'algèbre ». Il est faux parce qu'il n'utilise que les **parties hermitiennes** a + a†, qui donnent 2 Re λ sur les places λ et λ̄ et ne les séparent pas. Avec 3 ou 6 éléments, parties hermitiennes seules : toujours 42 ; avec 2 éléments et leurs deux parties : 12 [E].

### 6b. §6, « le trilemme est un théorème » [T + E] : erreur partielle
* **Juste [T]** : pour un terme γ_F ⊗ (w + σ·s) invariant par translation, avec J = (J_Fγ_F) ⊗ (iσ_y ∘ conj) et k ↦ −k, KO-1 ⇔ s impaire et w paire (calcul direct). Sur (ℤ/N)³ avec N pair, s(k*) = 0 aux 8 TRIM, quelle que soit la localité.
* **Faux : le pas (i)**. « Un terme qui relève ces nœuds sans créer de surface nodale doit anticommuter avec γ_F ⊗ σ_j » est une condition **suffisante** (E² = |s|² + M²), pas nécessaire. **Contre-exemple [E]** (`rel_V149_trilemme.py`) : le terme 1_F ⊗ u(k)·σ avec u **paire**, par exemple u = (0, 0, r Σ_j(1 − cos k_j)). Il est hermitien, local, invariant par translation et invariant de jauge (identité sur F, ordre un exact) ; il respecte KO-1 exactement (‖JD(k)J⁻¹ + D(−k)‖ = 0 avec le F complet). Sur N = 8, il fait passer les nœuds de **8 à 2 par espèce** : k = 0 de chiralité +1 et k = (0, 0, −π/2) de chiralité −1. Il ne crée **pas** de surface nodale (2 points seulement avec |f| < 0,05 sur une grille 64³). Avec la dérivée SLAC (KO-1 exact, N = 8 et 16), le même terme donne **1 seul nœud par espèce**.
* **Conséquences** :
  - « Toute espèce chargée garde 8 nœuds » n'est **pas un théorème** ; ce qui reste est Nielsen–Ninomiya : pour un D local, chaque espèce a au moins un nœud de chiralité opposée, et le spectre est vectoriel. « Donc 4 ont l'hélicité inversée » n'est vrai que pour le symbole sin k.
  - Le « trilemme à N pair » (V149 §3.3) ne porte que sur l'ansatz où toute la partie spinorielle est multipliée par γ_F. Si l'on autorise un terme 1_F ⊗ σ·u pair, KO-1 et N pair ne sont plus un obstacle supplémentaire : il reste le dilemme usuel localité (commutateurs bornés) / absence de doubleur.
  - Les déductions qui s'appuient sur §6 doivent être revues dans ce sens : V149 §6 (correction de V144 §0) et V150 §1 (« théorème des TRIM »).
* **Réserve sur le contre-exemple [C]** : il brise la symétrie cubique (direction z privilégiée). Je n'ai pas cherché de version invariante par le groupe cubique ; si l'on exige cette invariance, une obstruction pourrait revenir. C'est à énoncer comme hypothèse explicite si c'est l'intention.
* (ii) (Wilson scalaire γ_F ⊗ w ⇒ surface nodale) : correct par le théorème des valeurs intermédiaires [T]. (iii) repose sur Ec4, confirmé ; mais Ec4 ne restreint que les termes **impairs**, et le contre-exemple est pair.

## 7. Autres remarques
* « Partenaires vers la profondeur 3 », « hélicité − hors racine : 0 » : vrais à N = 5, ils doivent être datés en N. Le second ne tient qu'à N = 5 pour la fenêtre choisie (§1.4).
* V152 §3 : « Nielsen–Ninomiya n'est pas contredit » [T] : le système complet est un hamiltonien hermitien de dimension finie ; la chiralité manquante est portée par les partenaires (énergie ≈ 0,73, unités de l'arbre) et la racine.
* Le critère de classe « bord » (n̄ ≥ N − 1,5) range les partenaires de N = 4 en « bord » ; il vaudrait mieux classer par χ et par w(N).
* Contrôle trivial m = +1 : premier état à 1,18978, intérieur (journal) [T] ; non recalculé.
* **Non vérifiés** : B0 (V152 §2) ; pochoir « roi » (V154 §2.3) ; N = 5 complet (6 secteurs sur 8 non recalculés) ; V149 (c), (d), Ea, Eb, Ec1–Ec3, Ec5–Ec6 hors §3.2 et §6.

## 8. Fichiers
`rel_tree.py` (arbre et K, W, Γ propres) ; `rel_B1.py` (+ `rel_B1_N3_tree.log`, `rel_B1_N4_tree.log`) ; `rel_flat.py` (+ .log) ; `rel_B1_N5_sect.py` (+ .log) ; `rel_F.py` (F reconstruit) ; `rel_R4.py` (+ .log) ; `rel_Ec4.py` (+ .log) ; `rel_B2_toy1d.py`, `rel_B2_reduction.py`, `rel_B2_partner_profile.py`, `rel_B2_k0pair.py` (+ .log) ; `rel_V149_trilemme.py` (+ .log).
