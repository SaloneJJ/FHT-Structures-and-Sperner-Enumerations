# Relecture indépendante de V179–V183 et V192–V201 (règle D-E13)

Relecteur : agent séparé. Il n'a écrit ni les notes ni les calculs relus. Date : 2026-11-12.
Matériel : notes dans `relec4/notes/`, scripts et journaux dans `relec4/code/`.
Méthode : tous les points clés ont été **recalculés avec des scripts écrits par le relecteur**, dans `relec4/rv/` (`rv_common.py`, `rv1`…`rv8`, avec leurs `.log`). Aucune fonction des scripts relus n'est importée. Les relevés sont obtenus par un système linéaire (noyau de H·DP − DP·H), et non par propagation sur un arbre couvrant. Le champ ℤ₃ est construit en jauge axiale, et non par élimination de Gauss. L'arbre à deux niveaux utilise d'autres poids génériques, une jauge aléatoire et des bases aléatoires dans les niveaux propres. Les scripts relus n'ont pas été modifiés. Temps de calcul total : environ 1 min 30 s (dont 60 s pour l'arbre C0–C4).

Légende des verdicts : **confirmé**, **corrigé** (le calcul tient mais un énoncé est faux ou trop fort), **faux**, **non vérifié**.

---

## 1. Tableau des verdicts

| note | point | verdict | commentaire |
|---|---|---|---|
| V179/V180 | Q1 : octet ℤ/2, [U_x, U_y] = −1 | confirmé | −1,0000, scalaire (rv1) |
| V180 | Q1b : ±√3, deux quadruplets irréductibles | confirmé | Σ\|χ\|²/24 = 1 pour les deux |
| V180 | Q2 : U(1), m = 0…5 → (−1)^m | confirmé | +1, −1, +1, −1, +1, −1 |
| V180 | Q3 : γ_i² = 1, anticommutation, covariance | confirmé, **précisé** | U γ_i U† = **+**γ_{\|R\|(i)} pour les 24 rotations : les γ_i se transforment comme la permutation des axes (A₁ ⊕ E), **pas** comme un vecteur (T₁) |
| V180 | Q4 : K26 (26, 48, 24, χ = 2) ; ℤ/2 → +1 ; m = 1, 2, 3 → −1, +1, −1 | confirmé | m = 5 donne aussi −1 ; niveaux à m = 1 : 2, 4, 6 (Σ\|χ\|² = 1, 1, 2) |
| V180 | Q5 : tétraèdre π/2, −1, deux doublets ±√3 | confirmé | 3π/2 (m = 3) donne aussi −1 |
| V179/V180 | « il faut un flux U(1) d'ordre 24 » (K26) | corrigé (mineur) | m = 3 (π/4 par carré, ordre 8) donne −1 : ℤ₈ suffit |
| V180 §2 | « nombre impair de faces frustrées » | corrigé (mineur) | le cube a 6 faces à π (nombre pair). La règle porte sur le flux total/2π = m |
| V182 §0 | parité de m pour un champ ℤ/2 | confirmé | relevé non uniforme (5 faces à +π, 1 à −π, « m = 2 ») : même H à une jauge près, commutateur −1 |
| V183 | C0 : +1 | confirmé | rv2, autres poids, jauge aléatoire |
| V183 | C1 : aucun relevé | confirmé | 23/24 rotations sans relevé |
| V183 | C2 : relevé, −1, quadruplets | confirmé | dégénérescences basses 4, 4, 4, 4, 4, 4 |
| V183 | C3 : aucun relevé | confirmé | bloc racine de dimension 16, comme dans V183 |
| V183 | C4 : aucun relevé (prédiction −1 fausse) | confirmé | mais ce résultat était **prévisible** par V182 §1.1 (voir M4) |
| V183 §1 | « un seul twist par branche, au niveau le plus fin » | corrigé (moyen) | lecture trop forte (M4) |
| V183 | D-P1 : dim A₁ = 4, dim T₁ = 8 | confirmé | A₁ non commutative (≅ M₂) avec twist, commutative sans twist |
| V183 | D-P2/D-P3 : Cl(5) covariant, Γ₆ scalaire pour toutes les solutions | confirmé | 141 solutions sur 200 départs, toutes avec Γ₆ = ±1 ; Γ₁Γ₂Γ₃ invariant dans tous les cas |
| V183 | D-P4 : aucune solution sans twist | confirmé | **démontré** : dans une algèbre commutative, deux involutions ne peuvent pas anticommuter. Meilleur résidu numérique 7,5 |
| V183 | argument papier « Γ₆ scalaire car les scalaires forment M₂ » | confirmé, avec une étape manquante | il faut d'abord Γ₆ ∈ A₁, ce qui vient de Γ₁Γ₂Γ₃ invariant (det R = 1) (m4) |
| V183 | Γ₅ = H/√3 possible, famille à une phase | confirmé | H/√3 ∈ A₁ ; démonstration en §2.2 |
| V183 | D-P5 négatif (carré 2/3) | confirmé | **prévisible** : Σ_d γ⊥_d = 0, donc les γ⊥ n'engendrent que E (dimension 2) |
| V183 | D-P6 : doublet de goût = E | confirmé | caractères E 4, C₃ 1, C₂ 4, C₄ 0, C₂′ 0 |
| V183 §3 | « le bord natif est chiral » | non vérifié, énoncé trop fort | conditionnel, et repose sur un flot diabatique (V181) (M5) |
| V181 | intégration de la relecture de V167 | non vérifié | aucun script V166/V167 fourni |
| V192/V193 | Q1 : A²\|_R = 3I + 2(J − I) sans twist, 3I avec twist | confirmé | valeurs propres 9, 1, 1, 1 et 3, 3, 3, 3 |
| V193 | Q2-a : H(K) = (t_f − t_c) H_cube | confirmé | et H(0) = (t_f + t_c) H_cube |
| V193 | Q2-b : v = √(t_f t_c) | confirmé, **précisé** | spectre exact E² = Σ_d (t_f² + t_c² + 2 t_f t_c cos k_d), dégénéré 4 + 4 partout. v² = t_c² + m r_W. La vitesse de Dirac (coefficient de α) est t_c (m5) |
| V193 | décomposition avec m = √3(t_f − t_c), r_W = t_c/√3 | confirmé | **identité exacte pour tout q** (écart 3·10⁻¹⁰, venant de la dérivée numérique de α) |
| V193 | covariance des α_d | confirmé | α ∈ T₁ (résidu 10⁻¹⁵) ; U_g H(K+q) U_g† = H(K+gq) pour les 24 rotations |
| V193 | « frères → masse **et Wilson** » | faux (moyen) | le terme de Wilson vient des cousins (M2) |
| V193 | fenêtre de Kaplan 1/3 < t_f/t_c < 1 | faux (moyen) | artefact ; la phase est toute la région t_f/t_c < 1 (M1) |
| V193/V194 | identification au modèle octupolaire BBH | confirmé | six matrices 8 × 8 qui anticommutent, flux π par plaquette, γ = t_f, λ = t_c. Convention : vaut pour un amas fait d'octets complets |
| V183 ↔ V193 | les Γ natifs sont-ils les mêmes ? | confirmé, cohérent | les α_d sont orthogonaux aux sauts γ_d, comme le disait V183, et se complètent en un Cl(5) de type V183 avec un Γ₄ unique : Γ₄ = ±Γ₅α₁α₂α₃ |
| V195 | E1-a : 8 états dans le gap pour ε < 1, aucun pour ε > 1 | confirmé | y compris ε = 0,1 et 0,2, sous la borne 1/3 de V193 |
| V195 | E1-b : E = ±√3(1 − ε²)ε^n | corrigé (moyen) | loi **asymptotique**, pas exacte (M3) |
| V195 | E1-c : 4 coins pairs, 4 impairs | **tranché : oui** | dans la base localisée, chaque mode de coin est à 100 % sur le sous-réseau de son site de coin |
| V195 | « les 8 modes de coin forment un octet tordu effectif » | confirmé (l'argument donné ne suffit pas) | H_eff est exactement un cube de Kawamoto–Smit : π sur les 6 faces, couplages seulement entre coins voisins |
| V195 | E2 : 7 = A₂ + T₁ + T₂ ; pas de moyenne invariante avec twist | confirmé | caractères vérifiés à la main |
| V197 | échec de R1 et son explication | confirmé | le 65ᵉ niveau est égal au 64ᵉ ; l'explication est correcte |
| V198 | C1 (rapports 65ᵉ/64ᵉ) | confirmé | 8,62 … 1,88, identiques |
| V198 | C2 : résidu < 5·10⁻⁴ | confirmé, mais **ce n'est pas un test** | l'ajustement est exact par construction (G2) |
| V198 | C3 : t_f′ = E_bloc/√3, t_c′ = w τ | confirmé | t_f′ = δ = √(ε² + 1/4) − 1/2 ; w = 1/(1 + δ²/ε²) est le poids 1D du site d'extrémité (m8) |
| V198 | C4 : τ*(ε) ≈ ε | confirmé, **expliqué** | au premier ordre en τ, τ* = u(1 + u²) avec u = δ/ε, donc τ* = ε + O(ε⁵). Valeurs exactes : 0,1999 ; 0,299 ; 0,491 (et non 0,485) pour ε = 0,2 ; 0,3 ; 0,5 |
| V198 | « profil alterné » | **faux** | le calcul à deux pas donne t(4) ≈ τ², pas t(4) fort (G1) |
| V198 | carte ε ↦ ε′ bien définie | confirmé | ε′ = √(ab)/(a − b), solution unique, signes absorbés par la jauge |
| V199 | R = ε/τ | non dérivé (mineur) | m9 |
| V200 | cristal alterné ℤ₃ (Bianchi mod 3) | confirmé | construit indépendamment en jauge axiale ; toutes les plaquettes vérifiées |
| V201 | F1 : commutateur −1 pour φ = 0, ±1 | confirmé | |
| V201 | F2 : octet de charge 5 : ±√6 deux fois, 0 quatre fois ; H² non ∝ 1 | confirmé | lemme général : sur un cube à sauts unitaires, H² = 3·1 ⇔ flux π sur les 6 faces |
| V201 | F3 : 4n modes nuls, indépendants de ε | confirmé | 4, 8, 12, 16 pour ε = 0,3 ; 0,5 ; 1,0 ; spectre invariant de jauge |
| V201 | décision d'arrêt (règle de V200 §3) | défendable, mais hors de la lettre de la règle ; conclusion trop large (M6) | |
| toutes | horodatage « critères avant calcul » | incohérences (M8) | |

---

## 2. Chiffres recalculés et démonstrations

### 2.1 V180 (rv1)
* Octet, flux 2πm/6 : commutateurs +1, −1, +1, −1, +1, −1. Niveaux :
  - m = 0 : (−3)¹, (−1)³, (1)³, (3)¹ ;
  - m = 1 et m = 5 : (±√6)², (0)⁴ ;
  - m = 2 et m = 4 : (±2)³, (0)² ;
  - m = 3 : (±√3)⁴, irréductibles.
* K26 : 26 sites, 48 arêtes, 24 carrés, χ = 2. ℤ/2 → +1. U(1) avec m = 1, 2, 3, 5 → −1, +1, −1, −1. À m = 1 : −3,504 (2), −2,244 (4), −0,742 (6, Σ\|χ\|² = 2).
* Tétraèdre, groupe T : π/2 → −1, (±√3)². 0 → +1. π → +1. 3π/2 → −1.
* **Précision sur Q3.** Pour les 24 rotations, U γ_i U† = +γ_{|R|(i)}, toujours avec le signe +. Sous le demi-tour x, γ_y ↦ +γ_y, alors qu'un vecteur donnerait −γ_y. Les γ_i portent donc la représentation de permutation des axes, A₁ (= H/√3) ⊕ E (= γ⊥). Conséquences :
  - « covariants » est exact, mais ce **ne sont pas des Γ vectoriels** ;
  - la suggestion de V180 §4.2 (remplacer les Γ empruntés par les γ_i) ne pouvait pas marcher ;
  - l'échec de D-P5 (V183) découle de Σ_d γ⊥_d = 0.

### 2.2 V183 (rv2, rv3)
* **Arbre à deux niveaux**, avec les poids 1,0 / 0,83 / 0,17 r² / 0,61 / 0,29 / 0,93 / 0,57 / 0,41 / 1,3, une jauge aléatoire et des bases aléatoires :
  - C0 : +1 ;
  - C1 : aucun relevé (23/24) ;
  - C2 : −1, niveaux bas tous de dimension 4 ;
  - C3 : aucun relevé, bloc racine de 16 ;
  - C4 : aucun relevé.
  Les cinq résultats sont reproduits.
* **Cl(5)**. Avec twist, dim A₁ = 4 (non commutative), dim T₁ = 8. Résultats de la recherche (200 départs) :
  - 141 solutions de Clifford trouvées ;
  - pour toutes, Γ₆ = ±1 (les deux signes apparaissent) et Γ₁Γ₂Γ₃ est invariant ;
  - commutant de dimension 4, sur lequel les rotations agissent comme End(E).
  Sans twist : aucune solution (meilleur résidu 7,5).
* **Démonstration de D-P2 à D-P4 et D-P6** (vérifiée par le calcul) :
  - On a ℂ⁸ = G₃/₂ ⊗ ℂ²_mult et G₃/₂ = E₁/₂ ⊗ E. Donc ℂ⁸ = E₁/₂ ⊗ E ⊗ ℂ²_mult, avec H/√3 = 1 ⊗ 1 ⊗ σ_z.
  - Les choix Γ_i = σ_i ⊗ 1 ⊗ σ_y, Γ₄ = 1 ⊗ 1 ⊗ σ_x, Γ₅ = 1 ⊗ 1 ⊗ σ_z conviennent, ce qui prouve l'existence.
  - Pour toute solution, Γ₁Γ₂Γ₃ est invariant (det R = 1), donc Γ₆ ∈ A₁ = 1 ⊗ M₂. Γ₆ commute avec Γ₄ et Γ₅, qui engendrent M₂, donc Γ₆ = ±1.
  - Le commutant de l'algèbre engendrée est End(E). Le doublet de goût est donc le facteur E de G₃/₂ = E₁/₂ ⊗ E.
  - L'argument de V182 §2.2 est juste ; il lui manque seulement l'étape « Γ₁Γ₂Γ₃ invariant ».

### 2.3 V193 (rv4, rv8)
* **Spectre exact** : E(k)² = Σ_d (t_f² + t_c² + 2 t_f t_c cos k_d), dégénéré 4 + 4 en tout k (écart 6·10⁻¹⁵ sur 200 points). En K + q, cela donne E² = 3(t_f − t_c)² + 2 t_f t_c Σ(1 − cos q_d). Le gap ne se ferme qu'à t_f = t_c, en K.
* **Décomposition** : elle est exacte pour q arbitraire et pour des (t_f, t_c) aléatoires. Les six matrices {γ_d, α_d} anticommutent deux à deux : c'est exactement la forme de Benalcazar–Bernevig–Hughes, h(k) = Σ_d [(γ + λ cos k_d) Γ_{2d} + λ sin k_d Γ_{2d−1}], avec γ = t_f et λ = t_c. Les auteurs décrivent un modèle octupolaire où chaque plaquette porte un flux π ([arXiv:1611.07987](https://arxiv.org/abs/1611.07987)). La condition |γ/λ| < 1 y est énoncée explicitement pour le quadrupôle ; pour l'octupôle, elle se démontre ici directement par la séparabilité (§2.4).
* **Γ₄** : l'unique Γ₄ ∈ A₁ (au signe près) qui anticommute avec Γ₅ et les α_d est Γ₄ = Γ₅α₁α₂α₃. Il **commute** avec le reste γ⊥ (ordre q²). Pour la piste B, ce reste sort donc de l'algèbre de Wilson–Dirac à 5 Γ.
* **Pas de symétrie de goût exacte** : le commutant de {H(k)} pour tous les k est de dimension 1, car Cl(6) ≅ M₈. La dégénérescence 4 vient de la structure de Clifford, pas d'une symétrie de goût. Le « doublet de goût » n'est défini qu'en q → 0.

### 2.4 V195 (rv5) : séparabilité exacte
* En jauge de Kawamoto–Smit, sur un amas ouvert, **H₃D = h ⊗ 1 ⊗ 1 + Z ⊗ h ⊗ 1 + Z ⊗ Z ⊗ h**, où h est la chaîne SSH ouverte (sauts t_f et t_c alternés) et Z = (−1)^x. Donc H² = Σ_d h_d² et E² = e_a² + e_b² + e_c² (vérifié opérateur par opérateur).
* **Énergie des modes de coin** : E = √3 δ_n **exactement**, où δ_n est l'énergie du mode de bord d'une chaîne SSH ouverte de 2n sites. Rapport E / [√3(1 − ε²)ε^n] :

| ε | n = 2 | n = 3 | n = 4 | n = 5 | n = 8 |
|---|---|---|---|---|---|
| 0,3 | 1,015 | 1,002 | 1,0002 | 1,0000 | — |
| 0,5 | 1,105 | 1,034 | 1,011 | 1,003 | 1,0001 |
| 0,7 | 1,44 | 1,23 | 1,13 | 1,07 | 1,012 |

  La loi n'est vraie qu'asymptotiquement en n.
* **Octet effectif** : dans la base produit des modes de bord localisés, H_eff² ∝ 1. Les seuls couplages non nuls relient des coins à distance de Hamming 1. Le produit des sauts autour de chaque face vaut −1 (flux π). |t_f′| = δ_n (0,0970 pour n = 3, ε = 0,5).
* **Chiralité de sous-réseau** : chaque mode de coin localisé est à 1,0000 sur le sous-réseau de son site de coin, soit 4 coins pairs et 4 impairs. E1-c est donc vérifié.
* **États dans le gap** (amas 8³) : 8 pour ε = 0,1 ; 0,2 ; 0,3 ; 0,5. Aucun pour 1,25 et 2. À ε = 0,8, les modes de coin existent, mais ils ne sont pas séparés à n = 4 (taille finie).

### 2.5 V197–V198 (rv6)
* **V197** : chaîne uniforme ε, 1, ε, 1, ε, 1, ε → un seul couple de modes de bord. En 3D, le 64ᵉ et le 65ᵉ niveaux coïncident (1,1538 et 1,1538 à ε = 0,2). L'explication de l'échec est juste.
* **V198** :
  - Le modèle hiérarchique est lui aussi séparable. La chaîne de 8 sites (ε, 1, ε, τ, ε, 1, ε) a, à basse énergie, quatre niveaux ±a, ±b.
  - La chaîne effective (x, y, x) a pour niveaux |e| = √(x² + y²/4) ± y/2. **Toute paire (a, b) est reproduite exactement**, avec x = √(ab) et y = a − b.
  - Mon ajustement 3D indépendant redonne t_f′ et t_c′ à 10⁻⁹ près, avec des résidus de 10⁻⁹ à 10⁻¹⁵. Le résidu < 5·10⁻⁴ de V198 n'est qu'une tolérance numérique.
* **Calcul à deux pas** (chaîne de 16 sites, équivalente par séparabilité à un amas de 16³ = 4096 sites), ε = 0,3, τ = τ* = 0,2990 :

| t(4) | t′(1)/t′(2) | t′(3)/t′(2) | résidu | autosimilaire ? |
|---|---|---|---|---|
| τ² = 0,0894 | 0,3025 | 0,2971 | 1,7·10⁻³ | oui |
| 1 (profil « alterné ») | 0,3042 | **3,02** | 1,1·10⁻⁴ | non |
| τ | 0,3288 | 1,013 | 1,5·10⁻³ | non |

  La valeur autosimilaire est t(4) = 0,0900 ≈ τ². Le profil autosimilaire est donc **t(h) ≈ ε, 1, ε, ε², ε³ …**, c'est-à-dire t(h) ≈ ε^|h−2| : une décroissance géométrique avec la profondeur, le lien entre frères faisant exception. Ce n'est pas une alternance.

### 2.6 V201 (rv7)
* Champ ℤ₃ en jauge axiale (a_x = 0) : toutes les plaquettes sont vérifiées.
* Octet : φ = 0 → (±√3)⁴, H² = 3 ; φ = ±1 → (±√6)², (0)⁴, H² ∈ {0, 6}. Commutateur −1 dans les trois cas.
* Modes nuls (φ = 1) : 4, 8, 12, 16 pour n = 1…4, quel que soit ε ∈ {0,3 ; 0,5 ; 1,0}. Leur poids sur les octets de coin vaut 0,34, pour une fraction de sites de 0,30 : ils ne sont pas localisés aux coins.

---

## 3. Erreurs trouvées

### Graves
* **G1 — V198 §2, repris en V201 §4 : le « profil alterné » est faux.**
  - Le calcul à deux pas (§2.5) montre que l'autosimilarité impose t(4) ≈ τ² ≈ ε². Le profil est faible (h = 1), fort (h = 2), puis décroît géométriquement : ε, ε², …
  - Avec t(4) = 1 (alterné), t′(3)/t′(2) = 3,0 au lieu de τ.
  - La piste « une règle de contact qui distingue la **parité** de la profondeur » (V198 §2 et §4.2) part donc d'une prémisse fausse. Il faudrait chercher une amplitude qui **décroît** avec la profondeur de l'ancêtre commun, et non une amplitude qui alterne.
* **G2 — V198 C2 « oui, exactement » n'est pas un test.**
  - Par la séparabilité, les 64 niveaux se réduisent à deux nombres (a, b), que la chaîne effective à deux paramètres reproduit toujours exactement. C2 ne pouvait pas échouer.
  - « L'écorce d'une échelle est le feuillage de la suivante, **exactement** » est donc, à un pas, une propriété triviale de la chaîne SSH de 4 sites, et non une preuve d'autosimilarité.
  - Le contenu non trivial n'apparaît qu'à deux pas : il est vérifié à ~10⁻³ près avec le profil géométrique, pas avec le profil alterné.
  - L'énoncé « acquis » de V201 §4 doit être reformulé.

### Moyennes
* **M1 — V193 §1 : la fenêtre de Kaplan 1/3 < t_f/t_c < 1 est un artefact.**
  - Le modèle n'a pas de doubleurs : le gap ne se ferme qu'en t_f = t_c, en K.
  - Le critère −2r_W < m < 0 a été appliqué en ignorant le reste t_c Σ(1 − cos q)γ⊥, qui est du même ordre que le terme de Wilson. Par exemple, en q = (π, 0, 0), ce point n'est pas un « doubleur » : E² = (t_f + t_c)² + 2(t_f − t_c)².
  - Les modes de coin existent à ε = 0,1 et 0,2 (V195 et rv4).
  - La phase topologique, au sens de BBH, est toute la région 0 < t_f/t_c < 1, pour un amas fait d'octets complets. Couper à travers les octets inverse les rôles.
  - La note avertissait « à ne pas citer comme fenêtre établie », mais V194 reprend « la phase t_f/t_c < 1 (celle que V193 désignait… par le critère de Kaplan) ».
* **M2 — V193 §1 et §3.1 : « frères → masse et Wilson » est faux.**
  - Les sauts entre frères ne dépendent pas de k : ils ne donnent que √3 t_f Γ₅.
  - Les cousins donnent −√3 t_c Γ₅, **le terme de Wilson (r_W = t_c/√3)**, la partie cinétique et le reste γ⊥.
  - Énoncé correct : « la masse m = √3(t_f − t_c) dépend des deux ; seuls les frères n'entrent **que** dans la masse ». La conséquence de V199 §1 (une texture par secteur ne peut toucher que t_f sans changer la cinétique ni le Wilson) reste vraie, mais pour cette raison-là.
* **M3 — V195 : « exactement E = ±√3(1 − ε²)ε^n », « rapport constant à 4 chiffres pour n = 3, 4, 5 ».**
  - C'est faux à ε = 0,5 : 3,4 %, 1,1 %, 0,3 % d'écart.
  - La loi exacte est E = √3 δ_n (mode de bord SSH). La formule n'en est que l'asymptote, déjà connue pour la chaîne SSH : ce n'est pas un « fait nouveau ».
  - L'inférence « ce spectre ⇔ cube à flux π » est insuffisante : le spectre ne fixe pas la connectivité. La conclusion est néanmoins vraie (§2.4).
* **M4 — V183 §1 et §4 : « un seul twist par branche, et seulement au niveau le plus fin ».**
  - Les cinq configurations s'expliquent toutes par une règle plus simple, déjà contenue dans V182 §1.1 : **deux secteurs couplés doivent porter la même classe projective**.
  - En effet, si les relevés sont uniques à une phase près, le commutateur est un relevé de l'identité, donc un scalaire. Un couplage entre un secteur de classe +1 et un secteur de classe −1 est donc nécessairement nul.
  - C4 (feuilles de classe +1, niveau 1 tordu de classe −1) était donc prévisible : la prédiction « −1 » contredisait la propre logique de V182.
  - La lecture « seulement au niveau le plus fin » n'est établie que si les nœuds sont scalaires ou si leurs multiplets sont induits d'en dessous. Un twist au niveau 1 avec des nœuds porteurs d'un multiplet projectif interne n'a pas été testé.
  - V183 §1 écrit bien « [E, pour cette famille] », mais V183 §4 le répète sans réserve.
* **M5 — V183 §3 : « le bord natif est chiral » (en gras).**
  - L'équivalence unitaire cellule par cellule est juste : un Cl(5) avec Γ₆ = ±1 sur ℂ⁸ est un module isotypique.
  - Mais la conclusion dépend (i) de la réserve octet/quadruplet écrite dans la même section, et (ii) du flot −q² de V181, qui est **diabatique** et non topologique à taille finie.
  - À rétrograder en [C].
* **M6 — V201 §2 : décision d'arrêt.**
  - La règle de V200 §3 prévoyait deux cas : arrêter si le filtre « ne filtre rien » ou s'il détruit les modes de coin « sans rien mettre à la place » ; continuer s'il apparaît « une loi de hiérarchie différente ». Le résultat (4n modes nuls) ne tombe littéralement dans aucun des deux : la décision est une **interprétation**, qu'il faut signaler comme telle.
  - Elle est défendable par un argument plus fort que celui donné : tout flux ℤ₃ non nul à travers les faces de l'octet y détruit H² ∝ 1 (lemme : H² = 3·1 ⇔ π sur les 6 faces), donc détruit le Dirac natif.
  - En revanche, la conclusion « la géométrie (écorce, fissures, renormalisation, filtre par le centre) ne sépare pas les secteurs sans casser le Dirac » est trop large : seul le fond dense a été testé. L'écorce texturée (option 1 de V199) et les filtres dilués ne l'ont pas été (V201 §3 le reconnaît).
* **M7 — V193 §2 : la séparation espace/fibre présentée comme « contrainte dérivée [T + C] ».**
  - La partie [T] se limite à A²\|_R = 3I sous twist. Le reste dépend d'hypothèses du projet que je n'ai pas relues : facteur 3 de Georgi–Jarlskog = terme J − I (H-coin), R/L ↔ I/II, ordre un.
  - V192 l'étiquetait [C] ; la montée en « contrainte dérivée » n'est pas justifiée.
* **M8 — Horodatage** (pertinent pour la règle « critères avant calcul ») :
  - V180 est daté de 15 h 45, mais V182 (15 h 37) cite les résultats de V180, et JJ décide sur V180 §4 à 15 h 32. L'heure de V180 est donc fausse.
  - V194 (16 h 24) cite V193, daté de 16 h 25.
  - Des notes denses (V194 : 45 lignes) sont datées une minute après la question de JJ.
  - Les heures sont visiblement approximatives ou reconstruites. L'antériorité des critères ne peut pas être certifiée par elles. Dans les cas relus, les critères sont cohérents avec un ordre « avant » (les prédictions fausses sont honnêtement signalées : C4, Q2-b, R1, E1-a), mais l'horodatage ne le prouve pas.

### Mineures
* **m1** — V180 §2 : « nombre impair de faces frustrées ». C'est la charge m = flux total/2π qui doit être impaire. Le cube a 6 faces à π.
* **m2** — V179 et V180 : « flux U(1) d'ordre 24 » pour K26. m = 3 (π/4 par carré, ordre 8) suffit.
* **m3** — V180 Q3 : covariance avec le signe +, en représentation A₁ ⊕ E et non T₁ (§2.1). V180 §4.2 est à corriger.
* **m4** — V182 §2.2 : étape manquante (Γ₁Γ₂Γ₃ invariant), et « scalaire » est employé tantôt pour « invariant », tantôt pour « multiple de 1 ».
* **m5** — V193 : « v = √(t_f t_c) » est le coefficient de q² dans E². La vitesse de Dirac au sens de l'opérateur est t_c, et v² = t_c² + m r_W.
* **m6** — V193 : « brisure de goût » par γ⊥. Il n'existe aucune symétrie de goût exacte (commutant trivial), mais la dégénérescence 4 subsiste partout. Γ₄ est unique, Γ₄ = ±Γ₅α₁α₂α₃, et commute avec γ⊥ : à garder en tête pour la piste B.
* **m7** — V195 §2 : « couplages … le long de la diagonale (3n) ». Ils sont nuls exactement dans ce modèle (séparabilité) : seuls les coins voisins se couplent.
* **m8** — V198 :
  - w est le poids 1D du site d'extrémité. Le poids 3D sur le site de coin vaut w³ (0,90 et non 0,96 à ε = 0,2).
  - τ*(0,5) = 0,491, et non 0,485, qui est la valeur au premier ordre.
  - « τ* ≈ ε, remarquable » s'explique : τ* = u(1 + u²) avec u = δ/ε, soit ε + O(ε⁵).
* **m9** — V199 §2 : « R = ε/τ ». Non dérivé : il confond le rapport ε′ avec le facteur d'échelle des énergies (≈ wτ) et ignore t(4). Sans conséquence, car (1) n'a pas été testé.
* **m10** — V201 §3 : « filtre ℤ₂ exclu d'avance ». Seule la version dense l'est.
* **m11** — V194 : l'identification avec BBH est juste, mais elle dépend de la convention de bord : la phase t_f/t_c < 1 est la phase à coins pour une troncature par octets complets.

---

## 4. Ce qui n'a pas pu être vérifié
* **V181** en entier : les scripts de V166 et V167 ne sont pas fournis.
* Les renvois internes au projet : V105, V107, V117, V121, V141, V142, V145–V147, V151, V158, V169, V173–V178. Cela concerne en particulier le facteur 3 de Georgi–Jarlskog, l'ordre un, la lecture Majorana de Γ₅ (V192 §1.3) et la lecture « goût E = doublet faible » (V192 §1.4).
* La valeur du flot −2 conditionnelle à V151 ⊗ 1₂ (V183 §3) : aucun calcul de flot n'a été refait.
* La condition |γ/λ| < 1 pour l'octupôle dans le texte même de BBH : la page source n'a été lue qu'en partie ; la condition est citée pour le quadrupôle. Elle est démontrée ici indépendamment (§2.4).
* Le spectre de Pearson–Bellissard (V194 §2.1, « cité de mémoire ») : seuls les caractères des 7 ondelettes ont été vérifiés.
* L'existence d'une configuration « twist au niveau 1 avec multiplets internes » (M4) : non construite.
* L'heure réelle de rédaction des critères.

---

## 5. Corrections recommandées (courtes)
1. **V198 et V201 §4** : remplacer « profil alterné » par « t(h) ≈ ε^|h−2| (décroissance géométrique au-delà des cousins) ». Remplacer « exactement » par « exact à un pas par séparabilité (non probant) ; vérifié à deux pas à 10⁻³ près ».
2. **V193** : retirer la fenêtre 1/3 < t_f/t_c < 1. Corriger « frères → Wilson ». Ajouter le spectre exact et la séparabilité.
3. **V195** : loi exacte E = √3 δ_n et asymptote √3(1 − ε²)ε^n. E1-c : tranché (oui).
4. **V183 §1 et §4** : énoncer la règle de concordance des classes. Ramener « un seul twist par branche » au rang de [C], limité à cette famille. **V183 §3** : « chiral » en [C].
5. **V201 §2** : signaler que la décision interprète la règle, et restreindre la conclusion au fond ℤ₃ dense.
6. **V180** : corriger « nombre impair de faces », « ordre 24 », et l'heure (15 h 45).

---

## 6. Résumé
Les calculs de spin sont tous reproduits par un code indépendant : classe −1 de l'octet tordu, règle (−1)^m pour m = 0…5, K26, tétraèdre, sauts γ_i, configurations C0–C4 de l'arbre. Les γ_i se transforment toutefois comme la permutation des axes (A₁ ⊕ E), pas comme un vecteur.
Le Cl(5) covariant de l'octet tordu est confirmé et démontré : ℂ⁸ = E₁/₂ ⊗ E ⊗ ℂ², Γ₆ = ±1 pour toute solution, aucune solution sans twist, Γ₅ = H/√3 possible. L'argument papier est juste, à une étape près.
Le résultat C4 (prédiction fausse) découlait déjà de la règle « classes égales sur les secteurs couplés ». La lecture « un seul twist, au niveau le plus fin » est trop forte.
Le Dirac natif de V193 est exactement le modèle octupolaire de Benalcazar–Bernevig–Hughes. La décomposition est exacte pour tout q, et la vitesse vaut √(t_f t_c).
Deux erreurs dans V193 : la fenêtre de Kaplan 1/3 < t_f/t_c < 1 est un artefact (phase = tout t_f/t_c < 1), et le terme de Wilson vient des cousins, pas des frères.
Fait structurel non vu dans les notes : le modèle est séparable (H = h ⊗ 1 ⊗ 1 + Z ⊗ h ⊗ 1 + Z ⊗ Z ⊗ h). Tous les résultats de coins et de renormalisation se réduisent à des chaînes SSH.
Conséquences pour V195 : 8 modes de coin, octet tordu effectif et E1-c sont confirmés, mais E = √3(1 − ε²)ε^n n'est qu'asymptotique.
Conséquences pour V198 : le test C2 (« exactement ») ne pouvait pas échouer. τ* ≈ ε s'explique analytiquement (τ* = ε + O(ε⁵)).
Le « profil alterné » est faux : le calcul à deux pas impose t(4) ≈ τ², soit t(h) ≈ ε^|h−2|.
V201 est reproduit (cristal ℤ₃ cohérent, H² non ∝ 1, 4n modes nuls indépendants de ε). L'arrêt est défendable, mais il interprète la règle de V200, et sa conclusion générale dépasse ce qui a été testé.
Les heures sont incohérentes à plusieurs endroits (V180, V194) et ne certifient pas l'antériorité des critères.
Non vérifié : V181, les renvois internes au projet (Georgi–Jarlskog, ordre un, doublet faible) et la valeur conditionnelle du flot de V183 §3.
