# Relecture indépendante de V167 (flot spectral à N = 4) et de sa méthode (V166)

Date : 2026-11-12. Relecteur : Claude (agent de relecture). Je n'ai pas écrit les calculs relus.
Étiquettes : **[E]** recalculé par moi, avec mes propres scripts ; **[T]** vérifié par lecture du code ou par preuve ; **[C]** lecture ou jugement.
Scripts et journaux : `/home/claude/fht/relecture3/`. Aucun fichier existant n'a été modifié. Le modèle de référence reste `V162_FHT_R3_jauge.py`, chargé tel quel (comme le fait V166). Tout le reste est réécrit :
* `rv3_core.py` : chargement du modèle ; reconstruction de Fourier **indépendante** (moindres carrés sur des θ aléatoires, et non la DFT de V166) ;
* `rv3_check_fourier.py` (N = 3) et `rv3_check_fourier_N4.py` (N = 4) : point 1 ;
* `rv3_sym*.py` : point 3 ;
* `rv3_scan.py` : balayages fins en θ, avec E, pente et profil par couche des 16 à 24 états les plus proches de 0. Le shift-invert y utilise l'ordre `MMD_AT_PLUS_A`, environ 4 fois plus rapide et 3 fois moins gourmand en mémoire que l'ordre par défaut ;
* `rv3_an_N3.py` et `rv3_an_N4.py` : analyses (sommes de V166 reproduites, trapèze, trois fenêtres) ;
* `rv3_an_exact.py` : intégrale de la fenêtre carrée **sans erreur de discontinuité**. Chaque état est apparié d'un point au suivant, E(θ) est interpolé linéairement, et l'on intègre exactement l'indicatrice \|E\| < ε ;
* `rv3_degen.py` : doublets presque dégénérés.

Journaux : `rv3_fourier_*.log`, `rv3_sym_spectre_N3.log`, `an_N3_q1.log`, `an_N3_q2.log`, `an_N3_exact.log`, `an_N4_q1.log`, `an_N4_q2.log`, `an_N4_q1_exact.log`, `an_N4_q2_exact.log`, `scan_*.log` (données `scan_*.npz`).

---

## Tableau des verdicts

| point | verdict |
|---|---|
| 1. Reconstruction de Fourier (fréquences k/2L, 4qL + 1 points) et dérivée analytique | **confirmée, exacte et complète [E/T]** |
| 2a. « −3,99 hors θ ≡ 0 mod π » : le mécanisme invoqué | **en partie juste [E]**. La paire bord–racine hybridée à θ = 0 est réelle. « Pentes très grandes » et « estimateur mal défini » sont faux : le défaut est un **sous-échantillonnage**, sur une grille commensurable avec la fenêtre |
| 2b. La sélection | **tri commode [E]**. Sa coïncidence avec −4 est **fortuite** : elle repose sur des états de bord à \|E\| = ε − 7·10⁻⁵ exactement sur des points de la grille. Le total reste à −1,96 |
| 2c. Convergence vers −4 | **oui, avec la grille raffinée [E], à 2 % près**. q = 2 : couche fine **−3,93**, grossières +3,93, total **0,000** (fenêtre ε de V166, intégrale exacte en θ). −3,87 avec ε/2 et avec une fenêtre lisse. Le −3,99 de V167 est donc lui-même biaisé de 0,06 |
| 2d. q = 1 convergé | **−0,998**, grossières +0,997, total 0,000 [E]. Les mêmes pièges de grille expliquent les −0,92 et −1,00 de V167 |
| 3. Égalité des sous-grilles 1/4 et 3/4 par θ → −θ | **confirmée [E]**. La relation exacte est E(−θ) = −E(θ), avec pentes et profils égaux, et non une isospectralité. Opérateur explicite trouvé à N = 2 seulement |
| 4. Verdicts F-P1, F-P2, F-P3 | **mal fondés dans V167 [E/C]**, en partie justes après convergence (§4) |
| 5. Autres erreurs | référence N = 3 erronée (−3,21 vient du transport (a)) ; +0,78 au lieu de +0,80 ; « loi en q² à 1 % » trop fort ; flot strict nul (§5) |

---

## 1. Reconstruction de Fourier et dérivée analytique

### 1.1 Lecture du code [T]
* Dans V162, θ n'entre que par les liens z fins u_z = e^{iqθ/L}.
* Avec le transport (b) :
  - un lien vertical porte le chemin x, y, z dans le bloc, donc une phase e^{iqθ d_{n+1} c_z / L} ;
  - un lien latéral z grossier porte e^{iqθ d_n/L} fois une moyenne indépendante de θ (règle r1 : moyenne non projetée) ;
  - dans les liens latéraux x et y, les segments z de P(x, c) et de P(x + e_j, c) se compensent.
* Les seules fréquences, en unités de 1/2L, sont donc **k ∈ {0, ±2q·2^j, j = 0..N}**, soit 11 composantes à N = 4. C'est le « 11 non nulles » des journaux.
* La DFT sur 4qL + 1 points, sur la période 4πL, est exacte. Elle est même suréchantillonnée : seules 11 des 4qL + 1 fréquences servent. La dérivée Σ (ik/2L) e^{iθk/2L} C_k est correcte.

### 1.2 Vérification [E] (`rv3_fourier_N3.log`, `rv3_fourier_N4_q1.log`, `rv3_fourier_N4_q2.log`)
| contrôle | N = 3, q = 1 | N = 3, q = 2 | N = 4, q = 1 | N = 4, q = 2 |
|---|---|---|---|---|
| ajustement par moindres carrés, bande élargie (\|k\| ≤ 4qL à N = 3) ou 11 fréquences sur 30 θ aléatoires (N = 4) : résidu | 1,4·10⁻¹⁴ | 2,2·10⁻¹⁴ | 7,8·10⁻¹⁵ | 2,3·10⁻¹⁴ |
| plus grande composante hors {0, ±2q·2^j} | 6·10⁻¹⁵ | 6·10⁻¹⁴ | 4·10⁻¹⁵ (DFT V166) | 1,2·10⁻¹⁴ (DFT V166) |
| C_k de la DFT de V166 − mes C_k | 7·10⁻¹⁵ | 6·10⁻¹⁴ | 3·10⁻¹⁵ | 7·10⁻¹⁵ |
| H reconstruit − H(V162) (θ = 0, π, π/q, 1000π + 0,3, aléatoires) | 1,3·10⁻¹⁴ | 2,4·10⁻¹³ | 7·10⁻¹⁵ | 3·10⁻¹⁴ |
| ∂_θH analytique − différences à 5 points (h = 10⁻³) | 2·10⁻¹⁰ | 3·10⁻¹⁰ | 1,1·10⁻¹¹ | 1,1·10⁻¹¹ |

**Contrôle de bout en bout [E].** Ma chaîne diffère de celle de V166 par l'ajustement de Fourier, l'ordre de factorisation et l'usage des symétries du §3. Elle reproduit pourtant :
* le flot de validation à N = 3 (48 pas) **aux 6 décimales** : [0,008645 ; 0,456719 ; 0,483716 ; −0,917643] ;
* à N = 4, toutes les sommes de V167 au millième : grille complète, pairs et impairs, θ ∈ {0, π}, sous-grilles 1/4 et 3/4, pour q = 1 et q = 2 (`an_N4_q*.log`).

Les chiffres bruts de V167 sont donc bien ceux que produit le code. Ce qui suit porte sur leur interprétation.

## 2. La sélection « hors θ ≡ 0 mod π » et la convergence

### 2.1 Ce qui se passe près de θ = 0 [E] (N = 4, q = 2 ; `an_N4_q2.log`, balayages fins sur [0 ; 0,01] au pas 5·10⁻⁴ et sur [0 ; 0,2] au pas 5·10⁻³)
Notations : w_fine est le poids dans la couche fine, ε = π/(2L) la demi-largeur de la fenêtre, f_n(θ) l'intégrande de la fenêtre pour la couche n.
* La fenêtre contient 4 états : deux doublets presque dégénérés (écart ~10⁻⁵ ; D et P_n y sont ∝ 1 à 10⁻³ près, donc l'estimateur ne dépend pas de la base, cf. `rv3_degen.py`) :
  - **bord** : w_fine = 1,00, pente −0,125 = −q/L ;
  - **grossiers** (couches 1–3, surtout 2) : w_fine = 0,00, pente **+0,51**.
* Croisement évité à θ = 0, avec un écart 2 × 9·10⁻⁴ = **1,8·10⁻³**. La zone d'hybridation est |θ| ≲ 0,004. **Exactement à θ = 0** (point fixe de θ → −θ, §3), les 4 états sont mélangés 50/50, avec la même pente +0,194. On a alors f_fine(0) = +1,98, contre −1,26 juste à côté (θ = 0,01 : −1,02 ; θ = 0,05 : −1,25).
* Les états grossiers restent dans la fenêtre pour |θ| < **0,193** : environ 2 pas de la grille à 64 pas de chaque côté. L'intervalle n'est donc pas « étroit » ; ce qui est étroit, c'est la zone d'hybridation.
* Les pentes à θ = 0 valent 0,19. **Rien n'est « très grand »**. L'estimateur à fenêtre est **bien défini**, avec un intégrande borné.

### 2.2 Pourquoi la grille de 64 pas se trompe [E]
* **Grille paire** (pas 2π/32).
  - Les points θ = ±0,196 tombent **juste hors de la fenêtre** pour les états grossiers (\|E\| − ε = +1,7·10⁻³). La compensation n'y est vue qu'à θ = 0, où la moitié est attribuée à la couche fine.
  - Les états de bord sortent de la fenêtre à θ ≈ π/(2q). En effet, E_bord ≈ −qθ/L et ε = π/(2L), donc \|E\| = ε pour θ = π/(2q), qui est **un point de la grille** (π/4 pour q = 2, π/2 pour q = 1).
  - À θ = π/4, le doublet de bord est à **\|E\| = ε − 6,8·10⁻⁵** : dedans de justesse. Ces 4 points (π/4 mod π/2) apportent −0,50 à la couche fine.
* **Grille impaire** (décalée d'un demi-pas). Elle évite les deux pièges : fine −1,99, total +0,04 par demi-grille.
* **Conséquence.**
  - « Hors θ ≡ 0 » retire bien la contribution mal attribuée (+0,39). Mais il retire aussi la contribution régulière de ces points (−0,12 chacun), et la compensation reste perdue : le **total vaut −1,96**.
  - Le −3,987 ne tombe sur −4 que parce que les états de bord à θ = π/4 mod π/2 sont dedans de 7·10⁻⁵. Un décalage de ε de 10⁻⁴ donnerait environ −3,49.
  - **La sélection n'est pas une correction ; son accord avec −4 est fortuit.** De plus, la valeur convergée est −3,93 (§2.3) : −3,987 n'est donc même pas le bon chiffre à 0,06 près.

### 2.3 Intégrale convergée, N = 4, q = 2 [E]
65 points uniformes sur [0, π/2] (pas 0,0245) et les deux zooms. On utilise f(−θ) = f(θ) et la période π (§3), puis le trapèze.

| estimateur | n = 0 | n = 1 | n = 2 | n = 3 | n = 4 (fine) | total |
|---|---|---|---|---|---|---|
| V166, 64 pas (reproduit) | +0,003 | +0,259 | +1,825 | +0,330 | −3,599 | −1,181 |
| V166, 128 pas (Riemann) | +0,004 | +0,346 | +2,443 | +0,444 | −3,790 | −0,552 |
| fenêtre ε, trapèze sur les points (indicatif) | +0,005 | +0,418 | +2,960 | +0,540 | −3,990 | −0,066 |
| **fenêtre ε, intégrale exacte en θ** (`an_N4_q2_exact.log`) | +0,005 | +0,419 | +2,965 | +0,541 | **−3,930** | **−0,000** |
| fenêtre ε/2, intégrale exacte | +0,005 | +0,416 | +2,924 | +0,528 | −3,873 | −0,000 |
| fenêtre lisse cos², demi-largeur ε (trapèze, intégrande continu) | +0,005 | +0,416 | +2,924 | +0,529 | −3,874 | −0,001 |
| flot brut, couche fine (fenêtre ε, exact) | | | | | 4,001 | — |

* Le trapèze sur la fenêtre carrée est biaisé, et de nouveau par les états de bord limites à θ = π/4 (\|E\| = ε − 7·10⁻⁵, sur un point de mon propre balayage). L'intégration exacte corrige ce biais : les totaux deviennent nuls à 10⁻³ près, comme ils doivent l'être, puisque le flot total d'une famille de matrices hermitiennes de dimension finie est nul sur un cycle. Elle est aussi cohérente avec la fenêtre lisse. Je retiens **−3,93** (fenêtre ε) et −3,87 (ε/2, et fenêtre lisse).
* Les deux autres lignes du tableau (64 et 128 pas) sont des sommes de Riemann sur la grille de V166.
* À 128 pas, la grille n'a toujours pas convergé : le point θ = 0 y pèse encore dθ. Pour converger, il faut une grille **décalée** (sans θ ≡ 0 mod π/q) ou une grille raffinée.
* **Dépendance à l'estimateur.** La zone d'hybridation enlève à la couche fine une part qui croît quand ε diminue : 0,07 (ε), 0,13 (ε/2). Le « −4 » est un flot **diabatique**, valable tant que ε ≫ écart (ici ε / écart ≈ 55). Ce n'est pas un invariant.
* **Flot strict** (comptage des passages par E = 0, suivi adiabatique) : **nul couche par couche**. Près de θ ≡ 0 mod π, le doublet inférieur reste à E < 0 (\|E\| min 7,2·10⁻⁴), et ailleurs aucun niveau ne passe par 0 ; de même pour q = 1 (\|E\| min 6·10⁻⁵). C'est la même situation qu'à N = 3 (relecture précédente, §1b), avec un écart 36 fois plus petit.

### 2.4 Référence N = 3, transport (b/r1) [E] (`an_N3_q1.log`, `an_N3_q2.log` ; 513 points par demi-période, diagonalisation creuse)
| | couche fine, ε | fine, ε/2 | total, ε | écart du croisement évité près de θ = 0 |
|---|---|---|---|---|
| N = 3, q = 1 | −0,918 | −0,868 | 0,000 | 1,2·10⁻² |
| N = 4, q = 1 | −0,998 | −0,996 | 0,000 | ≈ 1,8·10⁻⁴ |
| N = 3, q = 2 | **−2,521** | −1,773 | 0,000 | 6,5·10⁻² |
| N = 4, q = 2 | −3,930 | −3,873 | 0,000 | 1,8·10⁻³ |

(Intégrales exactes en θ : `an_N3_exact.log`, `an_N4_q*_exact.log`.)

* À N = 3, la même sélection « hors θ ≡ 0 » ne donne pas −4 : −2,57 à 64 pas.
* Entre N = 3 et N = 4, le passage vers −q² est réel, et la dépendance en ε diminue. C'est bien la taille qui était en cause, mais pas avec le chiffre cité par V167 (§5).

### 2.5 N = 4, q = 1 [E] (`an_N4_q1.log`)
Balayages : 65 points sur [0, π] (pas 0,049), puis [0 ; 0,2] au pas 0,01 et [0 ; 0,002] au pas 10⁻⁴. On utilise f(−θ) = f(θ) et la période 2π.
* **Structure.**
  - Un état de bord : w_fine = 1,00, pente −0,0625 = −1/L.
  - Un état grossier : pente +0,367, dans la fenêtre pour |θ| < 0,27.
  - Croisement évité à θ = 0 : E = ±9·10⁻⁵ ; zone d'hybridation |θ| ≲ 5·10⁻⁴ ; à θ = 0 exactement, mélange 50/50 et f_fine = +0,78.
  - L'état de bord sort de la fenêtre à θ = π/2, **sur un point de la grille paire**, avec \|E\| = ε − 1,2·10⁻⁴. C'est le même piège que pour q = 2.
* **Sommes de V166 reproduites** : pairs −0,423, impairs −0,499, complète −0,922, total −0,117.
* **Intégrale exacte** (`an_N4_q1_exact.log`) :
  - fenêtre ε : couche fine **−0,998** ; grossières +0,009 / +0,472 / +0,491 / +0,025 (somme +0,997) ; total 0,000 ;
  - fenêtre ε/2 : fine −0,996 ;
  - fenêtre lisse : fine −0,996 ;
  - flot brut de la couche fine : 1,000.
* La dépendance en ε a pratiquement disparu (écart ≈ 1,8·10⁻⁴, soit ε / écart ≈ 500).
* Le −0,922 de V167 est un artefact de grille, et le « −1,00 sur la demi-grille impaire » est juste, mais par chance.
* La prédiction de V166 (« plus proche de −1 qu'à N = 3 ») est vérifiée : on passe de −0,918 à −0,998.

## 3. Symétrie θ → −θ et égalité des sous-grilles 1/4 et 3/4

* **Relation spectrale [E]** (`rv3_sym_spectre_N3.log`, diagonalisation complète, N = 3, q = 1 et 2, θ = 0,37 et 2,2). Pour tous les 2340 états :
  - E_j(−θ) = −E_i(θ) à 10⁻¹⁴ près ;
  - **pentes égales** (et non opposées) à 4·10⁻¹² près (q = 1) et 2·10⁻¹⁰ (q = 2) ;
  - profils par couche égaux à 4·10⁻¹² près.
* **Conséquence.** Chaque terme w·E′ de la fenêtre est identique en θ et en −θ, donc f_n(−θ) = f_n(θ). L'indice j ↦ 64 − j envoie la sous-grille 1 mod 4 sur 3 mod 4. L'égalité exacte des journaux 1a et 1b est donc expliquée.
  - Une isospectralité H(−θ) ≅ H(θ) aurait donné des contributions **opposées**, donc un flot nul. V167 dit « symétrie θ → −θ » sans préciser le signe. C'est E → −E qui est en jeu.
* **À N = 4 [E, indirect].** Mes sommes reconstruites sur [0, π/2] en supposant la symétrie reproduisent au millième les sommes de V166, calculées sur [0, 2π) sans l'utiliser.
* **Périodicité [E].** H(θ + 2π/q) = Λ H(θ) Λ†, avec Λ une jauge diagonale, à 4·10⁻¹⁴ près (N = 3 et 4). Les points θ = 0 et π sont donc équivalents pour q = 2.
* **Opérateur [E/T, partiel].**
  - À N = 2, U = Λ·(R_z ⊗ Γ₃), où R_z est la réflexion en z, est exact (2·10⁻¹⁵) : U H(θ) U† = −H(−θ). Ainsi Γ₃ anticommute avec Γ₁, Γ₂, Γ₄ et Γ₅ (masse et Wilson changent de signe), et R_z renverse θ.
  - À N ≥ 3, cet opérateur est brisé, à 1,7·10⁻³ près, sur les seuls liens racine → niveau 1. Avec la règle r1, les liens z du niveau 1 ont \|u\| = 0,99 < 1, et le module du transport vertical dépend de l'enfant représentant, que R_z change.
  - Je n'ai trouvé aucun opérateur dans la famille testée : 16 opérations spatiales (dont x ↔ y), conjugaison complexe ou non, produits de Γ et rotation (1 + Γ₁Γ₂)/√2, jauge diagonale.
  - **L'égalité est donc vérifiée numériquement, mais non démontrée à N ≥ 3.**

## 4. Les verdicts de V167 face aux critères de V166

V166 ne fixe pas le nombre de pas. Le critère porte donc sur l'estimateur, et une grille raffinée est une mesure conforme à V166. La sélection de points, elle, ne l'est pas.

* **F-P1** (fine −1,00 ± 0,08 ; grossières +1,00 ± 0,1 ; |Σ| < 0,1).
  - À 64 pas : fine −0,922 (dedans, à 0,002 près) ; grossières **+0,805** (V167 écrit +0,78) ; total −0,117. **Deux conditions sur trois échouent.** « Partiellement vérifiée » est acceptable ; « non vérifiée au sens de V166 » serait plus exact.
  - Le bilan cite « −1,00 sur la demi-grille impaire » sans le marquer comme sélection a posteriori. La demi-grille paire donne −0,85.
  - **Après convergence : fine −0,998, grossières +0,997, total 0,000 : vérifiée [E]** (§2.5). V167 ne pouvait pas le savoir avec 64 pas.
* **F-P2** (fine −4,0 ± 0,4 ; |Σ| < 0,2).
  - À 64 pas : fine −3,599, **juste hors tolérance** (de 0,001) ; total −1,18. **Non vérifiée.**
  - « Vérifiée hors des deux points pathologiques » est une reclassification a posteriori, que le §2 de V167 reconnaît lui-même.
  - **Après convergence (fenêtre ε, intégrale exacte) : −3,93 et total 0,00 : vérifiée [E].** Avec ε/2 et avec la fenêtre lisse (−3,87), elle reste dans la tolérance.
  - La conclusion de V167 est donc juste. Sa démonstration ne l'est pas : le bon chiffre est obtenu par un tri dont l'accord est fortuit (§2.2).
* **F-P3** (flot brut de la couche fine ≈ |flot net| ; « ≈ » n'est pas chiffré dans V166).
  - À 64 pas : 1,075 contre 0,922 (+17 %) pour q = 1 ; 4,386 contre 3,599 (+22 %) pour q = 2. V167 cite la demi-grille impaire (2,00 contre 1,99) : nouvelle sélection.
  - Convergé : 4,00 contre 3,93 (+1,8 %) pour q = 2 ; 1,000 contre 0,998 pour q = 1. **Vérifiée.** Le flot inverse de la couche fine ne vient que de la zone d'hybridation.
* **Clause d'échec** : non déclenchée. C'est correct (−0,922 ∈ [−1,1 ; −0,9] ; −3,60 à moins de 0,6 de −4).
* **Honnêteté.** L'« écart d'honnêteté » du §2 de V167 est franc : la sélection est dite a posteriori, et −3,60 est donné comme mesure conforme. Mais le bilan (§0) met en gras −3,99 et « vérifiée », et le §2 range « −4 (loi en q², à 1 %) » dans les acquis [E]. La hiérarchie du texte contredit l'aveu. De plus, la justification par le mécanisme est en partie fausse (§2.1–2.2).

## 5. Autres erreurs, ambiguïtés, conclusions trop fortes
1. **Référence N = 3 erronée.** « À N = 3 (V164), q = 2 donnait −3,21 » : ce chiffre vient du transport **(a)**, à 128 pas (`V164_q2_a_*.log` : −1,535 − 1,672). Pour (b/r1), aucun calcul q = 2 à N = 3 n'existait. Je trouve **−2,52** (convergé, ε) et −1,77 (ε/2). La comparaison de V167 mélange deux transports.
2. **Arithmétique.** Compensation pour q = 1 : +0,805, et non +0,78.
3. **« Loi en q², à 1 % »** : le −3,99 de V167 est biaisé. Convergé, on a −3,93 (ε) et −3,87 (ε/2) : la loi en q² tient à 2–3 % pour q = 2, et à 0,2 % pour q = 1 (−0,998). L'écart à −4 est réel et décroît avec N, avec l'écart du croisement évité.
4. **« Signature d'un fermion de Weyl chiral … avec l'anomalie attendue »** : le flot strict est nul ; seul le flot diabatique vaut −q². C'est une interprétation [C] raisonnable, mais pas un acquis [E]. La compensation (couches 1–3) n'était pas mesurée dans V167. Je la mesure : elle ferme le total une fois convergée.
5. **Mécanisme.** « Pentes très grandes » et « estimateur à fenêtre mal défini » ne correspondent pas aux données. Les causes sont un sous-échantillonnage, et une grille commensurable avec la fenêtre : E_bord ≈ −qθ/L et ε = π/2L placent les sorties de fenêtre sur les points de grille. Ajoutons le point fixe θ = 0, où la paire est hybridée 50/50.
6. **Traçabilité.** `V166_N4_q2_theta0pi.log` n'a ni en-tête ni commande : on ne sait pas quel script l'a produit. Je reproduis ses chiffres (+0,388 fine, +0,775 total). `V166_validation_N3.log` contient une ligne parasite (« 32:def prof… »).
7. **Recommandations.**
   - Grille décalée θ_j = (j + ½)·2π/NTH, avec NTH ≥ 256, ou une fenêtre lisse.
   - Ordre `MMD_AT_PLUS_A` dans le shift-invert : 23–30 s par θ à N = 4, contre ~100 s, et environ 0,5 Go au lieu de ~3 Go.
   - Annoncer le flot strict (nul) à côté du flot diabatique.

## 6. Ce que je n'ai pas pu vérifier
* Le pompage adiabatique de Kubo (M2) à N = 4 : il demande tout le spectre (dimension 18 724), ce qui est hors de portée ici. Je l'ai remplacé par la fenêtre lisse et par ε/2.
* L'intégration exacte suppose E(θ) linéaire entre deux points du balayage, et le poids et la pente interpolés. Pour la zone d'hybridation de q = 1 (pas 10⁻⁴) et de q = 2 (pas 5·10⁻⁴), cela laisse une incertitude que j'estime ≤ 0,01.
* q ≠ 1, 2 ; autres valeurs de m ; transport (a) à N = 4 : non examinés.
* L'opérateur explicite de la symétrie θ → −θ à N ≥ 3 (§3).
* La physique invoquée (anomalie, Callan–Harvey) : lecture seulement.
