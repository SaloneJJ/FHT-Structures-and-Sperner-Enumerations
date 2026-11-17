# Relecture indépendante de V159, V161 et V163 (piste B, R2–R3)

Date : 2026-11-12. Relecteur : Claude (agent de relecture). Je n'ai pas écrit les calculs relus.
Étiquettes : **[E]** recalculé par moi, avec mon propre code ; **[T]** vérifié par lecture du code ou par preuve ; **[C]** lecture ou jugement.
Scripts et journaux : `/home/claude/fht/relecture2/` (aucun fichier existant n'a été modifié).
* `rv_core.py` : géométrie, liens de jauge et modèles, réécrits sans reprendre le code de V158. Les liens sont construits par une règle générique : intégrale exacte d'un potentiel linéaire, plus fonction de transition de la section. Deux jauges sont disponibles : G1 = (0, BX, θ/L), celle de V158, et G2 = (−BY, 0, θ/L).
* `rv_test_gauge.py` : holonomies, lignes de Wilson, comparaison exacte avec les matrices de V158.
* `rv_flow.py` : flot par couche, avec trois estimateurs (voir ci-dessous).
* `rv_spectres.py` : unicité du singulet, et spectres avec et sans flux.
* `rv_lemme.py` : lemme de V158 §1.
* `rv_jauge.py` : transports (a) et (b) réécrits ; covariance ; facteurs de forme ; robustesse ; B-P4 ; flot.

Les trois estimateurs du flot dans `rv_flow.py` :
* **M1** : fenêtre de Hellmann–Feynman, comme V158, avec ε et ε/2 ;
* **M2** : pompage adiabatique de Kubo sur tout le spectre, F(n) = −∮ Σᵢ w′ᵢ(n) G(Eᵢ). Il est égal à M1 dans la limite continue en θ (intégration par parties), mais discrétisé autrement, et son total est nul par construction ;
* **M3** : suivi des niveaux par recouvrement et comptage des passages par E = 0.

Contrainte d'exécution : un autre processus du projet (`V162_FHT_R3_jauge.py`, deux instances) a occupé jusqu'à 1,5 cœur pendant toute la relecture. Tous mes calculs ont tourné sur un seul fil.

---

## Tableau des verdicts

| point | verdict |
|---|---|
| 1a. Construction du flux (A, translation magnétique, liens verticaux, H(2π) ≅ H(0)) | **confirmé [E]** |
| 1b. Formule du flot (fenêtre, Hellmann–Feynman, poids par couche) comme mesure de chiralité | **confirmé avec réserve [T/E]** : correct dans le plat ; dans l'arbre à N = 3, la répartition par couche dépend de l'estimateur à ~10 % |
| 1c. −0,99 (plat) et −0,93 (arbre) | **confirmé [E]** (reproduits avec mon code, dans une autre jauge) ; −0,93 n'est pas robuste (−0,84 à −0,94 selon l'estimateur) |
| 2. Lemme de V158 §1 | **confirmé, avec une erreur d'énoncé [T/E]** : l'opérateur qui commute est Ω = Γ₁Γ₂Γ₃Γ₄Γ₅, et non Γ₆ (Γ₆ anticommute avec Γ₁…Γ₅) |
| 3a. Implémentation de Kawamoto–Smit (η, η₄ sur l'enfant, masse ε) | **confirmé [T/E]** |
| 3b. Bord gappé (0,670 à M = 1 ; 0,105 à M = 0,3) | **confirmé avec réserve [E]** : chiffres reproduits sous flux ; sans flux, à M = 0,3, un état de bord est à **0,009 (N = 3) et 0,003 (N = 4)** |
| 3c. Lecture « vectoriel, puis gappé » | **confirmé avec réserve [C]** : solide à M = 1 ; « gappé » est trop fort à M = 0,3 |
| 4a. Unicité de a = (1, 0, 1, 0) | **confirmé [E]** (L = 4 et L = 6), dans l'ansatz testé |
| 4b. Flot ∓1,98 dans le plat (bords ouverts) | **confirmé [E]** : ∓1,983 (ε) et ∓1,994 (ε/2). Multiplicité 2 expliquée (§4) |
| 4c. « Volume de l'arbre sans gap » | **confirmé [E]**, et aussi sans flux. L'argument « opérateur affaibli » (élément max 0,042 contre 0,062) est un artefact du flux |
| 4d. « Même obstacle que V119–V120 » | **lecture plausible, non démontrée [C]** |
| 5a. Transport (b) : covariance, chemins, moyenne | **confirmé [E/T]** |
| 5b. Règle de la racine | **vérifié qu'elle est nécessaire [E]** (moyenne exactement nulle) ; prescription arbitraire [C] |
| 5c. f_n = 0 exactement pour j ≥ n | **confirmé avec réserve [E]** : vrai pour j < N et pour un fond uniforme seulement ; faux pour j = N ; seulement approché (5–17 %) pour un fond non uniforme |
| 5d. Test central 0,17 contre 3,19 | **reproduit [E], mais mal interprété** : ‖Vψ‖ n'est pas la mesure annoncée dans V162 ; avec la mesure annoncée, la conclusion s'inverse ; « 2 états par θ » est faux (2 états en tout) |
| 5e. Classes n̄ ≥ N − 0,5 et ≤ N − 1,5 | **fragile [E]** : à 16 θ, le rapport passe de 0,17 à 0,44 si l'on reclasse les 2 états intermédiaires |
| 5f. B-P1 « confirmée » | **trop fort [E]** : à 48 pas, le total vaut −0,15 et la compensation (+0,78) sort de la tolérance ±0,05 de V162 ; à 96 pas, −0,936 / +0,871 |
| 5g. Rapprochement avec Grabowska–Kaplan | **globalement honnête, avec réserves [C]** ; références non vérifiables (pas de réseau) |

---

## 1. Diagnostic de flot spectral (V158 §2, V159 §1.1)

### 1a. Construction du flux : confirmée [E]
* **Comparaison avec V158.** Mes matrices en jauge G1 coïncident exactement avec celles de V158 (|ΔH| ≤ 3·10⁻¹⁷ ; dérivée en θ, 7·10⁻¹⁸), pour l'arbre et le plat, et pour les modèles `bor` et `ksO`. La transition ψ(X + L) = e^{iBLY} ψ(X) a le bon signe.
* **Holonomies.** J'ai tiré environ 300 boucles fermées aléatoires par cas (marche aléatoire, puis retour par le plus court chemin). Elles comprennent des liens verticaux, des enroulements et les niveaux grossiers s = 1 et 2, que V158 ne contrôle pas (`plaquettes()` saute s < 3). Pour toutes les boucles contractiles, l'holonomie vaut e^{iB·aire projetée} à 2·10⁻¹⁵ près, dans les jauges G1 et G2, pour l'arbre et le plat.
* **Lignes de Wilson.** Les lignes de Wilson non contractiles sont identiques en G1 et en G2 (6·10⁻¹⁶).
* **Spectres.** Le spectre est le même en G1 et en G2 (2·10⁻¹³). H(2π) et H(0) ont le même spectre (10⁻¹³). L'équivalence est d'ailleurs évidente : A_z = θ/L est une jauge pure, λ = θZ/L, sauf pour l'holonomie en z [T].
* `rv_test_gauge_N3.log`.

### 1b. Ce que mesure la formule [T]
* Le flot fenêtré est F(n) = ∮ Σ_{|E|<ε} wᵢ(n) Eᵢ′ / 2ε. Par intégration par parties, il est égal à −∮ Σᵢ wᵢ′(n) G(Eᵢ), où G est une marche lissée de largeur ε : c'est la **charge pompée adiabatiquement par la mer de Dirac dans la couche n** (afflux de Callan–Harvey). C'est une mesure légitime de la chiralité (anomalie du niveau de Landau zéro). Le flot attendu vaut q² Φ par tour.
* **Réserve importante.** Dans l'arbre à N = 3, **aucun niveau ne traverse E = 0**. Le niveau chiral du bord et l'état grossier s'évitent avec un écart de 0,0063 (je le retrouve, ainsi que V159).
  - Le flot spectral strict (comptage des passages par zéro, M3) est donc **identiquement nul couche par couche**.
  - Le −0,93 est un flot « diabatique » qui n'existe que parce que ε ≫ 0,006. Ce n'est pas un défaut, car le pompage de caractère est physique. Mais le chiffre n'est pas un invariant topologique entier, et il dépend de ε.
* **Bases propres.** Le produit (pente diagonale) × (poids diagonal) dépend de la base choisie dans un sous-espace dégénéré. Ce n'est pas un problème ici : je n'ai pas vu de dégénérescence dans la fenêtre, et M2, qui ne dépend pas de la base, donne les mêmes chiffres.

### 1c. Chiffres reproduits [E]
| cas | M1, ε (pas) | M1, ε/2 | M2, Kubo (pas) | V159 |
|---|---|---|---|---|
| plat, Wilson m = −1 | −0,990 / +0,990 (48) | −0,995 | −0,955 (48) | −0,990 |
| arbre, Wilson m = −1, couche fine | −0,929 (64) ; −0,936 (128) | −0,770 (64) ; −0,842 (128) | −0,843 (32) | −0,929 |
| arbre, couches 1 et 2 | +0,514 / +0,414 (64) | +0,523 / +0,404 (128) | +0,467 / +0,377 | +0,514 / +0,414 |

* Le signe et l'ordre de grandeur sont confirmés, ainsi que l'absence de compensation à la racine (0,000) et la compensation dans les couches 1–2.
* Dans l'arbre à N = 3, la valeur de la couche fine se situe entre −0,84 et −0,94 selon l'estimateur. « −0,93 » est un chiffre, pas une mesure convergée. V159 ne le présente pas comme entier, mais l'écrit sans barre d'erreur.
* Le calcul à N = 4 (`rv_flow_wil_tree_N4.log`) n'a pas abouti : il a dépassé 25 min sur un seul fil, avec le processeur partagé (§7).
* M3 (comptage) échoue dans le plat : les deux bords traversent zéro au même θ, sur le même bloc k_z, et s'hybrident. M3 n'est donc pas utilisable comme contrôle dans ce cas ; je ne l'utilise pas comme preuve.

## 2. Lemme de V158 §1 [T, contrôle E dans `rv_lemme.log`]
* Les 6 générateurs de Jordan–Wigner (ΠZ)X_j et (ΠZ)Y_j anticommutent et engendrent M₈(ℂ) : exact.
* **Erreur d'énoncé.** « Γ₆ (∝ Γ₁Γ₂Γ₃Γ₄Γ₅) commute avec tout l'hamiltonien » est faux tel quel. Γ₆ **anticommute** avec Γ₁…Γ₅ (‖[Γ₆, Γⱼ]‖ = 2) et n'est pas proportionnel à Γ₁⋯Γ₅ (écart minimal √2).
  - L'opérateur qui commute est **Ω = Γ₁Γ₂Γ₃Γ₄Γ₅**. Il est central dans Cl(5), vérifie Ω² = 1 et a des valeurs propres ±1 de multiplicité 4 ; on a Ω = ±i ω Γ₆, avec ω la chiralité de Cl(6).
* **La conclusion reste juste.** Tout hamiltonien de la forme Σ (opérateur de réseau) ⊗ (polynôme en Γ₁…Γ₅) commute avec Ω. Sur le mode de paroi (−iΓ₄Γ₅ = +1), Γ₁Γ₂Γ₃ vaut −i dans le secteur Ω = +1 et +i dans le secteur Ω = −1 : les hélicités sont opposées, et le bord est vectoriel.
* **Portée [T].** Le lemme ne couvre que les hamiltoniens de Clifford sur la fibre ℂ⁸.
  - Avec Γ₆ en plus (d·Γ, d ∈ ℝ⁶), C₂ = 0 encore, car H⁴(S⁵) = 0.
  - Un hamiltonien ℂ⁸ local quelconque, non linéaire en Clifford, n'est pas couvert. « Fibre locale ⇒ vectoriel » n'est donc pas un théorème général. Il faudrait écrire « fibre locale ⇒ vectoriel, pour les hamiltoniens de Dirac–Wilson à coefficients dans Cl(6) ».

## 3. V159 : Kawamoto–Smit sur l'arbre

### 3a. Implémentation : confirmée [T/E]
* Phases η₁ = 1, η₂ = (−1)^{x₁}, η₃ = (−1)^{x₁+x₂}.
* η₄ = (−1)^{x₁+x₂+x₃} est évalué sur l'enfant (ligne de S). Dans l'arbre, il vaut donc (−1)^{Σc} : le parent se couple à la **combinaison alternée** des 8 enfants. C'est un choix fort, que V159 signale comme « autres phases verticales » non testées.
* Masse M·ε, avec ε = (−1)^{Σx+n}. Dans le plat, ε anticommute avec tous les sauts ; dans l'arbre, non.
* Mes matrices sont identiques à celles de V158.

### 3b. Bord gappé : reproduit, avec une réserve [E] (`rv_spectres_BC.log`, θ = 0)
| cas | Φ = 1 (V159) | Φ = 1 (relecture) | **Φ = 0 (relecture)** |
|---|---|---|---|
| arbre N = 4, M = 1, O ou W | 0,670 | 0,6704 | bord 0,624 ; global 0,527 (n̄ 1,7) |
| arbre N = 4, M = 0,3, O | 0,105 | 0,1048 | **0,0032 (n̄ 3,5)** |
| arbre N = 3, M = 0,3, O | — | 0,147 | **0,0090 (n̄ 2,5)** |
| plat N = 3, M = 1, O / W | 1,047 / 0,651 | 1,047 / 0,651 | 1,047 / 0,651 |

* À M = 1, le bord est bien gappé, avec ou sans flux.
* À M = 0,3, le « gap » de 0,105 est mesuré **sous flux, à θ = 0**. Sans flux, il existe un état proche du bord dont l'énergie décroît avec N (0,009 puis 0,003). C'est un seul état, et non une branche. Il est donc sans conséquence pour la chiralité : le flot net est nul. Mais « bord gappé à M = 0,3 » est trop fort. Un gap mesuré dans un fond de flux n'est pas le gap à champ nul.
* Le flot de l'arbre N = 3 (`ksO`), −0,132 « à la racine », est un artefact de discrétisation de M1. M2 donne un total nul et +0,010 à la racine, et ε/2 donne −0,002.
* (O) et (W) sont identiques dans l'arbre. Explication probable [C] : la masse ε, pondérée par la profondeur n, change de signe avec n, et le saut vertical alterné η₄ rend la paroi presque invisible. Non approfondi.

### 3c. Lecture [C]
« Le bord est vectoriel (flot net 0), puis gappé » est cohérent avec la classification. Dans le plat, la masse ε = γ₅ ⊗ ξ₅ donne C₂ = 2·(½) − 2·(½) = 0, donc aucun mode de bord protégé ; la paroi ne porte pas non plus de chiralité nette. S-P0 (« 4 goûts sans masse à la paroi ») était naïve, et V159 la déclare correctement réfutée.

## 4. V161 : masse singulet de goût

### 4a. Unicité : confirmée [E] (`rv_spectres_A.log`)
* Avec un code indépendant, sur les tores L = 4 et L = 6, un seul couple (a, phase) parmi 16 × 2 vérifie les trois conditions : a = (1, 0, 1, 0), phase 1.
* Cela correspond à la phase η₅ = (−1)^{x₁+x₃} de l'opérateur singulet γ₅ ⊗ 1 à 4 liens de la littérature des fermions décalés [C, de mémoire, non vérifié en ligne].
* Réserve : l'unicité ne vaut que dans l'ansatz f·C₁C₂C₃C₄, avec f = ±1 et une phase dans {1, i}.

### 4b. Flot ∓1,98 dans le plat : confirmé [E]
* M1 donne ∓1,983 (ε) et ∓1,994 (ε/2) à N = 3, avec 48 pas.
* Sans flux, on trouve 4 modes exactement nuls par bord, ce qui correspond à 2 fermions de Weyl à deux composantes en k = 0.
* **La multiplicité 2, que V161 dit « non expliquée », s'explique [T, argument standard] :**
  - chaque goût massif contribue ±½ à C₂ ;
  - avec le singulet, les 4 goûts ont le même signe, d'où C₂ = ±2 ;
  - un bord ouvert est une interface avec le vide (C₂ = 0) et porte |ΔC₂| = 2 Weyl ;
  - une paroi +M | −M porterait |ΔC₂| = 4, ce qui est la prédiction C-P1 de V160 ;
  - le modèle de Wilson de V151 a C₂ = 1, d'où 1 Weyl par bord.
* « Pas de matrice » : correct. La chiralité vient de la structure du modèle de Kawamoto–Smit. Le mécanisme est celui, connu, des fermions décalés à paroi ou à bord [C].

### 4c. « Volume de l'arbre sans gap » : confirmé [E]
* Arbre N = 3 et 4, θ = 0 : les spectres de V161 sont reproduits à 10⁻³ près.
* **Sans flux**, le constat est plus net encore :
  - à N = 4, 23 états ont |E| < 0,2, dont 13 avec n̄ ≤ 2,5 ;
  - il y a 6 modes exactement nuls sur la couche fine.
* Dans le plat ouvert sans flux, après les modes de bord, le premier état est à 0,708 : le volume plat est gappé, celui de l'arbre ne l'est pas.
* Flot de l'arbre N = 3, couche fine : −0,432 (O) et +0,109 (W) reproduits exactement. Les autres estimateurs donnent −0,55 et −0,46 (O), et +0,30 et +0,12 (W). La valeur n'est pas entière et dépend de l'estimateur, ce qui confirme qu'il n'y a pas de branche propre.
* **Réserve sur l'argument « affaibli » (norme max 0,042 contre 0,062).**
  - C'est le plus grand **élément de matrice**, mesuré **sous flux**. Sans flux, l'élément max de l'arbre vaut 0,177, plus que le plat (0,0625). Cet indicateur n'a pas de sens.
  - En norme d'opérateur, l'arbre est bien plus faible : 0,47 contre 0,81 sans flux, 0,27 contre 0,77 sous flux. L'idée survit, mais avec un autre chiffre.
* La valeur « ‖{Γ₅, K}‖ = 0,062 dans le plat » est une anticommutation brisée par le flux et les bords, pas par l'arbre. Je ne l'ai pas séparée davantage [non vérifié].

### 4d. « Même obstacle que V119–V120 » [C]
C'est plausible : la non-commutation des décalages verticaux et latéraux sur l'arbre empêche le saut diagonal d'anticommuter avec la cinétique. Mais rien n'est démontré : une seule famille (C₄ = moyenne des enfants, η₄ alterné, moyenne symétrisée sur 24 ordres) a été testée. La formulation de V161 (« semblent en tension… pas démontré en général ») est honnête.

## 5. V163 : jauge U(1), transport par blocs

### 5a. Implémentation de (b) : confirmée [E/T]
* Réécrit indépendamment (`rv_jauge.py`) : jauge du parent = enfant 0 ; liens verticaux = chemin x, puis y, puis z dans le bloc ; lien latéral = moyenne des 8 chemins 0 → c → c + 2eⱼ → 0′, puis projection sur U(1).
* Mon transport (a) coïncide exactement avec la géométrie G2 (|ΔT| = |ΔS| = 0).
* **Covariance** sous une jauge aléatoire du niveau fin : 4·10⁻¹⁴ (q = 1) et 1·10⁻¹³ (q = 2), avec a = 0,2.
* Flux de 1,0000 aux niveaux 1 à 3. Normes des moyennes brutes : 0,9952 (niveau 2) et 0,9239 (niveau 1) ; **0,0000 en x et y à la racine** (1,0 en z). Les chiffres de V163 sont reproduits.

### 5b. Règle de la racine [E/C]
* La moyenne est exactement nulle, pour la raison donnée : les holonomies en x aux hauteurs L/4 et 3L/4 valent ∓i.
* Mettre le lien à 0 supprime aussi le terme de Wilson latéral (x, y) de la racine. La masse effective de la racine passe alors de m + 4r − r(cos x + cos y + cos z) à m + 4r − r cos z. C'est donc un changement **du modèle** à la racine, et pas seulement du transport. Il faut le dire, et tester d'autres prescriptions : r1 (moyenne non projetée), ou un lien égal à −1 comme dans (a).

### 5c. f_n = 0 exactement pour j ≥ n [E] (`rv_jauge_base_N3.log`, `rv_jauge_ffgen_N3.log`)
* **Confirmé pour j < N et un fond uniforme**, en RMS par cellule comme en norme d'opérateur (la définition de V162 est la norme ; le code calcule une RMS) : résidus ≤ 6·10⁻¹⁰. Cela vaut aussi pour des perturbations A_x(Y) et A_z(X).
* **Faux pour j = N** (ℓ = maille fine) : f = [0 ; 4 ; 2 ; 1]. Le cosinus est déjà replié en une constante au niveau fin. Il faut écrire « j < N ».
* **Mécanisme exact [T].** Ce n'est pas seulement « un nombre entier de périodes par cellule ».
  - Au premier ordre, la phase projetée vaut Σ_c w_c δ_c, avec w_c = Re(e^{iφ_c}/Σ e^{iφ}).
  - Avec l'ordre de chemin x, y, z, les phases φ_c ne prennent que **deux valeurs**, selon c_x. Pour deux phases, l'identité Re(1/(1 + e^{iβ})) = Re(e^{iβ}/(1 + e^{iβ})) = ½ donne des poids égaux.
  - Le lien grossier est alors **exactement la somme** des liens fins de la colonne, et cette somme s'annule sur une période complète.
* **Avec un fond non uniforme** (bruit de A_x(X, Y), amplitude 0,3), l'annulation n'est plus exacte :
  - f_{0,1}/f_N = 0,07 / 0,05 (yX, j = 1) ;
  - 0,11 à 0,17 (xY, j = 1 et 2) ;
  - 0,04 (zX, j = 1).
* « Coupure nette » est donc une propriété du fond uniforme. En général, la suppression est forte mais approchée. C'est plus proche de Grabowska–Kaplan que ne le dit V163, mais dans le sens d'une suppression incomplète.

### 5d. Test central B-P4 [E] (`rv_jauge_p4_N3_16.log`, `_64.log`)
* **Reproduction.** Avec ‖Vψ‖ et 16 θ, je retrouve exactement les chiffres de V163 : 0,7069 / 0,1204 / 0,170 pour (b), j = 2, et 0,7126 / 2,2711 / 3,19 pour (a).
* **Effectifs.** « 2 états de compensation par θ » (V163 §1.3) est **faux** : il y a **2 états en tout** sur les 16 θ, à θ = 0,39 et 5,89 (les 8 états de bord sont aussi un total). Avec 64 θ, on trouve 8 à 10 états, et le rapport reste à 0,172.
* **Mesure.** V162 §2 annonçait « la somme des |⟨ψ′|V|ψ⟩|² sur les états ψ′ de la fenêtre ». Le code calcule ‖Vψ‖, qui somme sur **tous** les états, y compris ceux de haute énergie. L'écart n'est pas signalé dans V163.
* **Avec la mesure annoncée, la conclusion s'inverse :**

| j | transport | bord (médiane) | miroir (médiane) | rapport miroir/bord |
|---|---|---|---|---|
| 1 | (a) | 1,5·10⁻⁷ | 7,3·10⁻⁵ | ≈ 490 |
| 1 | (b) | 1,1·10⁻⁷ | 5,1·10⁻⁵ | ≈ 460 |
| 2 | (a) | 4,7·10⁻¹⁰ | 4,2·10⁻⁷ | ≈ 880 |
| 2 | (b) | 1,7·10⁻¹⁰ | 3,1·10⁻⁷ | ≈ 1 800 |

  - (64 θ.) Dans la fenêtre, le miroir est **plus** couplé que le bord, et (b) ne réduit son couplage de basse énergie que d'un facteur ~1,4 par rapport à (a).
  - Les nombres absolus sont minuscules, parce que la fenêtre ne contient qu'un ou deux états. Cette mesure est donc peu informative elle aussi. Mais elle ne soutient pas B-P4.
* **Ce que mesure ‖Vψ‖ [T].**
  - C'est la norme du sommet sur le support de ψ. Elle se calcule à partir du facteur de forme (B-P3) et du profil de profondeur de ψ (miroir : [0 ; 0,53 ; 0,44 ; 0,03]).
  - Le rapport 0,17 est donc une **reformulation de B-P3**, pondérée par le profil, et non un test indépendant du découplage.
  - Le 3,19 de (a) vient surtout de la longueur des liens grossiers : f_n ∝ 2^{N−n}, quelle que soit la longueur d'onde. Ce n'est pas une absence de suppression au sens physique : un A_y uniforme donnerait le même facteur.
  - Un indicateur pertinent serait l'effet de V sur les observables de basse énergie : déplacement d'énergie au second ordre des états miroirs, ou réponse en courant ou en charge à la longueur d'onde k. Il n'a pas été calculé [non vérifié].

### 5e. Classes [E]
* Les seuils n̄ ≥ 2,5 et ≤ 1,5 (N = 3) excluent les états intermédiaires de θ = 0 (n̄ ≈ 2,2, profil [0,01 ; 0,27 ; 0,21 ; 0,51]).
* Si on les compte comme miroir (n̄ < 2,5), le rapport de (b), j = 2, passe de **0,17 à 0,44** à 16 θ. Il reste à 0,17 à 64 θ, car les états intermédiaires y sont minoritaires.
* Avec 16 θ et 2 états miroir, le chiffre central repose sur **deux états** (une paire ±E à deux θ). Ce n'est pas assez pour en tirer « 6 fois ».

### 5f. B-P1 et B-P2 [E]
* **B-P1** (transport b, q = 1, N = 3) :
  - 96 pas, fenêtre ε : couche fine −0,936 ; couches grossières +0,871 ; total −0,065 ;
  - 96 pas, fenêtre ε/2 : couche fine −0,823 ;
  - 48 pas (V163) : −0,924 / +0,779, total −0,145.
  - Le total non nul de V163 est bien un défaut de discrétisation. Mais la compensation de V163 (+0,779) est à 0,15 de V159 (+0,929), hors de la tolérance ±0,05 fixée par V162. « Confirmée » devrait être « confirmée qualitativement ; non convergée à 48 pas ».
* **B-P2** (q = 2) : V163 la déclare non concluante, ce qui est honnête. Je ne l'ai pas recalculée [non vérifié].
  - Remarque [T] : le flux de 1,0000 affiché au niveau 1 pour q = 2 dans (a) est une ambiguïté de mesure (plaquette de phase π).
  - Remarque [C] : avec q = 2 à N = 3, la longueur magnétique vaut environ 2,3 mailles, et les états de Landau ne sont pas résolus.

### 5g. Grabowska–Kaplan [C] (référence non vérifiable ici, pas de réseau)
* Ce qui est honnête :
  - l'analogie annoncée comme hypothèse [C] ;
  - les différences reconnues : coupure nette contre suppression douce, et perte de topologie à la racine.
* Ce qui manque :
  1. Chez Grabowska et Kaplan, le miroir est une **seconde paroi** portant un champ fermionique complet et propageant. Ici, le « miroir » se réduit à quelques états grossiers (2 par cycle à N = 3). Le problème résolu n'est pas le même : la piste B n'a déjà presque pas de miroir (V152), et la question est surtout l'**anomalie** (où va la charge pompée), traitée par B-P1.
  2. Le transport par blocs est une transformation de type RG par blocs (moyenne de chemins et projection, type APE). L'analogie avec un flot de gradient est générique : ce sont deux lissages.
  3. L'argument central de Grabowska–Kaplan (un champ fermionique miroir qui ne voit qu'un champ lissé, avec la question de la localité et des vertex de 't Hooft) n'a ici aucun équivalent testé.
  4. La conservation de la topologie par le flot de Grabowska–Kaplan est affirmée de mémoire. La référence arXiv:1911.10925 citée dans V162 n'a pas pu être vérifiée.

## 6. Autres remarques
* V163 B-P3 : V162 définit f_n par une norme d'opérateur ; le code utilise une RMS par cellule. Les conclusions ne changent pas (je l'ai vérifié), mais l'écart n'est pas signalé.
* V159 §1.2, « −0,13 (racine) » à N = 3 : artefact du pas en θ (§3b).
* Les lignes « min|E| bord » de `V158_KS.log` ne viennent d'aucun script versionné. Je les ai reproduites, mais le script qui les a produites manque.
* V161 C-P2 est jugée confirmée sur la base de N = 3. À N = 4, seul le spectre a été obtenu, pas le flot. C'est honnêtement signalé.

## 7. Ce que je n'ai pas pu vérifier
* B-P2 (q = 2) et B-P4 à N = 4 : non recalculés, faute de temps de calcul (CPU partagé).
* Le flot de l'arbre de Wilson à N = 4 : **non obtenu** (arrêté à 25 min, journal vide). On ne sait donc pas si la dépendance en ε de la couche fine (−0,84 à −0,94 à N = 3) diminue avec N. C'est le contrôle prioritaire à faire, avec une version qui ne factorise qu'une fois par θ.
* Les références (Grabowska–Kaplan PRL 116, 211602 ; arXiv:1911.10925 ; phase η₅ = (−1)^{x₁+x₃}) : de mémoire seulement.
* Le comptage strict des passages par zéro (M3) n'est pas fiable dans le plat (§1c).
