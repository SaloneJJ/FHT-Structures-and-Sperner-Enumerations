# Critères P5 — couplage à la courbure du Dirac natif sur la racine fermée (oreiller), comparé à un fermion ordinaire — écrits le 16/11/2026 vers 8 h 15, AVANT calcul
Demande de JJ (16/11, 7 h 51 : point 5, « on va le faire »). Objet : problème ouvert 6 de A3 (facteur 1,46 du terme d'Einstein–Hilbert induit, mesuré sur le tore plat par
une modulation des sauts ; verdict « non résolu », lecture [C]). Lecture directe sur un espace **fermé et courbe**, sans modulation.

## Pourquoi une référence, et laquelle
La racine fermée de A3 (oreiller = double d'un cube, A3-D4 : racine de référence) est **plate par morceaux** : sa courbure est concentrée sur 12 arêtes (déficit π) ; aucune formule
continue simple ne donne son spectre. On compare donc le Dirac natif à un **fermion ordinaire sur exactement la même géométrie** : le Wilson–Dirac à 4 composantes (2 Weyl),
massless, r_W = 1, sur le niveau fin de l'oreiller, avec le même atlas orienté (transitions = demi-tours) et le même principe de signes que D2 (holonomie autour d'une arête =
relèvement continu de la rotation d'angle π). Le Wilson–Dirac a pour limite continue un Dirac standard, couplé minimalement à la courbure.

## Calcul
* Natif : D2 (FHT_28/09_racine_S3/D2_oreiller.py), point R-unif, μ = 1, s = +1, classe de signes de D2_N<N>.json ; K plus petites |E| (inversion spectrale) ; on garde les
  **survivants** du niveau fin (poids fin ≥ 0,5, poids racine < 0,5), comme D2.
* Référence : Wilson–Dirac sur le niveau fin de l'oreiller (2 · 8^N cellules), 2K plus petites |E|.
* Tailles : N = 3 (L = 8) et N = 4 (L = 16), sur le Mac.
* Grandeurs : niveaux distincts (valeurs propres regroupées à 10⁻⁴ près en relatif) et multiplicités ; **rapports sans vitesse** ρ_n = E_n / E_0 (n = 1 … 5) ;
  ajustement E_nat,n² = v² E_W,n² + c ⇒ ξ_eff = c / (v² E_W,0²) (couplage non minimal effectif) ; trace de la chaleur z(t) sur les modes calculés, dans la fenêtre où la troncature
  est négligeable.

## Contrôles (bloquants)
* V1 : Wilson sur le tore (8 cubes périodiques, mode « tore » de D2) = spectre analytique ±√(Σ sin²k + (Σ(1 − cos k))²), à 10⁻¹⁰.
* V2 : hermiticité ; symétrie ± des spectres ; invariance par le signe global des transports de B (jauge).
* V3 : le principe de signes donne, pour les transports de spineurs (σ_a σ_b), la même classe que D2 pour l'octet.

## Lecture fixée d'avance
* **Standard** si l'écart max |ρ_n(natif) − ρ_n(Wilson)| / ρ_n(Wilson) sur n = 1 … 5 est ≤ 3 % à N = 4 et ne croît pas de N = 3 à N = 4 : le survivant se couple à la courbure
  comme un fermion de Weyl ordinaire ; son coefficient d'Einstein–Hilbert induit à basse énergie est le coefficient standard, et **le facteur 1,46 du tore n'est pas confirmé**
  sur la racine fermée (il viendrait du protocole de modulation ou des modes à l'échelle de coupure).
* **Anomal** si l'écart est ≥ 10 % aux deux tailles, de même signe : le couplage effectif ξ_eff est donné ; comparaison qualitative au 1,46.
* Sinon : **non résolu** (le dire).
* Multiplicités : le natif (s = +1) doit avoir la moitié des multiplicités du Wilson (une Weyl contre deux) ; un autre rapport est rapporté tel quel.
* Limite déclarée : le coefficient d'Einstein–Hilbert de l'action spectrale intègre tous les modes jusqu'à la coupure ; ce test ne porte que sur l'espèce légère.

## Amendement déclaré (16/11, ≈ 8 h 35) — après l'essai de validation à N = 2 et la lecture des spectres déjà publiés de D2, AVANT tout calcul P5 à N ≥ 3
* Sur l'oreiller, le survivant natif n'a de niveaux **sous la bande de volume** qu'à N = 4 : D2 (N = 4) donne 0,2156 (×4), 0,3301 / 0,3303 (×8), 0,3929 (×4) sous
  E_bulk(tore) = 0,774 ; à N = 3 aucun niveau (S³ trop petit, rayon effectif ≈ 3,7). La comparaison porte donc sur **N = 4 seul**, sur les niveaux natifs sous **0,75 · E_bulk**
  (E_bulk lu dans le JSON de D2), au moins trois niveaux distincts.
* La condition « ne croît pas de N = 3 à N = 4 » est remplacée par un contrôle de la **référence** : Wilson seul à L = 32 (N = 5, niveau fin), pour mesurer ses propres
  artefacts de réseau ; les rapports ρ_n de Wilson doivent varier de moins de 3 % entre L = 16 et L = 32, sinon la comparaison est non résolue.
* Lecture : **standard** si max |écart| ≤ 5 % sur les niveaux disponibles ; **anomal** si ≥ 10 % ; sinon non résolu. ξ_eff est rapporté dans tous les cas.
* Le cas N = 3 est calculé quand même (quelques minutes) et rapporté, sans lecture.

## Amendement 2 (16/11, ~21 h 25 UTC) — procédure d'appariement réécrite. DÉCLARÉ : écrit APRÈS la lecture littérale de N = 4 et sa relecture post hoc (RESULTATS_P5.md)
Ce n'est donc pas un critère aveugle pour N = 4. Il le sera pour le contrôle de Wilson à L = 32 et pour tout calcul ultérieur.
- Niveaux : valeurs |E| triées, fusionnées en grappes quand l'écart relatif entre voisins est < 2×10⁻³ (au lieu de 10⁻⁴). Niveau = moyenne pondérée de la grappe ; multiplicité = taille de la grappe.
- Fenêtre : |E| < 0.75 E_bulk = 0.5806, la même que la lecture littérale.
- Appariement : niveaux natifs et de Wilson pris dans l'ordre. **Condition de validité** : chaque paire a un rapport de multiplicités natif/Wilson égal à 1/2 (une Weyl contre deux).
  Si une paire échoue, l'appariement s'arrête à la paire précédente ; s'il reste moins de trois paires valides, le verdict est « non résolu ».
- Grandeurs : ρ_n = E_n/E_0 sur les paires valides ; écarts relatifs ; ajustement E_nat² = v² E_W² + c par moindres carrés sur les paires valides ; ξ_eff = c/(v² E_W,0²).
- Seuils inchangés : standard si max |écart| ≤ 5 %, anomal si ≥ 10 %, non résolu sinon. Le contrôle de la référence (Wilson L = 16 contre L = 32, variation des ρ_n < 3 % avec la même procédure) reste requis pour un verdict définitif.
- Une sensibilité au seuil de fusion (10⁻³ et 5×10⁻³) est rapportée.

## Amendement 3 (16/11, ~21 h 45 UTC, avant tout calcul de ce type) : méthode du contrôle Wilson L = 32
- La factorisation LU de `wilson 5` (dimension 262 144) a saturé la mémoire du Mac (19 Go d'échange) et a été arrêtée sans résultat.
- Nouvelle méthode, sans factorisation : itération de sous-espace filtrée par polynômes de Tchebychev sur A = H², blocs de 160 vecteurs, 120 valeurs propres voulues, résidu ‖A x − λ x‖ < 10⁻⁷. Mémoire estimée ≈ 3 Go à N = 5. Script : P5_wilson_iter.py, qui importe build_W et classe_D2 de la v2 sans les modifier.
- **Validation bloquante (essai à N = 4)** : les 112 plus petites |E| de Wilson, dans la fenêtre |E| < 0.5806, doivent coïncider à 10⁻⁶ près avec celles de P5_N4.json (calculées par LU et eigsh). Sinon, on ne lance pas N = 5.
- Contrôle (amendement 1, procédure de l'amendement 2) : on compare les 10 premiers niveaux de Wilson à L = 32 à ceux de L = 16 (grappes à 2×10⁻³, multiplicités égales requises). Le contrôle passe si les ρ_n varient de moins de 3 %.
- Verdict final de P5 :
  - contrôle passé : la lecture « standard » de N = 4 est confirmée ;
  - contrôle échoué : « non résolu » ;
  - on rapporte aussi l'écart natif–Wilson réestimé en retranchant la dérive de Wilson entre L = 16 et L = 32.
