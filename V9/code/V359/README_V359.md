# V359 — scripts de relecture indépendante de V347–V350, V353, V354
Écrits par le relecteur, sans importer le code relu. Données : V37_FHT_K26volmin.json, V36_FHT_K26vol.json, V32_FHT_K26etoile.json (02_objets).
* V359_combinatoire.py (+ .log) : Betti (GF(32003)) des zones induites et de K26 ; caractère de S4_L sur H1(L-X) ; cône de R_z ; 3-cycles et aires ;
  secteurs I/II ; Fiedler ; transpositions impaires de S4_L, S4_R ; liens des coins ; homologie des sous-complexes à 1 g et à 2 g.  `python3 V359_combinatoire.py <dossier JSON>`
* V359_aires.py (+ .log) : aires des triangles des règles R1–R3 de V350 (V37, V36).
* V359_torsion.py (+ .log) : torsion par SVD (V37, V36), contrôles C3, simplexe, RP² ; superdéterminant ; Betti mod 2, 3, 5, 7. (V32 omis : coût.)
* V359_parites.py (+ .log) : parités (Γ4, FΓ5, S) sous les 48 relevés de O_h (code propre) ; F = −Γ4 S.
* V359_phi.py (+ .log) : dépendance en φ (énergie du vide et Tr e^{−tH²}) avec un Clifford différent (conjugaison unitaire aléatoire), grille 10³ non décalée,
  μ ∈ {−1, 0, 0,5, 1, 1,5, 2,5}, signe de Wilson, masse m0.
numpy seul (calculs faits dans l'espace de travail cloud, recopiés).
