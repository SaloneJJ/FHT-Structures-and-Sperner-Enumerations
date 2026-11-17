# AUDIT_FABLE — scripts de l'audit des trois articles (Fable 5.2, 13/11/2026)

Tous ces scripts sont **écrits de zéro** par l'auditeur : aucun n'importe ni ne recopie un script de CalculCodes. Ils ne lisent que l'énoncé des articles
(et, pour AF3, les objets de `04_objets`). Chacun affiche l'énoncé contrôlé et un verdict ; le journal de l'exécution du 13/11 est à côté (`.log`).
Environnement : Python 3 et numpy seulement (pas de scipy). Aucun n'écrit de fichier. Lancement : depuis ce dossier, `python3 AFn_….py`.

| script | article, registre | ce qu'il vérifie | durée | verdict du 13/11 |
|---|---|---|---|---|
| `AF1_enumeration_independante.py` | A1 ; F1–F2 | énumération des FHT à partir des seuls axiomes : force brute sur toutes les familles de parties (v ≤ 4), puis décomposition récursive (v ≤ 5) : 1, 1, 4, 23, 231 ; nombres étiquetés 8, 149, 8 488 | 4 s | conforme |
| `AF2_algebre_tordue_S4.py` | A1 ; F6 | ℝ^sgn[S₄] ≅ ℂ ⊕ ℍ ⊕ M₃(ℂ) par les constantes de structure (centre, idempotents, signature de la forme trace) ; les 15 groupes d'ordre 24 : seul S₄ a cinq classes ; seul (S₄, sgn) donne A_F | 5 s | conforme |
| `AF3_K26_homologie_automorphismes.py` | A2 ; T3–T4 | K26 reconstruit d'après la Définition : 716 arêtes, 47 167 simplexes, Betti (12, 36, 36, 12), \|Aut\| = 13 824 ; pièces A et B ; mêmes contrôles sur V37 (542), V36 (750), V32 (490) | 5 s | conforme |
| `AF4_clifford_octet.py` | A3 ; C6–C9, M1–M3 | 35 identités de l'octet (Lemme « Clifford identities », commutant, BBH, J_lat, O_φ, angle φ) | 1 s | conforme ; relève « F* = −F » (faux : F est hermitien) |
| `AF5_chaine_reduite_pont.py` | A3 ; C14, M1, M2, M4 | chaîne de Kaplan réduite : deux branches (goût et λ), no-go Γ₅ ⊗ D_F, masse par O₀ (uniforme : ε_D ; couche fine : ε_D(1 − λ²)), pont m(h) | 1 min | conforme ; relève l'énoncé (iii) de la Prop. « Insertion » |
| `AF6_flux_pi_energie.py` | A3 ; C4 (V205, jamais relu) | énergie à demi-remplissage : cube (2¹² configurations), amas 4³ et 6³, réseau périodique 8³ | 20 s | conforme |
| `AF7_action_spectrale_controles.py` | A3 ; C3, A3–A9 | niveau plat (8,04 ; 3,04), traces de jauge, anomalies, table des couplages à une boucle, facteur du temps | 1 s | conforme ; relève l'en-tête des deux dernières colonnes de la table |
| `AF8_bimodule_controles.py` | A2 ; T2, T5 | unicité de T = AB^t + BA^t, dimension 32, restrictions à D₄, 84 paramètres, forme d'intersection, H₁(K₄) = 3′ | 2 s | conforme |

Programmes de l'auteur relancés par l'audit (sur une copie, résultats dans les rapports 05_audit/AUDIT_FABLE_A*.md) : `FHT_Enumeration.py`, `FHT_Invariants.py`,
`autonomes/t3_saveur.py`, `t4_tm1_cp.py`, `t5_tm1_phases.py`, `V355/V355_controle.py`, `V356/V356_chaine_P0.py`, `V336_FHT_flot.py` (N = 3, R-unif),
`V342_FHT_masse_survivant.py`, `V343_FHT_flot_insertion.py 0.3 AY 64`, `V345_pont.py`, `V352/V352_noyau.py` (N = 2).
Pour lancer V335–V345 il faut `PYTHONPATH=V331_relecture` (module `rvA_common`) et scipy.
