# P5 version 2 (solveur eigsh). valider et N = 3 sont deja faits (P5.log, P5_N3.json).
cd "$(dirname "$0")"
date; python3 P5_courbure_oreiller_v2.py mesure 4 96
date; python3 P5_courbure_oreiller_v2.py wilson 5 240
date; echo FIN_P5
