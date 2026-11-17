# P5 : couplage a la courbure sur l'oreiller (criteres : CRITERES_P5.md). Tous les coeurs. Duree estimee : 30 min a 1 h 30.
cd "$(dirname "$0")"
date; python3 P5_courbure_oreiller.py valider
date; python3 P5_courbure_oreiller.py mesure 3 96
date; python3 P5_courbure_oreiller.py mesure 4 96
date; python3 P5_courbure_oreiller.py wilson 5 240
date; echo FIN_P5
