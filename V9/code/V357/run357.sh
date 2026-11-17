# V357 : flot sur l'arbre (V336_FHT_flot.py tel quel, mode unif, t_f = ε), à lancer depuis 05_code/CalculCodes (V336 lit V162 dans le dossier courant)
# ordre réel : ε = 0,2 et 0,1 (14 h 50), puis 0,6 et 0,0 (15 h 00 – 15 h 20), heure AST
for e in 0.2 0.1 0.6 0.0; do python3 V336_FHT_flot.py 3 0 1 1 64 1 unif $e > V357/log357_eps$e.txt 2>&1; done
# spectres de Bloch N = 2 pour V357_noyau_analyse.py (dans V357/v352run) :
#   python3 ../V352/V352_noyau.py 2 8 1 1.0 ; python3 ../V352/V352_noyau.py 2 8 1 1.0 0.1 AY ; python3 ../V352/V352_noyau.py 2 8 1 1.0 0.1 G5
