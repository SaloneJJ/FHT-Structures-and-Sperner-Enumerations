# AUDIT 14/11/2026 -- reproduction du noyau de la chaleur a N = 3 (A3, Prop. « Heat kernel »), grille m = 6, constructeur de V352 (V352_arbre.Hnat).
# Usage (depuis CalculCodes) : python3 AUDIT_FABLE/AF13_noyau_N3_calcul.py 0 32   (ecrit n3/*.npy dans le dossier courant ; ~ 14 s par point, 32 points)
import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import sys, os, itertools, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'V352')); import V352_arbre as M
m, N = 6, 3; a, b = int(sys.argv[1]), int(sys.argv[2])
def canon(j): return min([tuple(sorted(j)), tuple(sorted((-x) % m for x in j))])
orb = {}
for j in itertools.product(range(m), repeat=3): orb[canon(j)] = orb.get(canon(j), 0) + 1
reps = sorted(orb)
for r in reps[a:b]:
    f = 'n3/%d%d%d.npy' % r
    if os.path.exists(f): continue
    t = time.time(); E = np.linalg.eigvalsh(M.Hnat(N, 2 * np.pi * np.array(r) / m, +1, mu=1.0)); np.save(f, E)
    print(r, orb[r], round(time.time() - t, 1), flush=True)
