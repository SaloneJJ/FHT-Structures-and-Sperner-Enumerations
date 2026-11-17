#!/usr/bin/env python3
"""D6b -- marquage P (Pati-Salam) et M (type SM) sur les deux cotes de X_grand : Betti et |Aut|.  Criteres : CRITERES_D6_brisure_spin.md (D6b).
Usage : python3 D6b_Xgrand.py P | tee D6b_P.log   (et M) ; duree : quelques minutes (|Aut| par individualisation-raffinement)"""
import sys, os, json, itertools, time
H_ = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H_, '..', 'AF20_algebres_tordues'))
from fht_tools import betti_z2, aut_order, maximal
S = json.load(open(os.path.join(H_, '..', 'AF20_algebres_tordues', 'S_spin_min.json')))
deco = {'P': [['l0', 'l2'], ['l3', 'l4']], 'M': [['l0', 'm']], '0': []}[sys.argv[1]]
S = [list(f) for f in S] + deco
B = [['z', 'x0', 'x1', 'x2', 'x3'], ['w', 'z'], ['w', 'x0'], ['w', 'x1'], ['w', 'x2'], ['w', 'x3']]
ren = lambda f, side: [side + ('w' if x == 'c' else x) for x in f]
G0 = ['g0', 'g1', 'g2', 'g3']; K4 = [list(e) for e in itertools.combinations(G0, 2)]
X = [ren(a, 'L_') + ren(b, 'R_') + [g] for a in S for b in B for g in G0] + [ren(b, 'L_') + ren(a, 'R_') + e for b in B for a in S for e in K4]
X = [list(f) for f in maximal([frozenset(f) for f in X])]; V = sorted({x for f in X for x in f})
t = time.time(); print(f'X_{sys.argv[1]} : sommets {len(V)}, facettes {len(X)}, Betti Z/2 {betti_z2(X)} ({time.time() - t:.0f} s)', flush=True)
t = time.time(); print(f'|Aut| = {aut_order(V, X)} ({time.time() - t:.0f} s)', flush=True)
