#!/usr/bin/env python3
"""D8b -- assemblage X_chi avec la piece spinorielle a degres separes (suspension de la classe de parallelisme marquee) et la couleur marquee.
Criteres : CRITERES_D8_separation_degres.md (D8b).  Usage : python3 D8b_assemblage.py | tee D8b.log"""
import sys, os, json, itertools, time
H_ = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H_, '..', 'AF20_algebres_tordues'))
from fht_tools import betti_z2, aut_order, maximal
S0 = [list(f) for f in json.load(open(os.path.join(H_, '..', 'AF20_algebres_tordues', 'S_spin_min.json')))]
vecs = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
T = sorted({tuple(sorted(t)) for t in itertools.combinations(range(8), 3) if all(sum(vecs[i][k] for i in t) % 3 == 0 for k in range(2))})
ta = 't%d' % T.index(tuple(sorted(vecs.index((1, j)) for j in range(3)))); tb = 't%d' % T.index(tuple(sorted(vecs.index((2, j)) for j in range(3))))
W = {'l2', 'l3', 'l4', ta, tb} | {f'v{i}' for i in range(2, 8)}
faces = set()
for f in S0:
    for r in range(1, len(f) + 1):
        for s in itertools.combinations(sorted(f), r): faces.add(s)
Y = [list(s) for s in maximal([frozenset(s) for s in faces if set(s) <= W])]
Schi = S0 + [['N'] + y for y in Y] + [['S'] + y for y in Y] + [['N', 'nN']]
print('droites affines paralleles a l0 :', ta, tb, '; S_chi : Betti', betti_z2(Schi), '|Aut|', aut_order(sorted({x for f in Schi for x in f}), Schi))
json.dump(Schi, open(os.path.join(H_, 'S_chi.json'), 'w'))
Bm = [['z', 'x0', 'x1', 'x2', 'x3'], ['w', 'z'], ['w', 'x0'], ['w', 'x1'], ['w', 'x2'], ['w', 'x3'], ['z', 'n']]
ren = lambda f, side: [side + ('w' if x == 'c' else x) for x in f]
G0 = ['g0', 'g1', 'g2', 'g3']; K4 = [list(e) for e in itertools.combinations(G0, 2)]
X = [ren(a, 'L_') + ren(b, 'R_') + [g] for a in Schi for b in Bm for g in G0] + [ren(b, 'L_') + ren(a, 'R_') + e for b in Bm for a in Schi for e in K4]
X = [list(f) for f in maximal([frozenset(f) for f in X])]; V = sorted({x for f in X for x in f})
t = time.time(); print(f'X_chi : sommets {len(V)}, facettes {len(X)}, Betti {betti_z2(X)} ({time.time() - t:.0f} s)', flush=True)
t = time.time(); print(f'|Aut(X_chi)| = {aut_order(V, X)} (prevu 1990656) ({time.time() - t:.0f} s)', flush=True)
json.dump(X, open(os.path.join(H_, 'X_chi.json'), 'w'))
