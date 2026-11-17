#!/usr/bin/env python3
"""D7 -- brisure couleur (marquage du sommet z de B0 : S5 -> S4) ; X_LR et X_SM.  Criteres : 08_descente_Xgrand/CRITERES_D7_brisure_couleur.md.
Usage : python3 D7_brisure_couleur.py | tee D7.log   (< 1 min)"""
import sys, os, json, itertools, time
import numpy as np
H_ = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H_, '..', 'AF20_algebres_tordues'))
from fht_tools import betti_z2, aut_order, maximal
from twisted_alg import closure, char_table, mul, sign, decompose
B = [['z', 'x0', 'x1', 'x2', 'x3'], ['w', 'z'], ['w', 'x0'], ['w', 'x1'], ['w', 'x2'], ['w', 'x3']]
Bm = B + [['z', 'n']]
for name, F in (('B0', B), ('B0 + marque z', Bm)):
    V = sorted({x for f in F for x in f}); print(f'{name} : Betti {betti_z2(F)}, |Aut| {aut_order(V, F)}')
# S4 = stabilisateur de z (indice 0) dans S5 agissant sur {z, x0..x3}
S5 = closure([(1, 0, 2, 3, 4), (1, 2, 3, 4, 0)]); S4 = [p for p in S5 if p[0] == 0]
decompose(S5, sign, 'R^sgn[S5]'); decompose(S4, sign, 'R^sgn[S4] (stab. de z)')
# module H1 = 4 standard (permutation de 5 points moins le trivial) restreint a S4
C, cl, idx, reps, X = char_table(S4); k = len(C); sizes = np.array([len(c) for c in C]); N = len(S4)
chi4 = np.array([sum(1 for i in range(5) if r[i] == i) - 1 for r in reps])
e0 = [l for l in range(k) if reps[l] == tuple(range(5))][0]; eta = np.array([sign(r) for r in reps]); sq = [cl[idx[mul(r, r)]] for r in reps]
print('4|S4 : composantes (degre, type tordu par sgn, multiplicite) :')
for a in range(k):
    m = np.sum(sizes * chi4 * np.conj(X[a])).real / N
    if m < 1e-6: continue
    tw = np.conj(X[a]) * eta; b = [bb for bb in range(k) if np.allclose(X[bb], tw, atol=1e-6)][0]
    typ = 'C (apparie)' if b != a else {1: 'R', -1: 'H'}[int(round(np.sum(sizes * eta * X[a][sq]).real / N))]
    print(f'   degre {int(round(X[a][e0].real))}, type {typ}, multiplicite {m:.0f}')
# assemblages
S0 = json.load(open(os.path.join(H_, '..', 'AF20_algebres_tordues', 'S_spin_min.json')))
ren = lambda f, side: [side + ('w' if x == 'c' else x) for x in f]
G0 = ['g0', 'g1', 'g2', 'g3']; K4 = [list(e) for e in itertools.combinations(G0, 2)]
def assemble(S, Bc):
    X = [ren(a, 'L_') + ren(b, 'R_') + [g] for a in S for b in Bc for g in G0] + [ren(b, 'L_') + ren(a, 'R_') + e for b in Bc for a in S for e in K4]
    X = [list(f) for f in maximal([frozenset(f) for f in X])]; return X, sorted({x for f in X for x in f})
for name, deco, pred in (('X_LR (appariement + couleur)', [['l0', 'l2'], ['l3', 'l4']], (16 * 24) ** 2 * 24), ('X_SM (droite + couleur)', [['l0', 'm']], (12 * 24) ** 2 * 24)):
    X, V = assemble([list(f) for f in S0] + deco, Bm); t = time.time()
    print(f'{name} : sommets {len(V)}, facettes {len(X)}, Betti {betti_z2(X)}, |Aut| = {aut_order(V, X)} (prevu {pred}) ({time.time() - t:.0f} s)', flush=True)
