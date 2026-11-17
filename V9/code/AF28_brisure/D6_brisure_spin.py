#!/usr/bin/env python3
"""D6 -- premiere brisure sur le complexe : S_spin_min decoree (P : appariement des droites -> Pati-Salam ; M : droite marquee -> type SM).
Criteres : 08_descente_Xgrand/CRITERES_D6_brisure_spin.md (ecrits avant calcul).  Usage : python3 D6_brisure_spin.py | tee D6.log  (< 1 min)"""
import sys, os, json, itertools
import numpy as np
H_ = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H_, '..', 'AF20_algebres_tordues'))
from fht_tools import betti_z2, aut_order, maximal
from twisted_alg import closure, char_table, mul
S = json.load(open(os.path.join(H_, '..', 'AF20_algebres_tordues', 'S_spin_min.json')))
vecs = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
neg = lambda i: vecs.index(((-vecs[i][0]) % 3, (-vecs[i][1]) % 3))
line = lambda i: min(i, neg(i))
def mp(M): return tuple(vecs.index(((M[0][0] * a + M[0][1] * b) % 3, (M[1][0] * a + M[1][1] * b) % 3)) for a, b in vecs)
GL = closure([mp([[1, 1], [0, 1]]), mp([[0, 1], [2, 0]]), mp([[2, 0], [0, 1]])]); assert len(GL) == 48
def det(p):
    c1 = vecs[p[vecs.index((1, 0))]]; c2 = vecs[p[vecs.index((0, 1))]]
    return 1 if (c1[0] * c2[1] - c2[0] * c1[1]) % 3 == 1 else -1
def lines_img(p): return {line(i): line(p[i]) for i in range(8)}
def chi_odd(p): return (sum(1 for i in range(8) if p[i] == i) - sum(1 for i in range(8) if p[i] == neg(i))) / 2
zc_perm = mp([[2, 0], [0, 2]])
from collections import Counter
def order(p):
    q = p; k = 1
    while q != tuple(range(8)): q = mul(q, p); k += 1
    return k

def analyse(name, Hgrp, F):
    V = sorted({x for f in F for x in f})
    print(f'\n== {name} : sommets {len(V)}, facettes {len(F)}, Betti Z/2 {betti_z2(F)}, |Aut| (aut_order) = {aut_order(V, F)}, '
          f'stabilisateur dans GL(2,3) : {len(Hgrp)}, ordres {dict(sorted(Counter(order(p) for p in Hgrp).items()))}')
    C, cl, idx, reps, X = char_table(Hgrp); k = len(C); N = len(Hgrp); sizes = np.array([len(c) for c in C])
    e0 = [l for l in range(k) if reps[l] == tuple(range(8))][0]; zc = cl[idx[zc_perm]]
    eta = np.array([det(r) for r in reps]); sq = [cl[idx[mul(r, r)]] for r in reps]
    chi4 = np.array([chi_odd(r) for r in reps])
    mult = [np.sum(sizes * chi4 * np.conj(X[j])).real / N for j in range(k)]
    out = []; seen = set()
    for a in range(k):
        d = X[a][e0].real; z = X[a][zc].real / d
        tw = np.conj(X[a]) * eta; b = [bb for bb in range(k) if np.allclose(X[bb], tw, atol=1e-6)][0]
        if b != a: typ = 'C (apparie)'
        else:
            nu = np.sum(sizes * eta * X[a][sq]).real / N; typ = {1: 'R', -1: 'H'}.get(int(round(nu)), f'? {nu:.2f}')
        out.append((a, int(round(d)), int(round(z)), typ, b, round(mult[a], 6)))
    print('   irreductibles (degre, z, type pour la structure tordue par det, multiplicite dans 4_spin) :')
    for a, d, z, typ, b, m in out:
        if z == -1: print(f'     [z=-1] degre {d} type {typ:12s} multiplicite {m}' + (f' (partenaire {b})' if typ.startswith('C') else ''))
    fac = []; used = set()
    for a, d, z, typ, b, m in out:
        if z != -1 or a in used: continue
        used |= {a, b}
        fac.append({'R': f'M{d}(R)', 'H': f'M{d // 2}(H)' if d > 2 else 'H'}.get(typ[0], f'M{d}(C)'))
    print('   secteur spinoriel (z = -1) de R^det[Aut] :', ' + '.join(fac))
    comp = [(d, typ) for a, d, z, typ, b, m in out if m > 1e-6]
    print('   4_spin|Aut =', ' + '.join(f'{d}({t[0]})' for d, t in comp))

Fbase = [list(f) for f in S]
analyse('S_spin_min (controle)', GL, Fbase)
pairing = [{0, 2}, {3, 4}]
HP = [p for p in GL if sorted(sorted({lines_img(p)[x] for x in b}) for b in pairing) == sorted(sorted(b) for b in pairing)]
analyse('P : appariement {l0,l2}|{l3,l4}', HP, Fbase + [['l0', 'l2'], ['l3', 'l4']])
HM = [p for p in GL if lines_img(p)[0] == 0]
analyse('M : droite l0 marquee', HM, Fbase + [['l0', 'm']])
