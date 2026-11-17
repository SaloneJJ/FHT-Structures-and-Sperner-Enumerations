#!/usr/bin/env python3
"""D10 -- conditions sur les J' de X_chi : blocs, J'^2 = +-1, ordre zero, equivalence.  Criteres : CRITERES_D10_conditions_Jprime.md.
Usage : python3 D10_conditions_Jprime.py | tee D10.log   (~ 8 min : refait Aut et l'homologie de AF9)"""
import json, time, os, sys, random
import numpy as np
H_ = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H_, '..', 'AF9_test_triplet'))
import AF9_porteurs_triplet as AF9
from AF21b_aut_raffine import aut_sgs_raffine
P = AF9.P; random.seed(1)
for p in ('X_chi.json', os.path.join('..', '..', '06_objets', 'X_chi.json')):
    if os.path.exists(os.path.join(H_, p)): X = json.load(open(os.path.join(H_, p))); break
t = time.time(); C = AF9.Cx('X_chi', X); gens, order = aut_sgs_raffine(C); nv = len(C.V); V = C.V
sig = [C.ix[AF9.sw(v)] for v in V]; gidx = [C.ix[v] for v in ('g0', 'g1', 'g2', 'g3')]
H = AF9.Homology(C); print('|Aut| =', order, '; generateurs', len(gens), '; Betti', {d: b for d, b in H.betti.items() if b}, f'({time.time()-t:.0f} s)', flush=True)
def sgn_g(g):
    pos = {p: k for k, p in enumerate(gidx)}; perm = [pos[g[p]] for p in gidx]
    return -1 if sum(1 for a in range(4) for b in range(a + 1, 4) if perm[a] > perm[b]) % 2 else 1
pairs = []; Ms = []; Ns = []
for g in gens:
    M, labE = H.block(0, g); h = tuple(sig[g[sig[i]]] for i in range(nv)); N, labO = H.block(1, h)
    pairs.append((M, (sgn_g(g) * N) % P)); Ms.append(M); Ns.append(H.block(1, g)[0])
Bas = AF9.intertwiners(pairs, labE, labO); print('dimC1 =', Bas.shape[0])
# 1. blocs
labE = np.array(labE); labO = np.array(labO)
for tgt in (4, 6):
    sub = [k for k in range(Bas.shape[0]) if np.any(Bas[k][labE == tgt] % P)]
    rows = (Bas[:, labE == tgt, :] % P).reshape(Bas.shape[0], -1); r = AF9.rank_mod(rows.T) if rows.size else 0
    print(f'1. solutions vers le degre {tgt} : dimension {r}')
def inv_mod(A):
    n = A.shape[0]; M = np.concatenate([A % P, np.eye(n, dtype=np.int64)], axis=1)
    for c in range(n):
        r = c + np.nonzero(M[c:, c])[0][0]; M[[c, r]] = M[[r, c]]; M[c] = M[c] * pow(int(M[c, c]), P - 2, P) % P
        o = np.nonzero(M[:, c])[0]; o = o[o != c]; M[o] = (M[o] - np.outer(M[o, c], M[c])) % P
    return M[:, n:]
def rand_A():
    c = np.array([random.randrange(P) for _ in range(Bas.shape[0])], dtype=np.int64); return np.tensordot(c, Bas, axes=1) % P
A1 = rand_A(); A2 = rand_A(); Ai = inv_mod(A1)
assert np.all((A1 @ Ai) % P == np.eye(len(A1), dtype=np.int64))
# 2. J'^2 = +-1 : A' = +-A^-1 verifie A' s(g) rho_pair(sigma g sigma) = rho_impair(g) A'
ok2 = True
for g, (M, Ns_) in zip(gens, pairs):
    h = tuple(sig[g[sig[i]]] for i in range(nv)); Mh, _ = H.block(0, h); No, _ = H.block(1, g)
    for sgn in (1, -1):
        Ap = (sgn * Ai) % P; ok2 &= bool(np.all((Ap @ (sgn_g(g) * Mh) - No @ Ap) % P == 0))
print('2. A\' = +-A^-1 verifie l\'intertwining inverse (J\'^2 = +1 et -1 realisables) :', ok2)
# 3. ordre zero : generateurs gauches (fixent R_ et g)
left = [g for g in gens if all(g[C.ix[v]] == C.ix[v] for v in V if not v.startswith('L_'))]
print('   generateurs purement gauches :', len(left))
def J_rho_Jinv(g):   # sur H = pair (+) impair : J' = [[0, A], [A^-1, 0]] (lineaire sur F_p)
    Me, _ = H.block(0, g); Mo, _ = H.block(1, g)
    # J rho(g) J^-1 : sur pair : A Mo A^-1 ; sur impair : A^-1 Me A
    return (A1 @ Mo @ Ai) % P, (Ai @ Me @ A1) % P
ok3 = True
for g in left:
    for hgen in left:
        Me, _ = H.block(0, g); Mo, _ = H.block(1, g); Je, Jo = J_rho_Jinv(hgen)
        ok3 &= bool(np.all((Me @ Je - Je @ Me) % P == 0)) and bool(np.all((Mo @ Jo - Jo @ Mo) % P == 0))
print('3. ordre zero [rho(g), J\' rho(h) J\'^-1] = 0 (g, h gauches) :', ok3)
# 4. equivalence : U = A2 A1^-1 commute avec rho_pair(Aut) et preserve les degres
U = (A2 @ Ai) % P
ok4 = all(np.all((U @ M - M @ U) % P == 0) for M in Ms)
deg_ok = bool(np.all(U[np.ix_(labE == 4, labE == 6)] % P == 0) and np.all(U[np.ix_(labE == 6, labE == 4)] % P == 0))
print('4. U = A2 A1^-1 dans le commutant de rho_pair(Aut) :', ok4, '; preserve les degres :', deg_ok)
print(f'({time.time()-t:.0f} s)')
