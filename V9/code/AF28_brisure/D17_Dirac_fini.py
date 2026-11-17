#!/usr/bin/env python3
"""D17 -- operateurs de Dirac finis admissibles (gamma-impair, J'D = DJ', ordre un) le long de la descente.
Criteres : 08_descente_Xgrand/CRITERES_D17_Dirac_fini.md (ecrits avant calcul).  Usage : python3 D17_Dirac_fini.py [noms] | tee D17.log"""
import sys, os, json, itertools, time, random
import numpy as np
H_ = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.join(H_, '..', 'AF9_test_triplet'), os.path.join(H_, '..', 'AF20_algebres_tordues')): sys.path.insert(0, p)
import AF9_porteurs_triplet as AF9
from fht_tools import maximal
P = AF9.P; random.seed(7)
from collections import defaultdict, Counter

def aut_sgs_base(C, order_list):
    """comme AF21b, mais base dans l'ordre donne ; renvoie (generateurs, ordre, niveau de chaque generateur)"""
    n = len(C.V); F = C.F; m = len(F); N = n + m
    adj = [[] for _ in range(N)]
    for j, f in enumerate(F):
        for v in f: adj[v].append(n + j); adj[n + j].append(v)
    Fset = C.Fset
    def refine_pair(c1, c2):
        while True:
            s1 = [(c1[x], tuple(sorted(c1[y] for y in adj[x]))) for x in range(N)]
            s2 = [(c2[x], tuple(sorted(c2[y] for y in adj[x]))) for x in range(N)]
            keys = sorted(set(s1) | set(s2)); k = {s: i for i, s in enumerate(keys)}
            n1 = [k[s] for s in s1]; n2 = [k[s] for s in s2]
            if len(set(n1)) == len(set(c1)) and len(set(n2)) == len(set(c2)): return n1, n2
            c1, c2 = n1, n2
    def search(c1, c2):
        c1, c2 = refine_pair(c1, c2)
        if Counter(c1) != Counter(c2): return None
        cells = defaultdict(list); cells2 = defaultdict(list)
        for x in range(n): cells[c1[x]].append(x); cells2[c2[x]].append(x)
        nontriv = [c for c in cells if len(cells[c]) > 1]
        if not nontriv:
            mp = {cells[c][0]: cells2[c][0] for c in cells}; perm = [mp[v] for v in range(n)]
            return perm if all(frozenset(perm[v] for v in f) in Fset for f in F) else None
        c = min(nontriv, key=lambda c: (len(cells[c]), c)); v = cells[c][0]; big = max(max(c1), max(c2)) + 1
        for w in cells2[c]:
            d1 = list(c1); d2 = list(c2); d1[v] = big; d2[w] = big
            r = search(d1, d2)
            if r is not None: return r
        return None
    def orbit(pt, gens):
        seen = {pt}; fr = [pt]
        while fr:
            nf = []
            for x in fr:
                for g in gens:
                    y = g[x]
                    if y not in seen: seen.add(y); nf.append(y)
            fr = nf
        return seen
    col, _ = refine_pair([0] * n + [1] * m, [0] * n + [1] * m)
    gens_all = []; levels = []; order = 1
    for li, v in enumerate(order_list):
        cand = [w for w in range(n) if col[w] == col[v]]
        if len(cand) == 1: continue
        big = max(col) + 1; lev = []; orb = {v}
        for w in cand:
            if w in orb: continue
            d1 = list(col); d2 = list(col); d1[v] = big; d2[w] = big
            r = search(d1, d2)
            if r is not None: lev.append(tuple(r)); orb = orbit(v, lev)
        order *= len(orb); gens_all.extend(lev); levels.extend([li] * len(lev))
        d1 = list(col); d1[v] = big; col, _ = refine_pair(d1, d1)
    return gens_all, order, levels

def nullcomb(L):
    """combinaisons c (lignes) avec c L = 0 mod P"""
    return AF9.nullspace_mod(L.T)

def analyse(name, X):
    t = time.time(); C = AF9.Cx(name, X); nv = len(C.V); V = C.V
    nonL = [C.ix[v] for v in V if not v.startswith('L_')]; Lv = [C.ix[v] for v in V if v.startswith('L_')]
    gens, order, levels = aut_sgs_base(C, nonL + Lv)
    left = [g for g in gens if all(g[u] == u for u in nonL)]
    sig = [C.ix[AF9.sw(v)] for v in V]; gidx = [C.ix[v] for v in ('g0', 'g1', 'g2', 'g3')]
    H = AF9.Homology(C); betti = {d: b for d, b in H.betti.items() if b}
    print(f'== {name} : |Aut| = {order} ; generateurs {len(gens)} dont {len(left)} purement gauches ; Betti {betti} ({time.time()-t:.0f} s)', flush=True)
    def sgn_g(g):
        pos = {p: k for k, p in enumerate(gidx)}; perm = [pos[g[p]] for p in gidx]
        return -1 if sum(1 for a in range(4) for b in range(a + 1, 4) if perm[a] > perm[b]) % 2 else 1
    pairs = []
    for g in gens:
        M, labE = H.block(0, g); h = tuple(sig[g[sig[i]]] for i in range(nv)); N, labO = H.block(1, h)
        pairs.append((M, (sgn_g(g) * N) % P))
    Bas = AF9.intertwiners(pairs, labE, labO)
    nE, nO = len(labE), len(labO)
    if Bas.shape[0] == 0 or nE != nO: print('   pas de J\' (dimC1 = %d)' % Bas.shape[0]); return
    for _ in range(5):
        A = np.tensordot(np.array([random.randrange(P) for _ in range(Bas.shape[0])], dtype=np.int64), Bas, axes=1) % P
        if AF9.rank_mod(A) == nE: break
    # inverse mod P
    M_ = np.concatenate([A % P, np.eye(nE, dtype=np.int64)], axis=1)
    for c in range(nE):
        r = c + np.nonzero(M_[c:, c])[0][0]; M_[[c, r]] = M_[[r, c]]; M_[c] = M_[c] * pow(int(M_[c, c]), P - 2, P) % P
        o = np.nonzero(M_[:, c])[0]; o = o[o != c]; M_[o] = (M_[o] - np.outer(M_[o, c], M_[c])) % P
    Ai = M_[:, nE:]
    Me = {g: H.block(0, g)[0] for g in left}; Mo = {g: H.block(1, g)[0] for g in left}
    Se = {g: (A @ Mo[g] % P) @ Ai % P for g in left}; So = {g: (Ai @ Me[g] % P) @ A % P for g in left}   # J rho(g) J^-1 sur pair / impair
    oz = all(np.all((Me[g] @ Se[h] - Se[h] @ Me[g]) % P == 0) and np.all((Mo[g] @ So[h] - So[h] @ Mo[g]) % P == 0) for g in left for h in left)
    # inconnue Q (impair <- pair : nO x nE), P_ = A Q A (pair <- impair)
    k = nO * nE; Bq = np.eye(k, dtype=np.int64)
    for g in left:
        for h in left:
            if Bq.shape[0] == 0: break
            Qs = Bq.reshape(-1, nO, nE)
            Ps = np.einsum('ij,kjl->kil', A, Qs) % P; Ps = np.einsum('kij,jl->kil', Ps, A) % P
            T1 = (np.einsum('kij,jl->kil', Ps, Mo[g]) - np.einsum('ij,kjl->kil', Me[g], Ps)) % P      # bloc pair <- impair de [D, rho(g)]
            T2 = (np.einsum('kij,jl->kil', Qs, Me[g]) - np.einsum('ij,kjl->kil', Mo[g], Qs)) % P      # bloc impair <- pair
            C1 = (np.einsum('kij,jl->kil', T1, So[h]) - np.einsum('ij,kjl->kil', Se[h], T1)) % P
            C2 = (np.einsum('kij,jl->kil', T2, Se[h]) - np.einsum('ij,kjl->kil', So[h], T2)) % P
            L = np.concatenate([C1.reshape(Bq.shape[0], -1), C2.reshape(Bq.shape[0], -1)], axis=1) % P
            cmb = nullcomb(L); Bq = (cmb @ Bq) % P if cmb.shape[0] else np.zeros((0, k), dtype=np.int64)
    print(f'   ordre zero : {oz} ; dim (sur F_p) des D admissibles (gamma-impair, J\'D = DJ\', ordre un) = {Bq.shape[0]}  ({time.time()-t:.0f} s)', flush=True)
    # repartition par degres (bloc Q : impair <- pair)
    if Bq.shape[0]:
        labE = np.array(labE); labO = np.array(labO); out = {}
        for de in sorted(set(labE)):
            for do in sorted(set(labO)):
                sub = Bq.reshape(-1, nO, nE)[:, labO == do, :][:, :, labE == de].reshape(Bq.shape[0], -1)
                r = AF9.rank_mod(sub.T) if sub.size else 0
                if r: out[f'{de}->{do}'] = r
        print('   rang par bloc de degres (pair -> impair) :', out, flush=True)

d = json.load(open(os.path.join(H_, 'V25_K3_pieces.json'))); G4, E4 = d['generation']['feuilles'], d['generation']['K4']
def buildK(A, B):
    I = [frozenset(['L_' + v for v in a] + ['R_' + v for v in b] + [g]) for a in A for b in B for g in G4]
    II = [frozenset(['L_' + v for v in b] + ['R_' + v for v in a] + e) for b in B for a in A for e in E4]
    return [sorted(f) for f in maximal(I + II)]
S0 = [list(f) for f in json.load(open(os.path.join(H_, 'S_spin_min.json')))]
B0 = [['z', 'x0', 'x1', 'x2', 'x3'], ['w', 'z'], ['w', 'x0'], ['w', 'x1'], ['w', 'x2'], ['w', 'x3']]
ren = lambda f, side: [side + ('w' if x == 'c' else x) for x in f]
K4 = [list(e) for e in itertools.combinations(G4, 2)]
def buildG(S, Bc):
    X = [ren(a, 'L_') + ren(b, 'R_') + [g] for a in S for b in Bc for g in G4] + [ren(b, 'L_') + ren(a, 'R_') + e for b in Bc for a in S for e in K4]
    return [sorted(f) for f in maximal([frozenset(f) for f in X])]
loadj = lambda f: json.load(open(os.path.join(H_, f)))
cx = {'K26': lambda: buildK(d['A'], d['B']), 'X_comp': lambda: loadj('X_compresse.json'), 'X_chi': lambda: loadj('X_chi.json'),
      'X_SM': lambda: buildG(S0 + [['l0', 'm']], B0 + [['z', 'n']]), "X_LR'": lambda: buildG(loadj('S_PS_neutre.json'), B0 + [['z', 'n']]),
      "X_PS'": lambda: buildG(loadj('S_PS_neutre.json'), B0), 'X_grand': lambda: buildG(S0, B0)}
for name in (sys.argv[1:] or list(cx)):
    analyse(name, cx[name]())
