#!/usr/bin/env python3
"""D16' -- compression : piece isospin recollee au simplexe couleur avec un S4 diagonal ; recherche du motif A0.
Criteres : CRITERES_D14_D16.md (D16', declare apres l'echec de D15).  Usage : python3 D16p_compression.py | tee D16p.log"""
import sys, os, json, itertools, time
import numpy as np
H_ = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(H_, 'D11_D12.py')).read().split("t0 = time.time()")[0]; exec(src)
t0 = time.time()
B = [p for p in GL if line(p[0]) == 0]
AGL = [p for p in B if vecs[p[vecs.index((1, 0))]][0] == 1]
ta = 't%d' % T.index(tuple(sorted(vecs.index((1, j)) for j in range(3))))
linesb = ['l2', 'l3', 'l4']                     # directions (1,0), (1,1), (1,2) : b = 0, 1, 2
PAIRS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
pkey = lambda P: frozenset(frozenset(b) for b in P)
def spair(s):  # permutation induite sur les appariements
    return [[pkey(P) for P in PAIRS].index(frozenset(frozenset(s[x] for x in b) for b in P)) for P in PAIRS]
def lperm(a): return [linesb.index(f'l{line(a[vecs.index(v)])}') for v in [(1, 0), (1, 1), (1, 2)]]
S4 = list(itertools.permutations(range(4)))
diag = []
for s in S4:
    a = [a for a in AGL if lperm(a) == spair(s)]; assert len(a) == 1; diag.append((s, a[0]))
def vd(s, a):
    d = vperm(a); d.update({f'x{i}': f'x{s[i]}' for i in range(4)}); return d
VPd = [vd(s, a) for s, a in diag]; Glab = [s for s, a in diag]
C_, cl, idx, reps, X = char_table(Glab); k = len(C_); e0 = [l for l in range(k) if reps[l] == (0, 1, 2, 3)][0]
nm = {}
for a_ in range(k):
    d = int(round(X[a_][e0].real))
    nm[a_] = {1: ('1' if all(np.isclose(X[a_][l], 1) for l in range(k)) else 'sgn'), 2: 'E', 3: '3'}[d] + ('' if d != 3 or X[a_][cl[idx[(1, 0, 2, 3)]]].real > 0 else "'")
chi = lambda a_, s: X[a_][cl[idx[s]]]
def faces_of(F):
    out = set()
    for f in F:
        f = tuple(sorted(f))
        for r in range(1, len(f) + 1):
            for s in itertools.combinations(f, r): out.add(s)
    return out
def eh(F):
    fs = faces_of(F); byd = {}
    for s in fs: byd.setdefault(len(s) - 1, []).append(s)
    maxd = max(byd); byd[-1] = [()]
    ind = {d: {s: i for i, s in enumerate(sorted(byd[d]))} for d in byd}
    D = {}
    for d in range(0, maxd + 1):
        M = np.zeros((len(ind[d - 1]), len(ind[d])))
        for s, j in ind[d].items():
            for t in range(len(s)): M[ind[d - 1][s[:t] + s[t + 1:]], j] += (-1) ** t
        D[d] = M
    def rho(g, d):
        n = len(ind[d]); M = np.zeros((n, n))
        for s, j in ind[d].items():
            im = [g.get(x, x) for x in s]; perm = sorted(range(len(im)), key=lambda i: im[i]); sg = 1; q = perm[:]
            for i in range(len(q)):
                while q[i] != i: j2 = q[i]; q[i], q[j2] = q[j2], q[i]; sg = -sg
            M[ind[d][tuple(sorted(im))], j] = sg
        return M
    R = {d: [rho(g, d) for g in VPd] for d in range(-1, maxd + 1)}
    res = {}
    for a_ in range(k):
        dg = int(round(X[a_][e0].real)); Q = {}
        for d in range(-1, maxd + 1):
            Pm = sum(np.conj(chi(a_, s)) * Rm for s, Rm in zip(Glab, R[d])) * dg / len(Glab)
            u, sv, vh = np.linalg.svd(Pm); Q[d] = u[:, sv > 1e-8]
        for d in range(0, maxd + 1):
            A = D[d] @ Q[d]; ker = Q[d].shape[1] - (np.linalg.matrix_rank(A, 1e-8) if A.size else 0)
            nx = D[d + 1] @ Q[d + 1] if d + 1 <= maxd else np.zeros((1, 0))
            h = ker - (np.linalg.matrix_rank(nx, 1e-8) if nx.size else 0)
            if h: res[(d, nm[a_])] = h // dg
    return res
link = []
for kk, P in enumerate(PAIRS):
    for b in P: link.append([linesb[kk], 'z', f'x{b[0]}', f'x{b[1]}'])
C0 = S0 + [[ta, 'o']] + [['z', 'x0', 'x1', 'x2', 'x3']] + link
V0 = sorted({x for f in C0 for x in f})
print(f"C0 : sommets {len(V0)}, Betti {betti_z2(C0)}, |Aut| = {aut_order(V0, C0)} ; homologie sous S4 diagonal : {sorted(eh(C0).items())}", flush=True)
seen = set(); orbs = []
for v in V0:
    if v in seen: continue
    o = {g.get(v, v) for g in VPd}; seen |= o; orbs.append(sorted(o))
print('orbites :', [o[0] + f'({len(o)})' for o in orbs], flush=True)
tg0 = {(0, '1'): 1, (1, 'E'): 1, (2, 'sgn'): 1}; tg1 = {(1, '1'): 1, (2, 'E'): 1, (3, 'sgn'): 1}
hits = []; tested = 0
for r in range(0, 3):
    for I in itertools.combinations(orbs, r):
        Iset = set().union(*map(set, I)) if I else set()
        Fi = [f for f in C0 if not (set(f) & Iset)] + [[v] for v in sorted(Iset)]
        Vi = {x for f in Fi for x in f}
        cands = [('rien', [], [])]
        for rw in range(1, len(orbs) + 1):
            for W in itertools.combinations(orbs, rw):
                Wset = set().union(*map(set, W))
                if not Wset <= Vi: continue
                fsi = faces_of(Fi); Y = [list(s) for s in maximal([frozenset(s) for s in fsi if set(s) <= Wset])]
                if not Y: continue
                cands.append(('cone', [o[0] for o in W], [['N'] + y for y in Y]))
                cands.append(('susp', [o[0] for o in W], [['N'] + y for y in Y] + [['S'] + y for y in Y] + [['N', 'nN']]))
        for op, Wn, add in cands:
            G = Fi + add; b = betti_z2(G)
            if sum(b.values()) != 4: continue
            tested += 1; res = eh(G)
            for tn, tg in (('type 0', tg0), ('type 1', tg1)):
                if res == tg: hits.append((len({x for g in G for x in g}), tn, [o[0] for o in I], op, Wn, G))
print(f'D16\' : {tested} testes ; succes : {len(hits)} ({time.time() - t0:.0f} s)', flush=True)
hits.sort(key=lambda h: h[0])
for nv, tn, I, op, Wn, G in hits[:8]:
    print(f'   {tn} : {nv} sommets ; isoles {I} ; {op} sur {Wn} ; |Aut| = {aut_order(sorted({x for g in G for x in g}), G)}', flush=True)
if hits: json.dump(hits[0][5], open(os.path.join(H_, 'A_compresse.json'), 'w'))
