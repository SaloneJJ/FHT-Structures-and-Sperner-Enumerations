#!/usr/bin/env python3
"""D11 (motif A0 : isoler / cone / suspension, groupe Borel) et D12 (appariement neutre : cones, groupe SD16).
Criteres : 08_descente_Xgrand/CRITERES_D11_D12.md (ecrits avant calcul).  Usage : python3 D11_D12.py | tee D11_D12.log  (~ 10-20 min)"""
import sys, os, json, itertools, time
import numpy as np
H_ = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H_, '..', 'AF20_algebres_tordues'))
from fht_tools import aut_order, maximal, betti_z2
from twisted_alg import closure, char_table, mul
S0 = [list(f) for f in json.load(open(os.path.join(H_, '..', 'AF20_algebres_tordues', 'S_spin_min.json')))]
vecs = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
neg = lambda i: vecs.index(((-vecs[i][0]) % 3, (-vecs[i][1]) % 3)); line = lambda i: min(i, neg(i))
T = sorted({tuple(sorted(t)) for t in itertools.combinations(range(8), 3) if all(sum(vecs[i][k] for i in t) % 3 == 0 for k in range(2))})
def mp(M): return tuple(vecs.index(((M[0][0] * a + M[0][1] * b) % 3, (M[1][0] * a + M[1][1] * b) % 3)) for a, b in vecs)
GL = closure([mp([[1, 1], [0, 1]]), mp([[0, 1], [2, 0]]), mp([[2, 0], [0, 1]])])
def vperm(p):
    d = {f'v{i}': f'v{p[i]}' for i in range(8)}
    d.update({f'l{line(i)}': f'l{line(p[i])}' for i in range(8)})
    d.update({f't{n}': f't{T.index(tuple(sorted(p[i] for i in t)))}' for n, t in enumerate(T)})
    return d
def tools(G):
    C_, cl, idx, reps, X = char_table(G); k = len(C_); e0 = [l for l in range(k) if reps[l] == tuple(range(8))][0]
    zp = mp([[2, 0], [0, 2]]); zc = cl[idx[zp]] if zp in idx else None
    irr = [(a, int(round(X[a][e0].real)), (int(round(X[a][zc].real / X[a][e0].real)) if zc is not None else 0)) for a in range(k)]
    chi = lambda a, p: X[a][cl[idx[p]]]
    VP = [vperm(p) for p in G]
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
        R = {d: [rho(g, d) for g in VP] for d in range(-1, maxd + 1)}
        res = {}
        for a, dg, z in irr:
            Q = {}
            for d in range(-1, maxd + 1):
                Pm = sum(np.conj(chi(a, p)) * Rm for p, Rm in zip(G, R[d])) * dg / len(G)
                u, s, vh = np.linalg.svd(Pm); Q[d] = u[:, s > 1e-8]
            for d in range(0, maxd + 1):
                A = D[d] @ Q[d]; ker = Q[d].shape[1] - (np.linalg.matrix_rank(A, 1e-8) if A.size else 0)
                nx = D[d + 1] @ Q[d + 1] if d + 1 <= maxd else np.zeros((1, 0))
                h = ker - (np.linalg.matrix_rank(nx, 1e-8) if nx.size else 0)
                if h: res[(d, a)] = h // dg
        return res
    def orbits(F):
        allv = sorted({x for f in F for x in f}); seen = set(); out = []
        for v in allv:
            if v in seen: continue
            o = {g.get(v, v) for g in VP}; seen |= o; out.append(sorted(o))
        return out
    return irr, eh, orbits, faces_of
def induced(F, W, faces_of):
    return [list(s) for s in maximal([frozenset(s) for s in faces_of(F) if set(s) <= W])]

t0 = time.time()
# ---------------------------------------------------------------- D11
B = [p for p in GL if line(p[0]) == 0]; irr, eh, orbits, faces_of = tools(B)
spin = [a for a, d, z in irr if z == -1]; d2 = [a for a in spin if irr[a][1] == 2][0]; sing = [a for a in spin if irr[a][1] == 1]
hits = []; tested = 0
for bname, F in (('S', S0), ('S+M', S0 + [['l0', 'm']])):
    orbs = orbits(F)
    for r in range(0, 3):
        for I in itertools.combinations(orbs, r):
            Iset = set().union(*map(set, I)) if I else set()
            Fi = [f for f in F if not (set(f) & Iset)] + [[v] for v in sorted(Iset)]
            Vi = {x for f in Fi for x in f}
            cands = [('rien', [], None)]
            for rw in range(1, len(orbs) + 1):
                for W in itertools.combinations(orbs, rw):
                    Wset = set().union(*map(set, W))
                    if not Wset <= Vi: continue
                    Y = induced(Fi, Wset, faces_of)
                    if not Y: continue
                    cands.append(('cone', [o[0] for o in W], [['N'] + y for y in Y]))
                    cands.append(('susp', [o[0] for o in W], [['N'] + y for y in Y] + [['S'] + y for y in Y] + [['N', 'nN']]))
            for op, Wn, add in cands:
                G = Fi + (add or [])
                b = betti_z2(G)
                if b.get(0, 0) < 1 or b.get(1, 0) < 2 or b.get(2, 0) < 1: continue
                tested += 1; res = eh(G)
                sp = {(d, a): m for (d, a), m in res.items() if a in spin}
                if sum(sp.values()) != 3: continue
                if sp.get((1, d2)) != 1: continue
                if sorted(d for (d, a) in sp if a in sing) != [0, 2]: continue
                zp = sum(m * irr[a][1] for (d, a), m in res.items() if a not in spin)
                Vg = sorted({x for g in G for x in g}); hits.append((zp, len(Vg), bname, [o[0] for o in I], op, Wn, res, G))
print(f'D11 : {tested} complexes testes (apres filtre) ; succes : {len(hits)} ({time.time()-t0:.0f} s)', flush=True)
hits.sort(key=lambda h: (h[0], h[1]))
for zp, nv, bname, I, op, Wn, res, G in hits[:8]:
    Vg = sorted({x for g in G for x in g})
    print(f'   z=+1 : {zp} ; sommets {nv} ; base {bname} ; isoles {I} ; {op} sur {Wn} ; |Aut| = {aut_order(Vg, G)} ; homologie {sorted(res.items())}', flush=True)
if hits: json.dump(hits[0][7], open(os.path.join(H_, 'S_A0spin.json'), 'w'))
# ---------------------------------------------------------------- D12
pairing = [{0, 2}, {3, 4}]
def limg(p): return {line(i): line(p[i]) for i in range(8)}
SD = [p for p in GL if sorted(sorted({limg(p)[x] for x in b}) for b in pairing) == sorted(sorted(b) for b in pairing)]
irr2, eh2, orbits2, faces2 = tools(SD); orbs = orbits2(S0); hits2 = []
base = eh2(S0)
for rw in range(1, len(orbs) + 1):
    for W in itertools.combinations(orbs, rw):
        Wset = set().union(*map(set, W)); Y = induced(S0, Wset, faces2)
        if not Y: continue
        G = S0 + [['N'] + y for y in Y]
        if betti_z2(G) != {1: 4}: continue
        if eh2(G) != base: continue
        Vg = sorted({x for g in G for x in g}); au = aut_order(Vg, G)
        if au == 16: hits2.append((len(Wset), [o[0] for o in W], G))
print(f'D12 : homologie de base sous SD16 {sorted(base.items())} ; succes (|Aut| = 16, homologie inchangee) : {len(hits2)} ({time.time()-t0:.0f} s)')
hits2.sort(key=lambda h: h[0])
for n, Wn, G in hits2[:5]: print(f'   cone sur {Wn} ({n} sommets)')
if hits2: json.dump(hits2[0][2], open(os.path.join(H_, 'S_PSneutre.json'), 'w'))
