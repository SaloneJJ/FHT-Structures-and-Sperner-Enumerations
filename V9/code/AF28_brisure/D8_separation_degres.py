#!/usr/bin/env python3
"""D8 -- separation des degres (chiralite) : cones et suspensions sur les sous-complexes induits B-invariants de la piece spinorielle marquee.
Criteres : 08_descente_Xgrand/CRITERES_D8_separation_degres.md (ecrits avant calcul).  Usage : python3 D8_separation_degres.py | tee D8.log"""
import sys, os, json, itertools, time
import numpy as np
H_ = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H_, '..', 'AF20_algebres_tordues'))
from fht_tools import aut_order, maximal
from twisted_alg import closure, char_table, mul
S0 = [list(f) for f in json.load(open(os.path.join(H_, '..', 'AF20_algebres_tordues', 'S_spin_min.json')))]
vecs = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
neg = lambda i: vecs.index(((-vecs[i][0]) % 3, (-vecs[i][1]) % 3)); line = lambda i: min(i, neg(i))
T = sorted({tuple(sorted(t)) for t in itertools.combinations(range(8), 3) if all(sum(vecs[i][k] for i in t) % 3 == 0 for k in range(2))})
def mp(M): return tuple(vecs.index(((M[0][0] * a + M[0][1] * b) % 3, (M[1][0] * a + M[1][1] * b) % 3)) for a, b in vecs)
GL = closure([mp([[1, 1], [0, 1]]), mp([[0, 1], [2, 0]]), mp([[2, 0], [0, 1]])])
B = [p for p in GL if line(p[0]) == 0]; assert len(B) == 12
def vperm(p):
    d = {f'v{i}': f'v{p[i]}' for i in range(8)}
    d.update({f'l{line(i)}': f'l{line(p[i])}' for i in range(8)})
    d.update({f't{n}': f't{T.index(tuple(sorted(p[i] for i in t)))}' for n, t in enumerate(T)})
    return d
C_, cl, idx, reps, X = char_table(B); k = len(C_); e0 = [l for l in range(k) if reps[l] == tuple(range(8))][0]
zc = cl[idx[mp([[2, 0], [0, 2]])]]
irr = [(a, int(round(X[a][e0].real)), int(round(X[a][zc].real / X[a][e0].real))) for a in range(k)]
chi_of = lambda a, p: X[a][cl[idx[p]]]
# V0 : caracteres de f_l0 et f_t (t_a = droite affine x = 1, t_b = x = 2)
ta = [vecs.index((1, j)) for j in range(3)]
def act_f(f, p):   # (g.f)(v) = f(g^-1 v)
    inv = [0] * 8
    for i, j in enumerate(p): inv[j] = i
    return np.array([f[inv[i]] for i in range(8)])
fl0 = np.zeros(8); fl0[0] = 1; fl0[1] = -1
ft = np.zeros(8); ft[ta] = 1; ft[[neg(i) for i in ta]] = -1
for name, f in (('f_l0', fl0), ('f_t', ft)):
    ch = []
    for p in B:
        g = act_f(f, p); c = np.dot(g, f) / np.dot(f, f); assert np.allclose(g, c * f), name; ch.append(c)
    match = [a for a, d, z in irr if d == 1 and all(np.isclose(chi_of(a, p), c) for p, c in zip(B, ch))]
    print(f'V0 : {name} engendre une droite invariante ; irreductible {match}')
names = {}
for a, d, z in irr:
    names[a] = f'{"2" if d == 2 else "1"}[{a}]' + ('(spin)' if z == -1 else '')
print('irreductibles de B (indice, degre, z) :', irr)

def faces_of(F):
    out = set()
    for f in F:
        f = tuple(sorted(f))
        for r in range(1, len(f) + 1):
            for s in itertools.combinations(f, r): out.add(s)
    return out
def equiv_homology(F):
    fs = faces_of(F); byd = {}
    for s in fs: byd.setdefault(len(s) - 1, []).append(s)
    maxd = max(byd); byd[-1] = [()]
    ind = {d: {s: i for i, s in enumerate(sorted(byd[d]))} for d in byd}
    def bd(d):   # C_d -> C_{d-1}
        M = np.zeros((len(ind[d - 1]), len(ind[d])))
        for s, j in ind[d].items():
            for t in range(len(s)): M[ind[d - 1][s[:t] + s[t + 1:]], j] += (-1) ** t
        return M
    D = {d: bd(d) for d in range(0, maxd + 1)}
    def rho(p, d):
        g = vperm(p); n = len(ind[d]); M = np.zeros((n, n))
        for s, j in ind[d].items():
            im = [g.get(x, x) for x in s]; order = sorted(range(len(im)), key=lambda i: im[i])
            sgn = 1
            perm = order[:]
            for i in range(len(perm)):
                while perm[i] != i: j2 = perm[i]; perm[i], perm[j2] = perm[j2], perm[i]; sgn = -sgn
            M[ind[d][tuple(sorted(im))], j] = sgn
        return M
    R = {d: [rho(p, d) for p in B] for d in range(-1, maxd + 1)}
    res = {}
    for a, dg, z in irr:
        Q = {}
        for d in range(-1, maxd + 1):
            P = sum(np.conj(chi_of(a, p)) * Rm for p, Rm in zip(B, R[d])) * dg / len(B)
            u, s, vh = np.linalg.svd(P); Q[d] = u[:, s > 1e-8]
        for d in range(0, maxd + 1):
            A = D[d] @ Q[d]; kerdim = Q[d].shape[1] - (np.linalg.matrix_rank(A, 1e-8) if A.size else 0)
            nxt = D[d + 1] @ Q[d + 1] if d + 1 <= maxd else np.zeros((len(ind[d]), 0))
            h = kerdim - (np.linalg.matrix_rank(nxt, 1e-8) if nxt.size else 0)
            if h: res[(d, a)] = h // dg
    return res
def fmt(res): return ', '.join(f'{names[a]} en degre {d} (x{m})' for (d, a), m in sorted(res.items()))

base = {'S+M': S0 + [['l0', 'm']], 'S': S0}
print('controle S+M :', fmt(equiv_homology(base['S+M'])))
# orbites de sommets sous B
allv = sorted({x for f in base['S+M'] for x in f}); orbs = []
seen = set()
for v in allv:
    if v in seen: continue
    o = {vperm(p).get(v, v) for p in B}; seen |= o; orbs.append(sorted(o))
print('orbites de sommets :', orbs)
spin = [a for a, d, z in irr if z == -1]; d2 = [a for a in spin if irr[a][1] == 2][0]; singlets = [a for a in spin if irr[a][1] == 1]
hits = []
for bname, F in base.items():
    Vb = {x for f in F for x in f}; fsb = faces_of(F)
    for r in range(1, len(orbs) + 1):
        for W in itertools.combinations(orbs, r):
            Wset = set().union(*map(set, W))
            if not Wset <= Vb: continue
            Y = [list(s) for s in maximal([frozenset(s) for s in fsb if set(s) <= Wset])]
            for op in ('cone', 'susp'):
                G = F + [['N'] + y for y in Y] + ([['S'] + y for y in Y] + [['N', 'nN']] if op == 'susp' else [])
                res = equiv_homology(G)
                sp = {(d, a): m for (d, a), m in res.items() if a in spin}
                deg2 = [d for (d, a), m in sp.items() if a == d2 for _ in range(m)]
                degs = {a: [d for (d, b), m in sp.items() if b == a for _ in range(m)] for a in singlets}
                ok = len(deg2) == 1 and all(len(degs[a]) == 1 for a in singlets) and all((degs[a][0] - deg2[0]) % 2 == 1 for a in singlets)
                if ok:
                    zplus = sum(m * irr[a][1] for (d, a), m in res.items() if a not in spin)
                    hits.append((bname, op, [o[0] for o in W], len(Wset), zplus, fmt(res), G))
print(f'\nsucces : {len(hits)}')
hits.sort(key=lambda h: (h[4], h[3]))
for bname, op, W, nw, zp, f, G in hits[:12]:
    Vg = sorted({x for g in G for x in g}); au = aut_order(Vg, G)
    print(f'  {bname:4s} {op:4s} Y = orbites de {W} ({nw} sommets) ; |Aut| = {au} ; z=+1 : {zp} ; homologie : {f}', flush=True)
