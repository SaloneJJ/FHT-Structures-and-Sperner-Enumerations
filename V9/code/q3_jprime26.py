"""Q3 : J' explicite sur K26 (trois generations).
I = A^L * B^R * {g}, II = B^L * A^R * K4.  R : a^L*b^R*c  ->  b^L*a^R*phi(c),
phi(g_i) = (-1)^i d[face opposee a g_i] (dans le tetraedre des feuilles de generation)."""
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__)) + '/'

import contextlib, io, sys, itertools, time, json
sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(HERE + 'triplet.py').read())
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import lsqr
from simplicial import harmonic, chain_action, perm_sign
from fht_tools import build_K3, act3, actside, maximal as mx
SP = HERE
d = json.load(open(SP + 'V25_K3_pieces.json'))
A_ = [list(f) for f in mx(d['A'])]
B_ = [list(f) for f in mx([b for b in d['B'] if not any(v.startswith('y') for v in b)])]
K, V = build_K3(SP + 'V25_K3_pieces.json', with_pendants=False)
U = sorted(set(v for f in A_ + B_ for v in f))
rank = {}
for k, v in enumerate(U): rank['L_' + v] = k
for k, v in enumerate(U): rank['R_' + v] = 100 + k
for k in range(4): rank['g%d' % k] = 200 + k
ordK = lambda v: rank[v]
t0 = time.time()
faces = set()
for f in K:
    f = tuple(sorted(f, key=ordK))
    for k in range(1, len(f) + 1):
        for s in itertools.combinations(f, k): faces.add(s)
byd = {}
for s in faces: byd.setdefault(len(s) - 1, []).append(s)
for dd in byd: byd[dd].sort(key=lambda s: tuple(map(ordK, s)))
ind = {dd: {s: i for i, s in enumerate(byd[dd])} for dd in byd}
def bdry(dd):
    r, c, v = [], [], []
    for j, s in enumerate(byd[dd]):
        for t in range(len(s)):
            r.append(ind[dd - 1][s[:t] + s[t + 1:]]); c.append(j); v.append((-1) ** t)
    return sp.csr_matrix((v, (r, c)), shape=(len(byd[dd - 1]), len(byd[dd])))
Bd = {dd: bdry(dd) for dd in range(1, max(byd) + 1)}
print("K26 : simplexes", {dd: len(byd[dd]) for dd in sorted(byd)}, f"({time.time()-t0:.0f} s)", flush=True)
def chain_vec(dd, cdict):
    x = np.zeros(len(byd[dd]))
    for s, c in cdict.items(): x[ind[dd][s]] += c
    return x
def harm_part(dd, z):
    assert np.abs(Bd[dd] @ z).max() < 1e-9
    y = lsqr(Bd[dd + 1], z, atol=1e-14, btol=1e-14, iter_lim=40000)[0]
    return z - Bd[dd + 1] @ y
# harmonic cycles of the pieces
def actpiece(g, v): return actside(g, v)
bA, iA, HA = harmonic(A_, repr); bB, iB, HB = harmonic(B_, repr)
gB = {g: {v: actpiece(g, v) for v in set(x for f in B_ for x in f)} for g in G}
RBm = {g: chain_action(bB, iB, 1, gB[g], repr) for g in G}
def isoproj(Rg, H, k):
    Pm = sum(dim[k] / 24 * chi[k][idx[g]] * (H.T @ Rg[g] @ H) for g in G)
    w, Uu = np.linalg.eigh((Pm + Pm.T) / 2)
    return H @ Uu[:, w > 0.5]
B1 = {'1': isoproj(RBm, HB[1], '1'), '3': isoproj(RBm, HB[1], '3')}
Acl = {'1': (0, HA[0]), '2': (1, HA[1]), 's': (2, HA[2])}
def side_cycle(byd_, col, dd, side):
    D = {}
    for i in range(len(byd_[dd])):
        c = col[i]
        if abs(c) < 1e-12: continue
        img = [side + v for v in byd_[dd][i]]
        o = sorted(range(len(img)), key=lambda t: ordK(img[t]))
        key = tuple(img[t] for t in o); D[key] = D.get(key, 0) + c * perm_sign(o)
    return D
def join(parts):
    out = {(): 1.0}
    for Pp in parts:
        new = {}
        for s1, c1 in out.items():
            for s2, c2 in Pp.items(): new[s1 + s2] = new.get(s1 + s2, 0) + c1 * c2
        out = new
    return out
# generation : basis of H~0(4 points) (sum zero, orthonormal) and phi
Qg = np.linalg.qr(np.eye(4)[:, :3] - np.eye(4)[:, [3]])[0]
def gen0(m): return {('g%d' % i,): Qg[i, m] for i in range(4) if abs(Qg[i, m]) > 1e-12}
def face_bd(i):
    a, b, c = [j for j in range(4) if j != i]
    return {('g%d' % b, 'g%d' % c): 1.0, ('g%d' % a, 'g%d' % c): -1.0, ('g%d' % a, 'g%d' % b): 1.0}
def gen1(m):
    D = {}
    for i in range(4):
        for k, v in face_bd(i).items(): D[k] = D.get(k, 0) + (-1) ** i * Qg[i, m] * v
    return {k: v for k, v in D.items() if abs(v) > 1e-12}
classes = []
for ka, (dA, HAk) in Acl.items():
    for kb, HBk in B1.items():
        for a in range(HAk.shape[1]):
            for b in range(HBk.shape[1]):
                aL = side_cycle(bA, HAk[:, a], dA, 'L_'); aR = side_cycle(bA, HAk[:, a], dA, 'R_')
                bL = side_cycle(bB, HBk[:, b], 1, 'L_'); bR = side_cycle(bB, HBk[:, b], 1, 'R_')
                for m in range(3):
                    dI, dII = dA + 1 + 0 + 2, 1 + dA + 1 + 2
                    zI = chain_vec(dI, join([aL, bR, gen0(m)]))
                    zII = chain_vec(dII, join([bL, aR, gen1(m)]))
                    classes.append(((ka, kb), dI, harm_part(dI, zI), dII, harm_part(dII, zII)))
print("classes :", len(classes), f"({time.time()-t0:.0f} s)", flush=True)
deg_vecs = {}
for n, (lab, dI, hI, dII, hII) in enumerate(classes):
    deg_vecs.setdefault(dI, []).append(('I', n, hI)); deg_vecs.setdefault(dII, []).append(('II', n, hII))
print("dimensions par degre :", {dd: len(v) for dd, v in sorted(deg_vecs.items())},
      " independantes :", all(np.linalg.matrix_rank(np.column_stack([h for _, _, h in deg_vecs[dd]]), tol=1e-8) == len(deg_vecs[dd]) for dd in deg_vecs), flush=True)
N = len(classes)
Mdeg = {dd: np.column_stack([h for _, _, h in lst]) for dd, lst in deg_vecs.items()}
def chain_act_K(dd, gL, gR, gg):
    r, c, val = [], [], []
    for j, s in enumerate(byd[dd]):
        img = [act3(gL, gR, gg, v) for v in s]
        o = sorted(range(len(img)), key=lambda t: ordK(img[t]))
        r.append(ind[dd][tuple(img[t] for t in o)]); c.append(j); val.append(perm_sign(o))
    return sp.csr_matrix((val, (r, c)), shape=(len(byd[dd]),) * 2)
def rho_matrix(gL, gR, gg):
    Mx = np.zeros((2 * N, 2 * N))
    for dd, lst in deg_vecs.items():
        A = chain_act_K(dd, gL, gR, gg) @ Mdeg[dd]
        X, *_ = np.linalg.lstsq(Mdeg[dd], A, rcond=None)
        assert np.abs(Mdeg[dd] @ X - A).max() < 1e-7
        for col, (t, n, _) in enumerate(lst):
            for row, (t2, n2, _) in enumerate(lst):
                Mx[(0 if t2 == 'I' else N) + n2, (0 if t == 'I' else N) + n] = X[row, col]
    return Mx
Rm = np.zeros((2 * N, 2 * N))
for n in range(N): Rm[N + n, n] = 1; Rm[n, N + n] = 1
e4 = (0, 1, 2, 3); gens = [(1, 0, 2, 3), (1, 2, 3, 0)]
res = {}
for lab, (gL, gR, gg) in [('gauche (01)', (gens[0], e4, e4)), ('gauche (0123)', (gens[1], e4, e4)),
                          ('droite (01)', (e4, gens[0], e4)), ('gen (01)', (e4, e4, gens[0])), ('gen (0123)', (e4, e4, gens[1]))]:
    A1 = rho_matrix(gL, gR, gg); A2 = rho_matrix(gR, gL, gg)
    e_plus = np.abs(A2 @ Rm - Rm @ A1).max(); e_minus = np.abs(A2 @ Rm + Rm @ A1).max()
    res[lab] = (e_plus, e_minus)
    print(f"   {lab:14s} : |R rho - rho' R| = {e_plus:.1e} ; |R rho + rho' R| = {e_minus:.1e}", flush=True)
# Schur : R*R constant on each component (i,j) x gen
def gram_block(t, ns):
    out = np.zeros((len(ns), len(ns)))
    for a_, na in enumerate(ns):
        for b_, nb in enumerate(ns):
            ca, cb = classes[na], classes[nb]
            da = ca[1] if t == 'I' else ca[3]; db = cb[1] if t == 'I' else cb[3]
            if da != db: continue
            ha = ca[2] if t == 'I' else ca[4]; hb = cb[2] if t == 'I' else cb[4]
            out[a_, b_] = ha @ hb
    return out
for comp in [('1', '1'), ('1', '3'), ('2', '1'), ('2', '3'), ('s', '1'), ('s', '3')]:
    ns = [n for n, cl in enumerate(classes) if cl[0] == comp]
    ev = np.linalg.eigvals(np.linalg.solve(gram_block('I', ns), gram_block('II', ns))).real
    print(f"   composante {comp} x gen (dim {len(ns)}) : valeurs propres de R*R dans [{ev.min():.6f}, {ev.max():.6f}]", flush=True)
print(f"total {time.time()-t0:.0f} s")
