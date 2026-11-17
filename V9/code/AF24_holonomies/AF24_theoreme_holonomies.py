#!/usr/bin/env python3
"""AF24 -- verification du theoreme des holonomies (version 2, 14/11/2026, soir).
Enonces verifies (02_theoreme_holonomies/THEOREME_HOLONOMIES_V2.md) :
  L1  [algebre]  pour tout g = (pi, s) de O_h : V(g) G5 V(g)^+ = G5 et V(g) FG V(g)^+ = sgn(pi) FG ; V(g) G4 V(g)^+ = det(g) G4.
  L1' [AF23]     sur le bord du 4-cube, les passages ou FG change de signe (transport de AF23) sont exactement ceux ou la permutation des axes est impaire.
  L2  [arete]    W0 = faces a permutation impaire ; dW0 = Omega = aretes de degre impair (mod 2).
  T   [obstr.]   masse de gout sans paroi  <=>  W0 est un bord (W0 dans B2)  ; si Omega = 0, obstruction = [W0] dans H2(M ; Z/2).
Complexes : bord du 4-cube (raffinement 2^n, comme AF23), S^3 « oreiller » (double du bloc 2^n), variete plate a quart de tour (G4), tore T^3 (controle).
Duree : < 1 min.  Usage : python3 AF24_theoreme_holonomies.py | tee AF24.log
"""
import sys, os, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, '..', 'AF23_dirac_S3'), os.path.join(HERE, '..', 'V352')): sys.path.insert(0, p)

def sgn(p):
    p = list(p); s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s

# ---------------------------------------------------------------- L1, L1'
def test_algebre():
    import AF23_dirac_natif_S3 as A
    ok = True; n = 0
    for pi in itertools.permutations(range(3)):
        for s in itertools.product((1, -1), repeat=3):
            V = A.transport(None, list(pi), list(s)); Vd = V.conj().T; n += 1
            ok &= np.allclose(V @ A.G5 @ Vd, A.G5)
            ok &= np.allclose(V @ A.FG @ Vd, sgn(pi) * A.FG)
            ok &= np.allclose(V @ A.G4 @ Vd, sgn(pi) * s[0] * s[1] * s[2] * A.G4)
    print(f'L1  : {n} elements de O_h ; G5 invariant, FG -> sgn(pi) FG, G4 -> det(g) G4 : {ok}')
    # L1' : passages de AF23 (meme construction que AF23_diag_holonomie)
    facets = [(i, s) for i in range(4) for s in (0, 1)]; tang = {f: [a for a in range(4) if a != f[0]] for f in facets}
    agree = True; nflip = 0; npas = 0
    for f in facets:
        i, s = f; tA = tang[f]
        for d in range(3):
            for e in (1, -1):
                j = tA[d]; g = (j, 1 if e > 0 else 0); tB = tang[g]
                sj = e * (1 if s == 0 else -1)
                cmap = [tB.index(i) if a == j else tB.index(a) for a in tA]; sig = [sj if a == j else 1 for a in tA]
                V = A.transport(None, cmap, sig); flip = np.allclose(V @ A.FG @ V.conj().T, -A.FG)
                agree &= (flip == (sgn(cmap) == -1)); nflip += flip; npas += 1
    print(f"L1' : bord du 4-cube : FG change de signe a {nflip}/{npas} passages ; flip <=> permutation impaire : {agree}")
    return ok and agree

# ---------------------------------------------------------------- complexes cubiques combinatoires
class Cx:
    """cubes : liste de dicts {corner(bits) -> id global de sommet} ; repere local = ordre des bits (axes 0, 1, 2)"""
    def __init__(self, name, cubes):
        self.name = name; self.cubes = cubes
        self.E = {}; self.F = {}; self.Einc = []; self.Finc = []
        for q, cor in enumerate(cubes):
            for a in range(3):                       # aretes locales : direction a
                for b in itertools.product((0, 1), repeat=2):
                    c0 = list(b); c0.insert(a, 0); c1 = list(b); c1.insert(a, 1)
                    key = frozenset((cor[tuple(c0)], cor[tuple(c1)]))
                    assert len(key) == 2, 'arete en boucle'
                    self.Einc.append((q, a, self.E.setdefault(key, len(self.E))))
            for a in range(3):                       # faces locales : normale a, cote t
                for t in (0, 1):
                    vs = frozenset(cor[c] for c in itertools.product((0, 1), repeat=3) if c[a] == t)
                    assert len(vs) == 4
                    self.Finc.append((q, a, t, self.F.setdefault(vs, len(self.F))))
        self.V = sorted({v for cor in cubes for v in cor.values()})
        self.vidx = {v: k for k, v in enumerate(self.V)}
        self.deg = np.zeros(len(self.E), int)
        for _, _, e in self.Einc: self.deg[e] += 1
        self.fcubes = {}
        for q, a, t, f in self.Finc: self.fcubes.setdefault(f, []).append((q, a, t))
        assert all(len(v) == 2 for v in self.fcubes.values()), 'face non partagee par exactement deux cotes de cubes'

    def face_edges(self, q, a, t):
        cor = self.cubes[q]; out = []
        for c in range(3):
            if c == a: continue
            for b in (0, 1):
                c0 = [0, 0, 0]; c1 = [0, 0, 0]; c0[a] = c1[a] = t
                o = [x for x in range(3) if x not in (a, c)][0]; c0[o] = c1[o] = b; c1[c] = 1
                out.append((c, frozenset((cor[tuple(c0)], cor[tuple(c1)]))))
        return out

    def transition(self, f):
        """permutation des axes locaux du cube A vers ceux du cube B a travers la face f"""
        (qA, aA, tA), (qB, aB, tB) = self.fcubes[f]
        eB = {k: c for c, k in self.face_edges(qB, aB, tB)}
        pi = [None] * 3; pi[aA] = aB
        for c, k in self.face_edges(qA, aA, tA): pi[c] = eB[k]
        assert sorted(pi) == [0, 1, 2]
        return qA, qB, pi

    def boundaries(self):
        nV, nE, nF, nC = len(self.V), len(self.E), len(self.F), len(self.cubes)
        d1 = np.zeros((nV, nE), np.uint8); d2 = np.zeros((nE, nF), np.uint8); d3 = np.zeros((nF, nC), np.uint8)
        for key, e in self.E.items():
            for v in key: d1[self.vidx[v], e] ^= 1
        seen = set()
        for q, a, t, f in self.Finc:
            if f in seen: continue
            seen.add(f)
            for _, k in self.face_edges(q, a, t): d2[self.E[k], f] ^= 1
        for q, a, t, f in self.Finc: d3[f, q] ^= 1
        return d1, d2, d3

def rank2(M):
    M = M.copy() % 2; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0: continue
        p = r + piv[0]; M[[r, p]] = M[[p, r]]
        nz = np.nonzero(M[:, c])[0]; nz = nz[nz != r]
        M[nz] ^= M[r]; r += 1
        if r == rows: break
    return r

def solve2(M, b):
    """une solution x de M x = b (mod 2) ou None"""
    A = np.concatenate([M % 2, (b % 2)[:, None]], axis=1).astype(np.uint8); rows, cols = M.shape; r = 0; piv = []
    for c in range(cols):
        p = np.nonzero(A[r:, c])[0]
        if len(p) == 0: continue
        p = r + p[0]; A[[r, p]] = A[[p, r]]
        nz = np.nonzero(A[:, c])[0]; nz = nz[nz != r]; A[nz] ^= A[r]; piv.append(c); r += 1
        if r == rows: break
    if A[r:, -1].any(): return None
    x = np.zeros(cols, np.uint8)
    for i, c in enumerate(piv): x[c] = A[i, -1]
    return x

def link_ok(C):
    """liens des sommets : surfaces fermees connexes de caracteristique 2 (variete)"""
    ok = True
    for v in C.V:
        Ev = {e for k, e in C.E.items() if v in k}
        tri = [(q, c) for q, cor in enumerate(C.cubes) for c, w in cor.items() if w == v]
        Fv = 0; adj = {}
        for q, a, t, f in C.Finc:
            cor = C.cubes[q]
            for c, w in cor.items():
                if w == v and c[a] == t: adj.setdefault((f, w), []).append((q, c))
        Fv = len(adj)
        if len(Ev) - Fv + len(tri) != 2: ok = False; continue
        parent = {x: x for x in tri}
        def find(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for lst in adj.values():
            if len(lst) != 2: ok = False; break
            parent[find(lst[0])] = find(lst[1])
        ok &= len({find(x) for x in tri}) == 1
    return ok

def analyse(C, enumerate_min=False):
    print(f'\n== {C.name} : {len(C.V)} sommets, {len(C.E)} aretes, {len(C.F)} faces, {len(C.cubes)} cubes')
    from collections import Counter
    print('   degres des aretes :', dict(sorted(Counter(C.deg.tolist()).items())))
    d1, d2, d3 = C.boundaries()
    assert not ((d1.astype(int) @ d2) % 2).any() and not ((d2.astype(int) @ d3) % 2).any()
    r1, r2, r3 = rank2(d1), rank2(d2), rank2(d3)
    b = (len(C.V) - r1, len(C.E) - r1 - r2, len(C.F) - r2 - r3, len(C.cubes) - r3)
    print('   Betti mod 2 :', b, '; liens de sommets = spheres :', link_ok(C))
    W0 = np.zeros(len(C.F), np.uint8); cube_pairs = []
    for f in C.fcubes:
        qA, qB, pi = C.transition(f); odd = sgn(pi) == -1; W0[f] = odd; cube_pairs.append((qA, qB, odd))
    Omega = (C.deg % 2).astype(np.uint8)
    dW0 = (d2.astype(int) @ W0) % 2
    print(f'   |W0| (faces a permutation impaire) = {int(W0.sum())} ; |Omega| (aretes impaires) = {int(Omega.sum())} ; L2 dW0 = Omega : {np.array_equal(dW0, Omega)}')
    # T : sans paroi <=> W0 = d3 U
    U = solve2(d3, W0); eps_bfs = bfs_signs(len(C.cubes), cube_pairs)
    print(f'   masse de gout sans paroi : {U is not None} (controle par propagation des signes : {eps_bfs})')
    if not Omega.any():
        print('   Omega = 0 : obstruction residuelle = [W0] dans H2(M ; Z/2) ; [W0] = 0 :', U is not None)
    if enumerate_min:
        x0 = solve2(d2, Omega); assert x0 is not None
        # noyau de d2 : base par elimination
        K = kernel2(d2); best = None; cnt = 0
        for bits in itertools.product((0, 1), repeat=K.shape[0]):
            w = (x0 + (np.array(bits, int) @ K)) % 2 if K.shape[0] else x0
            s = int(w.sum())
            if best is None or s < best: best, cnt = s, 1
            elif s == best: cnt += 1
        print(f'   parois de bord Omega : 2^{K.shape[0]} choix ; aire minimale {best} faces, atteinte par {cnt} choix ; |W0| = {int(W0.sum())}')
    return U is not None

def kernel2(M):
    A = (M % 2).astype(np.uint8).copy(); rows, cols = A.shape; r = 0; piv = []
    for c in range(cols):
        p = np.nonzero(A[r:, c])[0]
        if len(p) == 0: continue
        p = r + p[0]; A[[r, p]] = A[[p, r]]; nz = np.nonzero(A[:, c])[0]; nz = nz[nz != r]; A[nz] ^= A[r]; piv.append(c); r += 1
        if r == rows: break
    free = [c for c in range(cols) if c not in piv]; K = []
    for fc in free:
        v = np.zeros(cols, np.uint8); v[fc] = 1
        for i, c in enumerate(piv): v[c] = A[i, fc]
        K.append(v)
    return np.array(K, np.uint8).reshape(len(K), cols)

def bfs_signs(n, pairs):
    adj = {}
    for a, b, odd in pairs: adj.setdefault(a, []).append((b, odd)); adj.setdefault(b, []).append((a, odd))
    eps = {0: 1}; st = [0]
    while st:
        a = st.pop()
        for b, odd in adj.get(a, []):
            w = -eps[a] if odd else eps[a]
            if b in eps:
                if eps[b] != w: return False
            else: eps[b] = w; st.append(b)
    return len(eps) == n

# ---------------------------------------------------------------- constructions
def bord_4cube(n):
    m = 2 ** n; cubes = []
    for i in range(4):
        for s in (0, 1):
            tA = [a for a in range(4) if a != i]
            for co in itertools.product(range(m), repeat=3):
                cor = {}
                for c in itertools.product((0, 1), repeat=3):
                    x = [0] * 4; x[i] = s * m
                    for k, a in enumerate(tA): x[a] = co[k] + c[k]
                    cor[c] = tuple(x)
                cubes.append(cor)
    return Cx(f'bord du 4-cube, niveau {n} (AF23)', cubes)

def oreiller(n):
    m = 2 ** n; cubes = []
    def vid(cp, x): return ('b',) + x if any(t in (0, m) for t in x) else ('i', cp) + x
    for cp in (0, 1):
        for co in itertools.product(range(m), repeat=3):
            cubes.append({c: vid(cp, tuple(co[k] + c[k] for k in range(3))) for c in itertools.product((0, 1), repeat=3)})
    return Cx(f'S^3 oreiller (double du bloc {m}^3), niveau {n}', cubes)

def tore(m):
    cubes = [{c: tuple((co[k] + c[k]) % m for k in range(3)) for c in itertools.product((0, 1), repeat=3)} for co in itertools.product(range(m), repeat=3)]
    return Cx(f'tore T^3 ({m}^3)', cubes)

def quart_de_tour(m):
    def vid(x):
        x, y, z = x
        if z == m: x, y, z = -y, x, 0               # (x, y, m) ~ (R(x, y), 0), R = rotation de 90 degres
        return (x % m, y % m, z)
    cubes = [{c: vid(tuple(co[k] + c[k] for k in range(3))) for c in itertools.product((0, 1), repeat=3)} for co in itertools.product(range(m), repeat=3)]
    return Cx(f'variete plate a quart de tour G4 ({m}^3)', cubes)

if __name__ == '__main__':
    print(__doc__.strip().splitlines()[0])
    test_algebre()
    analyse(bord_4cube(0), enumerate_min=True)
    for n in (1, 2): analyse(bord_4cube(n))
    for n in (0, 1, 2): analyse(oreiller(n))
    analyse(tore(3)); analyse(quart_de_tour(3)); analyse(quart_de_tour(4))
