"""Outils : complexes simpliciaux par facettes, homologie harmonique avec caractere,
Betti mod p (creux), automorphismes d'hypergraphe (individualisation-raffinement)."""
import itertools, numpy as np
from collections import defaultdict

def maximal(facets):
    fs = sorted(set(frozenset(f) for f in facets), key=len, reverse=True)
    out = []
    for f in fs:
        if not any(f < g for g in out): out.append(f)
    return out

def all_faces(facets):
    faces = set()
    for f in facets:
        f = tuple(sorted(f, key=repr))
        for k in range(1, len(f) + 1):
            for s in itertools.combinations(f, k): faces.add(s)
    return faces

# ---------- dense harmonic homology with group character ----------
def harmonic(facets, order):
    """order : key function giving a total order on vertices"""
    faces = set()
    for f in facets:
        f = tuple(sorted(f, key=order))
        for k in range(1, len(f) + 1):
            for s in itertools.combinations(f, k): faces.add(s)
    maxd = max(len(s) for s in faces) - 1
    byd = {d: sorted([s for s in faces if len(s) == d + 1], key=lambda s: tuple(map(order, s))) for d in range(maxd + 1)}
    byd[-1] = [()]
    ind = {d: {s: i for i, s in enumerate(byd[d])} for d in byd}
    def bd(d):
        M = np.zeros((len(byd[d-1]), len(byd[d])))
        for j, s in enumerate(byd[d]):
            for t in range(len(s)): M[ind[d-1][s[:t] + s[t+1:]], j] = (-1) ** t
        return M
    B = {d: bd(d) for d in range(0, maxd + 1)}
    H = {}
    for d in range(-1, maxd + 1):
        L = (B[d+1] @ B[d+1].T if d + 1 <= maxd else 0) + (B[d].T @ B[d] if d >= 0 else 0)
        L = np.atleast_2d(L) if not np.isscalar(L) else np.zeros((len(byd[d]),) * 2)
        w, U = np.linalg.eigh(L)
        H[d] = U[:, w < 1e-8]
    return byd, ind, H

def perm_sign(seq):
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]: s = -s
    return s

def chain_action(byd, ind, d, vmap, order):
    M = np.zeros((len(byd[d]),) * 2)
    for j, s in enumerate(byd[d]):
        img = [vmap[v] for v in s]
        keys = [order(v) for v in img]
        o = sorted(range(len(img)), key=lambda t: keys[t])
        M[ind[d][tuple(img[t] for t in o)], j] = perm_sign(o)
    return M

def character_by_degree(facets, order, group_maps):
    """group_maps : list of vertex dicts. returns {d: array of traces}"""
    byd, ind, H = harmonic(facets, order)
    out = {}
    for d in H:
        if d < 0 or H[d].shape[1] == 0: continue
        out[d] = np.array([np.trace(H[d].T @ chain_action(byd, ind, d, g, order) @ H[d]) for g in group_maps])
    return out

# ---------- sparse Betti numbers mod p ----------
def betti_modp(facets, order, p=1000003):
    faces = set()
    for f in facets:
        f = tuple(sorted(f, key=order))
        for k in range(1, len(f) + 1):
            for s in itertools.combinations(f, k): faces.add(s)
    maxd = max(len(s) for s in faces) - 1
    byd = defaultdict(list)
    for s in faces: byd[len(s) - 1].append(s)
    byd[-1] = [()]
    ind = {d: {s: i for i, s in enumerate(byd[d])} for d in byd}
    def rank_bd(d):
        rows = []
        for s in byd[d]:
            row = {}
            for t in range(len(s)):
                row[ind[d-1][s[:t] + s[t+1:]]] = (-1) ** t % p
            rows.append(row)
        # gaussian elimination on sparse rows mod p
        pivots = {}
        r = 0
        for row in rows:
            row = dict(row)
            while row:
                c = min(row)
                if c in pivots:
                    prow = pivots[c]; f = row[c]
                    for k, v in prow.items():
                        nv = (row.get(k, 0) - f * v) % p
                        if nv: row[k] = nv
                        else: row.pop(k, None)
                else:
                    inv = pow(row[c], p - 2, p)
                    row = {k: v * inv % p for k, v in row.items()}
                    pivots[c] = row; r += 1
                    break
        return r
    rk = {d: rank_bd(d) for d in range(0, maxd + 1)}
    rk[maxd + 1] = 0
    betti = {d: len(byd[d]) - rk.get(d, 0) - rk[d + 1] for d in range(-1, maxd + 1)}
    return betti, {d: len(byd[d]) for d in range(-1, maxd + 1)}

def lefschetz(facets, order, vmap):
    """sum_k (-1)^k tr(g | C~_k), reduced (includes empty simplex in degree -1)"""
    faces = all_faces(facets)
    L = -1.0  # empty simplex, degree -1, trace +1, sign (-1)^-1
    for s in faces:
        img = tuple(vmap[v] for v in s)
        if set(img) == set(s):
            keys_s = sorted(s, key=order)
            pos = {v: i for i, v in enumerate(keys_s)}
            seq = [pos[vmap[v]] for v in keys_s]
            L += (-1) ** (len(s) - 1) * perm_sign(seq)
    return L

# ---------- automorphism group order of a hypergraph ----------
class Hyper:
    def __init__(self, vertices, facets):
        self.V = list(vertices); self.vi = {v: i for i, v in enumerate(self.V)}
        self.F = [frozenset(self.vi[v] for v in f) for f in facets]
        self.n, self.m = len(self.V), len(self.F)
        self.adj = [[] for _ in range(self.n + self.m)]
        for j, f in enumerate(self.F):
            for v in f:
                self.adj[v].append(self.n + j); self.adj[self.n + j].append(v)
        self.Fset = set(self.F)

    def refine(self, col):
        N = len(col)
        while True:
            sig = [(col[x], tuple(sorted(col[y] for y in self.adj[x % (self.n + self.m)] ))) for x in range(N)]
            # works on union copies: adj index modulo base size, but neighbors must be shifted
            raise NotImplementedError

def aut_order(vertices, facets, verbose=False):
    """|Aut| of the hypergraph (vertex permutations preserving the facet set),
    by individualization-refinement on the vertex/facet incidence graph."""
    V = list(vertices); vi = {v: i for i, v in enumerate(V)}
    F = [frozenset(vi[v] for v in f) for f in facets]
    n, m = len(V), len(F); N = n + m
    adj = [[] for _ in range(N)]
    for j, f in enumerate(F):
        for v in f: adj[v].append(n + j); adj[n + j].append(v)
    Fset = set(F)

    def refine_pair(c1, c2):
        """refine two colourings of the same graph jointly (shared colour names)"""
        while True:
            s1 = [(c1[x], tuple(sorted(c1[y] for y in adj[x]))) for x in range(N)]
            s2 = [(c2[x], tuple(sorted(c2[y] for y in adj[x]))) for x in range(N)]
            keys = sorted(set(s1) | set(s2))
            k = {s: i for i, s in enumerate(keys)}
            n1 = [k[s] for s in s1]; n2 = [k[s] for s in s2]
            if len(set(n1)) == len(set(c1)) and len(set(n2)) == len(set(c2)):
                return n1, n2
            c1, c2 = n1, n2

    def compatible(c1, c2):
        from collections import Counter
        return Counter(c1) == Counter(c2)

    def search(c1, c2):
        c1, c2 = refine_pair(c1, c2)
        if not compatible(c1, c2): return False
        # discrete on vertices?
        cells = defaultdict(list)
        for x in range(n): cells[c1[x]].append(x)
        cells2 = defaultdict(list)
        for x in range(n): cells2[c2[x]].append(x)
        nontriv = [c for c in cells if len(cells[c]) > 1]
        if not nontriv:
            mp = {cells[c][0]: cells2[c][0] for c in cells}
            return all(frozenset(mp[v] for v in f) in Fset for f in F)
        c = min(nontriv, key=lambda c: (len(cells[c]), c))
        v = cells[c][0]
        big = max(max(c1), max(c2)) + 1
        for w in cells2[c]:
            d1 = list(c1); d2 = list(c2)
            d1[v] = big; d2[w] = big
            if search(d1, d2): return True
        return False

    base = [0] * n + [1] * m
    order = 1
    fixed = []
    col, _ = refine_pair(base, base)
    stab_col = col
    for v in range(n):
        # orbit of v under stabilizer of fixed points
        cand = [w for w in range(n) if stab_col[w] == stab_col[v]]
        if len(cand) == 1:
            continue
        big = max(stab_col) + 1
        orbit = 0
        for w in cand:
            d1 = list(stab_col); d2 = list(stab_col)
            d1[v] = big; d2[w] = big
            if search(d1, d2): orbit += 1
        order *= orbit
        if verbose: print("   fix", V[v], "orbite", orbit)
        d1 = list(stab_col); d1[v] = max(stab_col) + 1
        stab_col, _ = refine_pair(d1, d1)
    return order
