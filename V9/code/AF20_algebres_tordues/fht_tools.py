"""FHT_21 / code_v27 : outils communs (session du 2026-11-08).
- complexes par facettes, Betti sur Z/2 rapides (masques binaires), automorphismes d'hypergraphe,
- action de S4^3 sur les feuilles de K3 / K26,
- mouvements elementaires sur un FHT plat (V28, P0).
Python 3, numpy seulement pour rien d'essentiel ; aucune dependance externe."""
import itertools, json
from collections import defaultdict

# ---------------------------------------------------------------- complexes
def maximal(facets):
    fs = sorted(set(frozenset(f) for f in facets), key=len, reverse=True)
    out = []
    for f in fs:
        if not any(f < g for g in out): out.append(f)
    return out

def betti_z2(facets):
    """Betti reduits sur Z/2 du complexe engendre par les facettes (listes de sommets hashables)."""
    faces = set()
    for f in facets:
        f = tuple(sorted(f, key=str))
        for k in range(1, len(f) + 1):
            for s in itertools.combinations(f, k): faces.add(s)
    byd = defaultdict(list)
    for s in faces: byd[len(s) - 1].append(s)
    byd[-1] = [()]
    maxd = max(byd)
    ind = {d: {s: i for i, s in enumerate(sorted(byd[d], key=lambda x: tuple(map(str, x))))} for d in byd}
    def rank(d):
        piv = {}
        r = 0
        for s in ind[d]:
            col = 0
            for t in range(len(s)):
                col ^= 1 << ind[d - 1][s[:t] + s[t + 1:]]
            while col:
                h = col.bit_length() - 1
                if h in piv: col ^= piv[h]
                else: piv[h] = col; r += 1; break
        return r
    rk = {d: rank(d) for d in range(0, maxd + 1)}
    rk[maxd + 1] = 0; rk[-1] = 0
    return {d: len(ind[d]) - rk[d] - rk[d + 1] for d in range(-1, maxd + 1) if len(ind[d]) - rk[d] - rk[d + 1]}

def connected_cover(facets, V):
    V = set(V)
    if set().union(*facets) != V: return False
    adj = defaultdict(set)
    for f in facets:
        f = list(f)
        for a in f: adj[a] |= set(f)
    start = next(iter(V)); seen = {start}; st = [start]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return seen == V

def is_clutter(facets):
    fs = list(facets)
    for i, a in enumerate(fs):
        for b in fs[i + 1:]:
            if a <= b or b <= a: return False
    return True

# ---------------------------------------------------------------- automorphismes (individualisation-raffinement)
def aut_order(vertices, facets):
    V = list(vertices); vi = {v: i for i, v in enumerate(V)}
    F = [frozenset(vi[v] for v in f) for f in facets]
    n, m = len(V), len(F); N = n + m
    adj = [[] for _ in range(N)]
    for j, f in enumerate(F):
        for v in f: adj[v].append(n + j); adj[n + j].append(v)
    Fset = set(F)
    def refine_pair(c1, c2):
        while True:
            s1 = [(c1[x], tuple(sorted(c1[y] for y in adj[x]))) for x in range(N)]
            s2 = [(c2[x], tuple(sorted(c2[y] for y in adj[x]))) for x in range(N)]
            keys = sorted(set(s1) | set(s2)); k = {s: i for i, s in enumerate(keys)}
            n1 = [k[s] for s in s1]; n2 = [k[s] for s in s2]
            if len(set(n1)) == len(set(c1)) and len(set(n2)) == len(set(c2)): return n1, n2
            c1, c2 = n1, n2
    from collections import Counter
    def search(c1, c2):
        c1, c2 = refine_pair(c1, c2)
        if Counter(c1) != Counter(c2): return False
        cells = defaultdict(list); cells2 = defaultdict(list)
        for x in range(n): cells[c1[x]].append(x); cells2[c2[x]].append(x)
        nontriv = [c for c in cells if len(cells[c]) > 1]
        if not nontriv:
            mp = {cells[c][0]: cells2[c][0] for c in cells}
            return all(frozenset(mp[v] for v in f) in Fset for f in F)
        c = min(nontriv, key=lambda c: (len(cells[c]), c)); v = cells[c][0]
        big = max(max(c1), max(c2)) + 1
        for w in cells2[c]:
            d1 = list(c1); d2 = list(c2); d1[v] = big; d2[w] = big
            if search(d1, d2): return True
        return False
    base = [0] * n + [1] * m
    col, _ = refine_pair(base, base); order = 1
    for v in range(n):
        cand = [w for w in range(n) if col[w] == col[v]]
        if len(cand) == 1: continue
        big = max(col) + 1; orbit = 0
        for w in cand:
            d1 = list(col); d2 = list(col); d1[v] = big; d2[w] = big
            if search(d1, d2): orbit += 1
        order *= orbit
        d1 = list(col); d1[v] = max(col) + 1
        col, _ = refine_pair(d1, d1)
    return order

# ---------------------------------------------------------------- S4^3 sur les feuilles de K3 / K26
PAIR = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
S4 = list(itertools.permutations(range(4)))
def s3img(g):
    out = []
    for P in PAIR:
        img = frozenset(frozenset(g[x] for x in b) for b in P)
        out.append([frozenset(frozenset(b) for b in Q) for Q in PAIR].index(img))
    return out
def actside(g, name):
    if name[0] in 'xy': return name[0] + str(g[int(name[1])])
    if name[0] == 'p': return 'p' + str(s3img(g)[int(name[1])])
    return name
def act3(gL, gR, gg, v):
    if v.startswith('L_'): return 'L_' + actside(gL, v[2:])
    if v.startswith('R_'): return 'R_' + actside(gR, v[2:])
    if v.startswith('g'): return 'g' + str(gg[int(v[1])])
    return v

def build_K3(pieces_json, with_pendants=True):
    d = json.load(open(pieces_json))
    A_, B_ = d["A"], d["B"]
    if not with_pendants: B_ = [b for b in B_ if not any(v.startswith('y') for v in b)]
    G4, E4 = d["generation"]["feuilles"], d["generation"]["K4"]
    I = [frozenset(['L_' + v for v in a] + ['R_' + v for v in b] + [g]) for a in A_ for b in B_ for g in G4]
    II = [frozenset(['L_' + v for v in b] + ['R_' + v for v in a] + e) for b in B_ for a in A_ for e in E4]
    K = maximal(I + II)
    return K, sorted(set(v for f in K for v in f))

def group_perms_K3(V):
    """liste des 13824 permutations des feuilles (dict) pour S4^3"""
    out = []
    for gL in S4:
        for gR in S4:
            for gg in S4:
                out.append({v: act3(gL, gR, gg, v) for v in V})
    return out

# ---------------------------------------------------------------- mouvements elementaires (V28, P0)
def neighbours(E, V, moves=('suppr', 'ajout', 'agrand', 'reduc', 'scission', 'jumeau')):
    """genere (nom, donnees, nouvelles_aretes, nouvelles_feuilles) ; E = liste de frozensets"""
    Es = set(E)
    adjpairs = set()
    for f in E:
        for a, b in itertools.combinations(sorted(f), 2): adjpairs.add((a, b))
    def valid(newE, newV):
        return is_clutter(newE) and connected_cover(newE, newV)
    for e in E:
        rest = [f for f in E if f != e]
        if 'suppr' in moves:
            if valid(rest, V): yield ('suppr', (e,), rest, V)
        if 'reduc' in moves and len(e) >= 3:
            for x in sorted(e):
                b = e - {x}
                nE = rest + [b]
                if valid(nE, V): yield ('reduc', (e, x), nE, V)
        if 'scission' in moves and len(e) >= 3:
            for x, y in itertools.combinations(sorted(e), 2):
                nE = rest + [e - {x}, e - {y}]
                if valid(nE, V): yield ('scission', (e, x, y), nE, V)
        if 'agrand' in moves:
            for x in V:
                if x in e: continue
                ne = e | {x}
                absorbed = [f for f in E if f <= ne]
                if len(absorbed) > 2: continue
                nE = [f for f in E if f not in absorbed] + [ne]
                if valid(nE, V): yield ('agrand', (e, x), nE, V)
    if 'ajout' in moves:
        for a, b in itertools.combinations(sorted(V), 2):
            if (a, b) in adjpairs: continue
            nE = E + [frozenset([a, b])]
            if valid(nE, V): yield ('ajout', (a, b), nE, V)
    if 'jumeau' in moves:
        for x in V:
            t = x + "'"
            nE = [f | {t} if x in f else f for f in E]
            yield ('jumeau', (x,), nE, V + [t])
