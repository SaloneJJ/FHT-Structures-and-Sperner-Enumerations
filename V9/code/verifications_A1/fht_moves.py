#!/usr/bin/env python3
"""fht_moves.py -- independent implementation of the moves of Article I (Fractal Hyper-Trees I, Section 5)
and of the checks reported in Section 8 ("Morphisms, classes and twisted algebras").

Written for the revision of Article I (November 2026) from the definitions of the article alone
(Definitions of fusion, fission, explosion, aggregation, leaf addition and leaf removal; twin blocks).
Pure Python 3, no dependency.

Representation.  A component is a Node: a dict  child-support -> child  (a child is a vertex, an int,
or a Node) and a set of host edges, each a frozenset of child supports.  The support of a child is the
frozenset of its vertices.  An FHT on V is its root Node (or the int x for the one-vertex FHT).

Library used by check_moves.py and census_s4.py.
"""
import itertools, random, sys, time, argparse
from collections import defaultdict

# ----------------------------------------------------------------------------- basic structure
def supp(c):
    return frozenset([c]) if isinstance(c, int) else c.V

class Node:
    __slots__ = ("ch", "E", "V")
    def __init__(self, ch, E):
        self.ch = dict(ch)          # support -> child
        self.E = set(E)             # host edges: frozensets of supports
        self.V = frozenset().union(*self.ch.keys())
    def copy(self):
        return Node({k: (v if isinstance(v, int) else v.copy()) for k, v in self.ch.items()}, set(self.E))

def own_edges(n):
    return {frozenset().union(*e) for e in n.E}

def realization(t):
    if isinstance(t, int):
        return frozenset()
    out = set(own_edges(t))
    for c in t.ch.values():
        out |= realization(c)
    return frozenset(out)

def components(t):
    """all non-leaf nodes with their distance to the root"""
    out, stack = [], [(t, 0)]
    while stack:
        n, r = stack.pop()
        if isinstance(n, int):
            continue
        out.append((n, r))
        for c in n.ch.values():
            stack.append((c, r + 1))
    return out

def level(c):
    return 0 if isinstance(c, int) else 1 + max(level(x) for x in c.ch.values())

# ----------------------------------------------------------------------------- reading a realization
def conn_components(points, edges):
    parent = {p: p for p in points}
    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for e in edges:
        e = list(e)
        for y in e[1:]:
            a, b = f(e[0]), f(y)
            if a != b: parent[a] = b
    cl = defaultdict(set)
    for p in points: cl[f(p)].add(p)
    return [frozenset(s) for s in cl.values()]

def tree_from_realization(V, E):
    """Reconstruction of Lemma 2.x (realization uniqueness): top edges = maximal edges, children of the
    root = connected components of the other edges.  Returns None if the family is not of this shape."""
    V = frozenset(V); E = set(E)
    if len(V) == 1:
        return next(iter(V)) if not E else None
    if not E: return None
    top = {e for e in E if not any(e < f for f in E)}
    rest = E - top
    kids = conn_components(V, rest)
    ch = {}
    for K in kids:
        if len(K) == 1:
            ch[K] = next(iter(K))
        else:
            sub = tree_from_realization(K, {e for e in rest if e <= K})
            if sub is None: return None
            ch[K] = sub
    hE = set()
    for e in top:
        he = frozenset(K for K in kids if K <= e)
        if not he or frozenset().union(*he) != e: return None
        hE.add(he)
    return Node(ch, hE)

# ----------------------------------------------------------------------------- literal axioms
def check_axioms(t, V):
    """Axioms (A1)-(A5) read literally on the family of components given by the tree, with
    level = 1 + max level of the children (the levels are recomputed, as in the article)."""
    V = frozenset(V)
    if isinstance(t, int):
        return V == frozenset([t])
    comps = []                      # (level, vertex set, own edges)
    def walk(n):
        if isinstance(n, int):
            comps.append((0, frozenset([n]), {frozenset([n])})); return 0
        l = 1 + max(walk(c) for c in n.ch.values())
        comps.append((l, n.V, own_edges(n))); return l
    N = walk(t)
    # A1
    if {c[1] for c in comps if c[0] == 0} != {frozenset([x]) for x in V}: return False
    for (l, S, E) in comps:
        if l == 0: continue
        # A2: connected Sperner hypergraph with vertex set S
        if any(e < f for e in E for f in E): return False
        if frozenset().union(*E) != S: return False
        if len(conn_components(S, E)) != 1: return False
    # A2: same level disjoint
    for a, b in itertools.combinations(comps, 2):
        if a[0] == b[0] and a[1] & b[1]: return False
    # A3
    tops = [c for c in comps if c[0] == N]
    if len(tops) != 1 or tops[0][1] != V: return False
    # A4
    for (l, S, E) in comps:
        for e in E:
            for (k, S2, _) in comps:
                if k < l and e & S2 and not (S2 < e): return False
    # A5
    for (l, S, E) in comps:
        if l >= 1 and not any(k == l - 1 and S2 <= S for (k, S2, _) in comps): return False
    return True

# ----------------------------------------------------------------------------- canonical form (isomorphism classes)
_CANON = {}
def canon(t):
    if isinstance(t, int):
        return "."
    key = (t.V, frozenset(frozenset(frozenset(k) for k in e) for e in t.E), frozenset((k, realization(c) if not isinstance(c, int) else c) for k, c in t.ch.items()))
    r = _CANON.get(key)
    if r is not None: return r
    r = _canon(t)
    if len(_CANON) > 2000000: _CANON.clear()
    _CANON[key] = r
    return r

def _canon(t):
    keys = list(t.ch.keys())
    cc = {k: canon(t.ch[k]) for k in keys}
    # isomorphism-invariant refinement of the colours (colour, sizes of the edges through k, then neighbours)
    inv = {k: (cc[k], tuple(sorted(len(e) for e in t.E if k in e))) for k in keys}
    for _ in range(2):
        inv = {k: (inv[k], tuple(sorted(tuple(sorted(inv[j] for j in e if j != k)) for e in t.E if k in e))) for k in keys}
    keys.sort(key=lambda k: inv[k])
    groups = [list(g) for _, g in itertools.groupby(keys, key=lambda k: inv[k])]
    best = None
    for perms in itertools.product(*[itertools.permutations(g) for g in groups]):
        order = [k for p in perms for k in p]
        pos = {k: i for i, k in enumerate(order)}
        es = tuple(sorted(tuple(sorted(pos[k] for k in e)) for e in t.E))
        if best is None or es < best: best = es
    return "(" + ",".join(cc[k] for k in keys) + "|" + ";".join("".join(map(chr, [65 + i for i in e])) for e in best) + ")"

# ----------------------------------------------------------------------------- admissible hosts, blocks
def admissible(X, E):
    X = frozenset(X)
    if not E: return False
    if any(len(e) < 2 for e in E): return False
    if any(e < f for e in E for f in E): return False
    if frozenset().union(*E) != X: return False
    return len(conn_components(X, E)) == 1

def blocks(n):
    sig = defaultdict(set)
    for k in n.ch:
        sig[frozenset(e for e in n.E if k in e)].add(k)
    return [frozenset(b) for b in sig.values()]

def subsets(s, kmin=0, kmax=None):
    s = list(s); kmax = len(s) if kmax is None else kmax
    for k in range(kmin, kmax + 1):
        for c in itertools.combinations(s, k):
            yield frozenset(c)

# ----------------------------------------------------------------------------- moves
# Every move returns (new_root, info).  The move is applied at the component `node` of the tree `t`
# (a path from the root identifies it); we rebuild the tree by copy.

def replace_node(t, target_V, new_node):
    """copy of t in which the node with support target_V is replaced by new_node (supports above unchanged)"""
    if isinstance(t, int): return t
    if t.V == target_V: return new_node
    ch = {}
    for k, c in t.ch.items():
        ch[k] = replace_node(c, target_V, new_node) if (not isinstance(c, int) and target_V <= c.V) else c
    return Node(ch, t.E)

def moves_th1p(t, n):
    X = list(n.ch); out = []
    for e in subsets(X, 2):
        if any(e <= f for f in n.E): continue
        A = {f for f in n.E if f <= e}
        if len(A) > 2: continue
        nn = Node(n.ch, (n.E - A) | {e})
        out.append((replace_node(t, n.V, nn), ("T1+", n.V, len(A))))
    return out

def moves_th1m(t, n):
    X = frozenset(n.ch); out = []
    for e in n.E:
        props = [b for b in subsets(e, 2) if b != e]
        for k in range(0, 3):
            for B in itertools.combinations(props, k):
                NE = (n.E - {e}) | set(B)
                if len(NE) != len(n.E) - 1 + k: continue      # B must be new distinct edges
                if admissible(X, NE):
                    nn = Node(n.ch, NE)
                    out.append((replace_node(t, n.V, nn), ("T1-", n.V, k)))
    return out

def pieces(c, eps):
    Y = frozenset(c.ch)
    return conn_components(Y, [f for f in c.E if f != eps])

def moves_th2m(t, n):
    out = []
    for ck, c in n.ch.items():
        if isinstance(c, int): continue
        for eps in c.E:
            P = pieces(c, eps)
            if len(P) < 2: continue
            newch = dict(n.ch); del newch[ck]
            bars = []
            for p in P:
                if len(p) == 1:
                    k = next(iter(p)); newch[k] = c.ch[k]; bars.append(k)
                else:
                    sub = Node({k: c.ch[k] for k in p}, {f for f in c.E if f != eps and f <= p})
                    newch[sub.V] = sub; bars.append(sub.V)
            NE = {(f - {ck}) | frozenset(bars) if ck in f else f for f in n.E}
            out.append((replace_node(t, n.V, Node(newch, NE)), ("T2-", n.V, ck, eps)))
    return out

def moves_th2p(t, n):
    out = []
    m = len(n.ch)
    for B in blocks(n):
        for S in subsets(B, 2, min(len(B), m - 1)):
            nonleaf = [s for s in S if not isinstance(n.ch[s], int)]
            for D in subsets(nonleaf):
                fixed = frozenset(s for s in S if s not in D)       # members kept as points
                choices = []
                for s in D:
                    sn = n.ch[s]
                    opts = [A for A in subsets(sn.ch, 1) if not any(f <= A for f in sn.E)]
                    choices.append(opts)
                for pick in itertools.product(*choices):
                    eps = fixed.union(*pick) if pick else fixed
                    pts = dict()
                    for s in S:
                        if s in D: pts.update(n.ch[s].ch)
                        else: pts[s] = n.ch[s]
                    cE = {eps}
                    for s in D: cE |= n.ch[s].E
                    c = Node(pts, cE)
                    newch = {k: v for k, v in n.ch.items() if k not in S}
                    newch[c.V] = c
                    NE = {(f - S) | {c.V} if S <= f else f for f in n.E}
                    out.append((replace_node(t, n.V, Node(newch, NE)), ("T2+", n.V, S, D, eps)))
    return out

def interior_moves(t, kinds=("T1+", "T1-", "T2-", "T2+"), maxrho=None):
    if isinstance(t, int): return []
    out = []
    for n, r in components(t):
        if maxrho is not None and r > maxrho: continue
        if "T1+" in kinds: out += [(x, i + (r,)) for x, i in moves_th1p(t, n)]
        if "T1-" in kinds: out += [(x, i + (r,)) for x, i in moves_th1m(t, n)]
        if "T2-" in kinds: out += [(x, i + (r,)) for x, i in moves_th2m(t, n)]
        if "T2+" in kinds: out += [(x, i + (r,)) for x, i in moves_th2p(t, n)]
    return out

def omega_plus(t, x):
    if isinstance(t, int):            # degenerate: F_1 -> St_2
        a = frozenset([t]); b = frozenset([x])
        return [(Node({a: t, b: x}, {frozenset([a, b])}), ("O+0",))]
    out = []
    for n, r in components(t):
        for B in blocks(n):
            xk = frozenset([x])
            nch = dict(n.ch); nch[xk] = x
            NE = {f | {xk} if B <= f else f for f in n.E}
            nn = Node(nch, NE)
            out.append((rebuild_up(t, n.V, nn), ("O+", n.V, B)))
    return out

def rebuild_up(t, target_V, new_node):
    """replace the node of support target_V by new_node, updating the supports of its ancestors"""
    if isinstance(t, int): return t
    if t.V == target_V: return new_node
    ch = {}; ren = {}
    for k, c in t.ch.items():
        if not isinstance(c, int) and target_V <= c.V:
            c2 = rebuild_up(c, target_V, new_node); ch[c2.V] = c2; ren[k] = c2.V
        else:
            ch[k] = c
    E = {frozenset(ren.get(k, k) for k in e) for e in t.E}
    return Node(ch, E)

def leaves_with_parent(t):
    out = []
    for n, r in components(t):
        for k, c in n.ch.items():
            if isinstance(c, int): out.append((c, n))
    return out

def omega_minus(t):
    """returns (allowed results, number of blocked removals)"""
    if isinstance(t, int): return [], 0
    allowed, blocked = [], 0
    for x, n in leaves_with_parent(t):
        xk = frozenset([x])
        bl = next(B for B in blocks(n) if xk in B)
        if len(bl) >= 2 and len(n.ch) >= 3:
            nch = {k: v for k, v in n.ch.items() if k != xk}
            NE = {f - {xk} for f in n.E}
            allowed.append((rebuild_up_remove(t, n.V, Node(nch, NE)), ("O-", x)))
        else:
            blocked += 1
    return allowed, blocked

def rebuild_up_remove(t, old_V, new_node):
    if isinstance(t, int): return t
    if t.V == old_V: return new_node
    ch = {}; ren = {}
    for k, c in t.ch.items():
        if not isinstance(c, int) and old_V <= c.V:
            c2 = rebuild_up_remove(c, old_V, new_node); ch[c2.V] = c2; ren[k] = c2.V
        else:
            ch[k] = c
    E = {frozenset(ren.get(k, k) for k in e) for e in t.E}
    return Node(ch, E)

def degenerate_minus(t):
    """St_2 -> F_1 (two choices of the removed leaf)"""
    if isinstance(t, int) or len(t.V) != 2: return []
    return [(min(t.V - {x}), ("O-0", x)) for x in t.V]

# ----------------------------------------------------------------------------- helpers
def star(V):
    V = sorted(V)
    return V[0] if len(V) == 1 else Node({frozenset([x]): x for x in V}, {frozenset(frozenset([x]) for x in V)})

def hierarchy(t):
    return frozenset(n.V for n, r in components(t))

def depth(t):
    return level(t)

def reduced_host_canon(n):
    bl = blocks(n)
    idx = {b: i for i, b in enumerate(bl)}
    bof = {k: idx[b] for b in bl for k in b}
    E = {frozenset(bof[k] for k in e) for e in n.E}
    best = None
    for p in itertools.permutations(range(len(bl))):
        es = tuple(sorted(tuple(sorted(p[i] for i in e)) for e in E))
        if best is None or es < best: best = es
    return best

def find_node(t, Vs):
    for n, r in components(t):
        if n.V == Vs: return n
    return None
