"""Verification des hypotheses du theoreme de rigidite (V24) sur K25, K33, K35."""
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__)) + '/'

import contextlib, io, sys, itertools
sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(HERE + 'bis4_pieces.py').read())

P = [('p', k) for k in range(3)]
Astar = [Z, ['f', P[0], P[1]], ['f', P[0], P[2]], ['f', P[1], P[2]], [P[0], P[1], P[2]]]
for k, ((a, b), (c, d)) in enumerate(PAIRINGS):
    Astar += [[('p', k), 'z', X[a], X[b]], [('p', k), 'z', X[c], X[d]]]
Astar += [['q']]
Bp = [Z] + [[v, 'w'] for v in Z]
Bpy = Bp + [[X[i], ('y', i)] for i in range(4)]
cases = {'K25': (Astar, Bp), 'K33': (Astar, Bpy), 'K35': (A, B)}

def faces_of(F):
    return set(all_faces([sorted(f, key=repr) for f in F]))
def minimal_nonfaces(F):
    V = sorted(set(v for f in F for v in f), key=repr)
    fc = set(frozenset(s) for s in faces_of(F))
    mnf = []
    for k in range(2, len(V) + 1):
        for S in itertools.combinations(V, k):
            S = frozenset(S)
            if S in fc: continue
            if all(S - {v} in fc for v in S): mnf.append(S)
        if k > 6: break   # facets have <= 5 vertices : minimal non-faces have <= 6
    return V, mnf
def join_irreducible(F):
    V, mnf = minimal_nonfaces(F)
    covered = set().union(*mnf) if mnf else set()
    # components of the minimal-non-face hypergraph
    parent = {v: v for v in V}
    def find(v):
        while parent[v] != v: parent[v] = parent[parent[v]]; v = parent[v]
        return v
    for S in mnf:
        S = list(S)
        for v in S[1:]: parent[find(v)] = find(S[0])
    comps = len(set(find(v) for v in V))
    return covered == set(V) and comps == 1, len(mnf), comps, set(V) - covered

def aut_pair(Af, Bf):
    """|Aut(U; A, B)| : permutations of V(U) preserving A and B (marker vertices)"""
    Fa = [list(f) + ['#A'] for f in maximal(Af)]
    Fb = [list(f) + ['#B'] for f in maximal(Bf)]
    Fm = Fa + Fb + [['#A', '#A2'], ['#A2', '#A3']]   # chemin pendant : distingue #A de #B sans jumeau
    V = sorted(set(v for f in Fm for v in f), key=repr)
    return aut_order(V, Fm)

for name, (Af, Bf) in cases.items():
    Am, Bm = [list(f) for f in maximal(Af)], [list(f) for f in maximal(Bf)]
    VA = set(v for f in Am for v in f); VB = set(v for f in Bm for v in f)
    ia, na, ca, ua = join_irreducible(Am); ib, nb, cb, ub = join_irreducible(Bm)
    ZA = set(frozenset(s) for s in faces_of(Am)) & set(frozenset(s) for s in faces_of(Bm))
    print(f"{name}:")
    print(f"   H1 irreductibles : A {ia} ({na} non-faces minimales), B {ib} ({nb}) ; A != B (homologies differentes)")
    print(f"   H2 A ∩ B = Z : {ZA == set(frozenset(s) for s in faces_of([Z]))} ; V(A) ∩ V(B) = V(Z) : {VA & VB == set(Z)} ;"
          f" A != Z : {VA != set(Z)} ; B != Z : {VB != set(Z)}")
    print(f"   H3 aretes de B de taille >= 2 : {min(map(len, Bm)) >= 2}")
    print(f"   H4 |Aut(A)| = {aut_order(sorted(VA, key=repr), Am)}, |Aut(U; A, B)| = {aut_pair(Am, Bm)}")
