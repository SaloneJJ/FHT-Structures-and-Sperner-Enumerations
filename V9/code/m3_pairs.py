"""M3 : classes d'isomorphisme des objets de type K26 = (A^L*B^R*{g}) u (B^L*A^R*K4(g))
pour A parcourant les 6 classes de pieces A rigides minimales (V23, 18 pieces) et B parcourant
toutes les pieces S4-invariantes sur Z+w a homologie 1+3 concentree en un degre.
Necessite landscape.pkl (produit par landscape.py). Ecrit m3_pairs.json."""
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__)) + '/'

import contextlib, io, sys, itertools, pickle, json
SP = HERE + ''
sys.path.insert(0, SP)
src = open(SP + 'landscape.py').read().split("# ---------- B minimal")[0]
with contextlib.redirect_stdout(io.StringIO()):
    exec(src)
import fht_tools as ft
def nm(v): return v[0] + str(v[1]) if isinstance(v, tuple) else str(v)
# --- A : 18 pieces, classes (m3_classes)
res = pickle.load(open(SP + 'landscape.pkl', 'rb')); orbsA, found = res['Z+X3+f+q']
Apieces = []
for r, combo, c, a, nv in found:
    if a != 24: continue
    Apieces.append(ft.maximal([frozenset(['z', 'x0', 'x1', 'x2', 'x3'])] + [frozenset(nm(v) for v in s) for i in combo for s in orbsA[i]] + [frozenset(['q'])]))
def iso(F1, F2, drop=frozenset(['q'])):
    U = [frozenset('a' + v for v in f) for f in F1 if f != drop] + [frozenset('b' + v for v in f) for f in F2 if f != drop]
    V = sorted(set(v for f in U for v in f)); a1 = ft.aut_order(sorted(set(v for f in U if next(iter(f))[0]=='a' for v in f)), [f for f in U if next(iter(f))[0]=='a'])
    return ft.aut_order(V, U) == 2 * a1 * a1
clsA = []
for i, F in enumerate(Apieces):
    for c in clsA:
        if iso(Apieces[c[0]], F): c.append(i); break
    else: clsA.append([i])
Areps = [Apieces[c[0]] for c in clsA]
# --- B : toutes les pieces sur Z+w, H~ = 1+3 dans un seul degre
V = Z + ['w']
seen, orbs = set(), []
for k in range(1, len(V) + 1):
    for S in itertools.combinations(V, k):
        S = frozenset(S)
        if S in seen: continue
        orb = frozenset(frozenset(act(g, v) for v in S) for g in G)
        seen |= orb; orbs.append(sorted(orb, key=lambda s: sorted(map(repr, s))))
Bs = []
for r in range(1, 5):
    for combo in itertools.combinations(range(len(orbs)), r):
        F = [list(s) for i in combo for s in orbs[i]]
        Fm = maximal(F)
        if len(Fm) != len(set(frozenset(f) for f in F)): continue
        if set(v for f in F for v in f) != set(V): continue
        c = content(F)
        nz = [(d, v) for d, v in c.items() if any(v)]
        if len(nz) == 1 and nz[0][1] == [1, 0, 0, 1, 0]:
            Bs.append((nz[0][0], [frozenset(nm(v) for v in f) for f in Fm]))
print(f"A : {len(clsA)} classes ; B : {len(Bs)} pieces (degre, kappa) :", [(d, len(F)) for d, F in Bs])
G4 = ['g0', 'g1', 'g2', 'g3']; E4 = [list(e) for e in itertools.combinations(G4, 2)]   # K4 = graphe complet
def K26(A_, B_):
    I = [frozenset(['L_' + v for v in a] + ['R_' + v for v in b] + [g]) for a in A_ for b in B_ for g in G4]
    II = [frozenset(['L_' + v for v in b] + ['R_' + v for v in a] + e) for b in B_ for a in A_ for e in E4]
    K = ft.maximal(I + II); return K, sorted(set(v for f in K for v in f))
rows = []
for ia, A_ in enumerate(Areps):
    for ib, (dB, B_) in enumerate(Bs):
        K, Vk = K26(A_, B_)
        a = ft.aut_order(Vk, K); b = ft.betti_z2(K)
        rows.append((ia, ib, dB, len(A_), len(B_), len(Vk), len(K), a, b))
        print(f"  A{ia} (kappa {len(A_)}) x B{ib} (deg {dB}, kappa {len(B_)}) : v={len(Vk)} kappa={len(K)} |Aut|={a} Betti={b}", flush=True)
json.dump(dict(A=[sorted(map(sorted, F)) for F in Areps], B=[(d, sorted(map(sorted, F))) for d, F in Bs], rows=rows), open(SP + 'm3_pairs.json', 'w'))
