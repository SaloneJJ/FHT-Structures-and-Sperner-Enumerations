"""H2 / P3+ : descente gloutonne par mouvements d'orbite sur K26 sous contrainte d'homologie.
A chaque pas : parmi les mouvements d'orbite (suppr, reduc, agrand, scission) qui conservent les Betti Z/2
et le caractere de Lefschetz virtuel sous S4^3 (125 classes), et gardent |Aut| >= 13824, on prend celui qui
minimise kappa. Arret quand aucun ne diminue kappa. Journal pas a pas dans h2_descent.log / .json."""
import sys, itertools, math, time, json
sys.path.insert(0, '.')
from fht_tools import *
from simplicial import perm_sign
K, V = build_K3('V25_K3_pieces.json', with_pendants=False)
grp = group_perms_K3(V)
# representants des 125 classes de conjugaison de S4^3
def cyc(g):
    seen, t = set(), []
    for i in range(4):
        if i in seen: continue
        n_, j = 0, i
        while j not in seen: seen.add(j); j = g[j]; n_ += 1
        t.append(n_)
    return tuple(sorted(t))
cl = {}
for g in S4: cl.setdefault(cyc(g), g)
reps = list(cl.values())
cmaps = [{v: act3(a, b, c, v) for v in V} for a in reps for b in reps for c in reps]
def lefschetz_vec(F):
    faces = set()
    for f in F:
        f = tuple(sorted(f))
        for k in range(1, len(f) + 1):
            for s in itertools.combinations(f, k): faces.add(s)
    out = []
    for m in cmaps:
        L = -1
        for s in faces:
            im = [m[v] for v in s]
            if set(im) == set(s):
                pos = {v: i for i, v in enumerate(s)}
                L += (-1) ** (len(s) - 1) * perm_sign([pos[x] for x in im])
        out.append(L)
    return tuple(out)
def img(g, x): return frozenset(g[v] for v in x) if isinstance(x, frozenset) else g[x]
b0 = betti_z2(K); L0 = lefschetz_vec(K)
log = open('h2_descent.log', 'w')
def P(*a):
    print(*a, flush=True); print(*a, file=log, flush=True)
P(f"depart K26 : kappa {len(K)}, Betti {b0}")
cur = K; step = 0; hist = []
while True:
    step += 1; t0 = time.time()
    seen, cand = set(), []
    reps_e, s = [], set()
    for e in cur:
        if e in s: continue
        o = {img(g, e) for g in grp}; s |= o; reps_e.append(e)
    moves = []
    for e in reps_e:
        moves.append(('suppr', (e,)))
        if len(e) >= 3:
            for x in sorted(e): moves.append(('reduc', (e, x)))
        for x in V:
            if x not in e: moves.append(('agrand', (e, x)))
    best = None
    for kind, data in moves:
        o = frozenset(tuple(img(g, x) for x in data) for g in grp)
        if (kind, o) in seen: continue
        seen.add((kind, o))
        if kind == 'suppr':
            rem = {d[0] for d in o}; nE = [f for f in cur if f not in rem]
        elif kind == 'reduc':
            rem = {d[0] for d in o}; nE = [f for f in cur if f not in rem] + [d[0] - {d[1]} for d in o]
        else:
            new = {d[0] | {d[1]} for d in o}
            if any(sum(1 for f in cur if f <= ne) > 2 for ne in new): continue   # agrand elementaire (V28)
            nE = list(set(cur) | new)
        nE = maximal(nE)
        if len(nE) >= len(cur) or (best and len(nE) >= len(best[2])): continue
        if not connected_cover(nE, V): continue
        if betti_z2(nE) != b0: continue
        best = (kind, data, nE)
    if best is None: P(f"pas {step} : aucun mouvement ne diminue kappa ; arret ({time.time()-t0:.0f} s)"); break
    kind, data, nE = best
    Lv = lefschetz_vec(nE); a = aut_order(V, nE)
    P(f"pas {step} : {kind} {[sorted(x) if isinstance(x, frozenset) else x for x in data]} -> kappa {len(nE)}, |Aut| {a}, Lefschetz = K26 : {Lv == L0} ({time.time()-t0:.0f} s)")
    if Lv != L0 or a < 13824:
        P("   rejete (caractere ou symetrie) ; arret prudent"); break
    cur = nE; hist.append(dict(step=step, kind=kind, kappa=len(nE), aut=a))
    json.dump(dict(hist=hist, facets=[sorted(f) for f in cur]), open('h2_descent.json', 'w'))
P("final : kappa", len(cur), "; |Aut|", aut_order(V, cur), "; dimensions des facettes", sorted(set(len(f) for f in cur)))
