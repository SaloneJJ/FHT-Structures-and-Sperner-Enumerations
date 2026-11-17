"""V37 / A3 : descente de compacite dans la classe rigide N3 : on minimise lexicographiquement
(rayon, diametre, excentricite moyenne) du graphe des aretes LG, puis du graphe d'incidence IG, puis le volume.
Contraintes : connexe, Betti, cone d'homologie, |Aut| = 13824, cones K_H (264). Usage : python3 a3_compact.py depart"""
import sys, itertools, time, json, collections
from fractions import Fraction
sys.path.insert(0, '.')
from fht_tools import *
from chain_inv import subgroup_list, factors_of, orbits
src = sys.argv[1]
K0, V = build_K3('V25_K3_pieces.json', with_pendants=False)
cur = K0 if src == 'K26' else [frozenset(f) for f in json.load(open(src))]
grp = group_perms_K3(V); ORB = [orbits(V, g) for _, g in subgroup_list(factors_of(V))]
def img(g, x): return frozenset(g[v] for v in x) if isinstance(x, frozenset) else g[x]
def vol(F):
    S = set()
    for f in F:
        f = tuple(sorted(f))
        for k in range(1, len(f) + 1): S.update(itertools.combinations(f, k))
    return len(S)
def ecc(adj):
    n = len(adj); e = []
    for s in range(n):
        dist = [-1] * n; dist[s] = 0; q = [s]
        for u in q:
            for w in adj[u]:
                if dist[w] < 0: dist[w] = dist[u] + 1; q.append(w)
        e.append(max(dist))
    return (min(e), max(e), Fraction(sum(e), n))
def key(F):
    vi = {v: i for i, v in enumerate(V)}; E = [frozenset(vi[v] for v in f) for f in F]; n, m = len(V), len(E)
    inc = collections.defaultdict(list)
    for j, e in enumerate(E):
        for a in e: inc[a].append(j)
    lg = [set() for _ in range(m)]
    for a, js in inc.items():
        for j in js: lg[j].update(js)
    for j in range(m): lg[j].discard(j)
    ig = [set() for _ in range(n + m)]
    for j, e in enumerate(E):
        for a in e: ig[a].add(n + j); ig[n + j].add(a)
    return ecc([list(s) for s in lg]) + ecc([list(s) for s in ig]) + (vol(F),)
def KH(F, O): return maximal([s for s in (frozenset(i for i, o in enumerate(O) if o <= f) for f in F) if s])
def acyc(small, big): return not betti_z2(list(big) + [f | {'*c*'} for f in small])
b0 = betti_z2(cur); kcur = key(cur)
tag = src.split('/')[-1].replace('.json', '')
log = open(f'a3_{tag}.log', 'w')
def P(*a): print(*a, flush=True); print(*a, file=log, flush=True)
fmt = lambda k: f"LG({k[0]},{k[1]},{float(k[2]):.4f}) IG({k[3]},{k[4]},{float(k[5]):.4f}) vol {k[6]}"
P(f"depart {src} : {fmt(kcur)}")
step = 0
while True:
    step += 1; t0 = time.time(); seen = set(); cands = []
    reps, s = [], set()
    for e in cur:
        if e in s: continue
        o = {img(g, e) for g in grp}; s |= o; reps.append(e)
    moves = []
    for e in reps:
        moves.append(('suppr', (e,)))
        if len(e) >= 3: moves += [('reduc', (e, x)) for x in sorted(e)]
        moves += [('agrand', (e, x)) for x in V if x not in e]
    for kind, data in moves:
        o = frozenset(tuple(img(g, x) for x in data) for g in grp)
        if (kind, o) in seen: continue
        seen.add((kind, o))
        if kind == 'suppr': nE = [f for f in cur if f not in {d[0] for d in o}]; sm = True
        elif kind == 'reduc': nE = [f for f in cur if f not in {d[0] for d in o}] + [d[0] - {d[1]} for d in o]; sm = True
        else:
            new = {d[0] | {d[1]} for d in o}
            if any(sum(1 for f in cur if f <= ne) > 2 for ne in new): continue
            nE = list(set(cur) | new); sm = False
        nE = maximal(nE)
        if not connected_cover(nE, V): continue
        k = key(nE)
        if k < kcur: cands.append((k, kind, data, nE, sm))
    cands.sort(key=lambda c: c[0]); chosen = None
    for k, kind, data, nE, sm in cands:
        if betti_z2(nE) != b0: continue
        small, big = (nE, cur) if sm else (cur, nE)
        if not acyc(small, big) or aut_order(V, nE) != 13824: continue
        if not all(acyc(KH(small, O), KH(big, O)) for O in ORB): continue
        chosen = (k, kind, nE); break
    if chosen is None: P(f"pas {step} : aucun voisin de la classe n'ameliore la cle ({len(cands)} candidats meilleurs hors classe ou invalides ; {time.time()-t0:.0f} s)"); break
    kcur, kind, cur = chosen
    P(f"pas {step} : {kind} -> {fmt(kcur)} ({time.time()-t0:.0f} s)")
    json.dump([sorted(f) for f in cur], open(f'a3_{tag}_final.json', 'w'))
P(f"final {src} : {fmt(kcur)}, kappa {len(cur)}")
json.dump([sorted(f) for f in cur], open(f'a3_{tag}_final.json', 'w'))
