"""H2 : mouvements appliques a une orbite entiere (V28, P0 : 'morphisme applique a une orbite') sur K26 / K3.
Chaque mouvement isole (suppr, reduc, scission, agrand, ajout, jumeau) est remplace par la reunion de ses
images sous le groupe connu, appliquees simultanement ; puis aretes maximales. Le groupe reste un sous-groupe
de Aut, donc S ne baisse jamais : le critere de Pareto nu degenere ; on rapporte les voisins qui conservent
l'homologie (puits contraint, F2)."""
import sys, itertools, math, time, json
sys.path.insert(0, '.')
from fht_tools import *
name = sys.argv[1]
K, V = build_K3('V25_K3_pieces.json', with_pendants=(name == 'K3'))
grp = group_perms_K3(V)
b0 = betti_z2(K); a0 = aut_order(V, K); k0, v0 = len(K), len(V)
print(f"{name} : v={v0}, kappa={k0}, |Aut|={a0}, Betti {b0}", flush=True)
def img(g, x): return frozenset(g[v] for v in x) if isinstance(x, frozenset) else g[x]
def orbit_of(data): return {tuple(img(g, x) for x in data) for g in grp}
seen = set(); moves = []
def add(kind, data):
    o = frozenset(orbit_of(data))
    if (kind, o) in seen: return
    seen.add((kind, o)); moves.append((kind, data, o))
reps_e, s = [], set()
for e in K:
    if e in s: continue
    o = {img(g, e) for g in grp}; s |= o; reps_e.append(e)
adj = set()
for f in K:
    for a, b in itertools.combinations(sorted(f), 2): adj.add(frozenset([a, b]))
for e in reps_e:
    add('suppr', (e,))
    if len(e) >= 3:
        for x in sorted(e): add('reduc', (e, x))
    for x in V:
        if x not in e: add('agrand', (e, x))
for p in itertools.combinations(V, 2):
    p = frozenset(p)
    if p not in adj: add('ajout', (p,))
for x in V: add('jumeau', (x,))
print("   classes de mouvements par orbite :", len(moves), flush=True)
t0 = time.time(); res = []
Kset = set(K)
for kind, data, o in moves:
    if kind == 'suppr':
        rem = {d[0] for d in o}; nE = [f for f in K if f not in rem]; nV = V
    elif kind == 'reduc':
        rem = {d[0] for d in o}; nE = [f for f in K if f not in rem] + [d[0] - {d[1]} for d in o]; nV = V
    elif kind == 'agrand':
        new = {d[0] | {d[1]} for d in o}
        if any(sum(1 for f in K if f <= ne) > 2 for ne in new): continue   # agrand elementaire : absorbe au plus 2 aretes (V28)
        nE = list(set(K) | new); nV = V
    elif kind == 'ajout':
        nE = K + [d[0] for d in o]; nV = V
    else:
        X = {d[0] for d in o}; nE = [f | {x + "'" for x in f & X} for f in K]; nV = V + [x + "'" for x in X]
    nE = maximal(nE)
    if not connected_cover(nE, nV): continue
    b = betti_z2(nE)
    if b != b0: res.append(dict(kind=kind, data=[sorted(x) if isinstance(x, frozenset) else x for x in data], betti=b, kappa=len(nE), v=len(nV))); continue
    a = aut_order(nV, nE)
    res.append(dict(kind=kind, data=[sorted(x) if isinstance(x, frozenset) else x for x in data], betti=b, kappa=len(nE), v=len(nV), aut=a))
kept = [r for r in res if r['betti'] == b0]
dom = [r for r in kept if r['aut'] >= a0 and r['kappa'] <= k0 and r['v'] <= v0 and (r['aut'] > a0 or r['kappa'] < k0 or r['v'] < v0)]
print(f"   voisins valides : {len(res)} ; a Betti conserves : {len(kept)} ; dominants a Betti conserves : {len(dom)} ({time.time()-t0:.0f} s)")
nodeg = [r for r in res if r['kappa'] < k0 and r['v'] <= v0]
print(f"   (Pareto nu, homologie libre : voisins a kappa < kappa0 et v <= v0 : {len(nodeg)} ; S >= S0 automatiquement)")
for r in dom[:20]: print("      ", r)
print("   Betti conserves, par type :", {k: sum(1 for r in kept if r['kind'] == k) for k in ['suppr', 'reduc', 'agrand', 'ajout', 'jumeau']})
print("   Betti conserves avec kappa < kappa0 :", [(r['kind'], r['kappa'], r['aut']) for r in kept if r['kappa'] < k0][:20])
json.dump(res, open(f'h2_orbit_{name}.json', 'w'))
