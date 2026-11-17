"""P2 pour un FHT plat dont le groupe connu est un produit de groupes agissant sur des ensembles
de feuilles disjoints (S4_L x S4_R x S_n pour la famille K^(n)). Orbites par forme canonique
facteur par facteur. Usage : python3 p2_product.py n"""
import sys, itertools, math, time, json
sys.path.insert(0, '.')
from fht_tools import *
from g_family import Kn
n = int(sys.argv[1])
K, V = Kn(n)
Lv = [v for v in V if v.startswith('L_')]; Rv = [v for v in V if v.startswith('R_')]; Gv = [v for v in V if v.startswith('g')]
def side_perms(vs, pre):
    return [{v: pre + actside(g, v[len(pre):]) for v in vs} for g in S4]
factors = [(set(Lv), side_perms(Lv, 'L_')), (set(Rv), side_perms(Rv, 'R_')),
           (set(Gv), [{('g%d' % i): ('g%d' % p[i]) for i in range(n)} for p in itertools.permutations(range(n))])]
def canon(S, facs):
    out = []
    for (Ls, H) in facs:
        part = [x for x in S if x in Ls]
        out.append(min(tuple(sorted(h[x] for x in part)) for h in H) if part else ())
    return tuple(out)
def stab(S, facs):
    return [(Ls, [h for h in H if {h[x] for x in S if x in Ls} == {x for x in S if x in Ls}]) for (Ls, H) in facs]
Kset = set(K)
t0 = time.time()
b0 = betti_z2(K); a0 = aut_order(V, K); S0, k0, v0 = math.log(a0), len(K), len(V)
print(f"K^({n}) : v={v0}, kappa={k0}, |Aut|={a0}, Betti {b0}", flush=True)
def reps(items, facs):
    seen, out = set(), []
    for it in items:
        c = canon(it if isinstance(it, (set, frozenset, tuple, list)) else [it], facs)
        if c in seen: continue
        seen.add(c); out.append(it)
    return out
edge_reps = reps(K, factors); leaf_reps = reps(V, factors)
pair_reps = reps([frozenset(p) for p in itertools.combinations(V, 2)], factors)
print(f"   orbites : aretes {len(edge_reps)}, feuilles {len(leaf_reps)}, paires {len(pair_reps)}", flush=True)
adj = set()
for f in K:
    for a, b in itertools.combinations(sorted(f), 2): adj.add(frozenset([a, b]))
results = []
def evaluate(kind, data, nE, nV):
    if not (is_clutter(nE) and connected_cover(nE, nV)): return
    results.append(dict(kind=kind, data=[sorted(x) if isinstance(x, frozenset) else x for x in data],
                        aut=aut_order(nV, nE), kappa=len(nE), v=len(nV), betti=betti_z2(nE)))
for e in edge_reps:
    st = stab(e, factors); rest = [f for f in K if f != e]
    evaluate('suppr', (e,), rest, V)
    if len(e) >= 3:
        for x in reps(sorted(e), st): evaluate('reduc', (e, x), rest + [e - {x}], V)
        for p in reps([frozenset(p) for p in itertools.combinations(sorted(e), 2)], st):
            x, y = sorted(p); evaluate('scission', (e, x, y), rest + [e - {x}, e - {y}], V)
    for x in reps([v for v in V if v not in e], st):
        ne = e | {x}; absorbed = [f for f in K if f <= ne]
        if len(absorbed) > 2: continue
        evaluate('agrand', (e, x), [f for f in K if f not in absorbed] + [ne], V)
for p in pair_reps:
    if p in adj: continue
    evaluate('ajout', (p,), K + [p], V)
for x in leaf_reps:
    t = x + "'"; evaluate('jumeau', (x,), [f | {t} if x in f else f for f in K], V + [t])
dom = [r for r in results if math.log(r['aut']) >= S0 - 1e-12 and r['kappa'] <= k0 and r['v'] <= v0
       and (math.log(r['aut']) > S0 + 1e-12 or r['kappa'] < k0 or r['v'] < v0)]
print(f"   {len(results)} classes de voisins ({time.time()-t0:.0f} s) ; dominants : {len(dom)} ; dont Betti conserves : {sum(1 for r in dom if r['betti'] == b0)}")
for r in dom[:20]: print("      ", r['kind'], r['data'], r['aut'], r['kappa'], r['v'], r['betti'])
gain = [r for r in results if r['aut'] > a0]
print("   voisins a |Aut| plus grand :", [(r['kind'], r['data'], r['aut'], r['v']) for r in gain][:12])
json.dump(dict(n=n, v=v0, kappa=k0, aut=a0, betti=b0, voisins=results), open(f'p2_K{n}.json', 'w'))
