"""V37 / A2 : invariants combinatoires de l'article Invariants (Part II) sur les objets du programme :
triplets (rayon, diametre, excentricite moyenne) du graphe primal PG, du graphe des aretes LG et du graphe
d'incidence IG ; profils de degres des feuilles, des tailles d'aretes ; |E(PG)|, |E(LG)|, |E(IG)|."""
import sys, json, collections, itertools
from fractions import Fraction
sys.path.insert(0, '.')
from chain_inv import OBJ, obj_json
OBJ['K26vol'] = obj_json('e4_K26_final.json'); OBJ['K26*vol'] = obj_json('e4_K26etoile_final.json')
def bfs_all(adj):
    n = len(adj); ecc = []
    for s in range(n):
        dist = [-1] * n; dist[s] = 0; q = [s]
        for u in q:
            for w in adj[u]:
                if dist[w] < 0: dist[w] = dist[u] + 1; q.append(w)
        ecc.append(max(dist))
    return min(ecc), max(ecc), Fraction(sum(ecc), n)
res = {}
for nm, (K, V) in OBJ.items():
    vi = {v: i for i, v in enumerate(V)}; E = [frozenset(vi[v] for v in f) for f in K]
    n, m = len(V), len(E)
    pg = [set() for _ in range(n)]
    for e in E:
        for a in e: pg[a] |= e - {a}
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
    r = dict(PG=bfs_all([list(s) for s in pg]), LG=bfs_all([list(s) for s in lg]), IG=bfs_all([list(s) for s in ig]),
             E_PG=sum(len(s) for s in pg) // 2, E_LG=sum(len(s) for s in lg) // 2, E_IG=sum(len(e) for e in E),
             degres=sorted(collections.Counter(len(inc[a]) for a in range(n)).items()),
             tailles=sorted(collections.Counter(len(e) for e in E).items()))
    res[nm] = r
    f = lambda t: f"({t[0]},{t[1]},{float(t[2]):.3f})"
    print(f"{nm:8s} PG{f(r['PG'])} LG{f(r['LG'])} IG{f(r['IG'])} |E(PG)|={r['E_PG']} |E(LG)|={r['E_LG']} |E(IG)|={r['E_IG']} "
          f"tailles={r['tailles']} degres(min,max)=({r['degres'][0][0]},{r['degres'][-1][0]})", flush=True)
json.dump({k: {kk: (str(vv) if kk in ('PG', 'LG', 'IG') else vv) for kk, vv in v.items()} for k, v in res.items()}, open('a2_invariants.json', 'w'))
