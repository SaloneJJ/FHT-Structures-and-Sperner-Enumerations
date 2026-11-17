"""Rejoue la descente h2_descent.log et teste a chaque pas si l'inclusion (K' c K pour suppr/reduc,
K c K' pour agrand) induit un isomorphisme en homologie Z/2 : H~(K u cone(K')) = 0 (cone de l'application)."""
import sys, re, ast
sys.path.insert(0, '.')
from fht_tools import *
K, V = build_K3('V25_K3_pieces.json', with_pendants=False)
grp = group_perms_K3(V)
def img(g, x): return frozenset(g[v] for v in x) if isinstance(x, frozenset) else g[x]
def mapping_cone_betti(small, big):
    return betti_z2(list(big) + [frozenset(f) | {'*c*'} for f in small])
cur = K
for line in open('h2_descent.log'):
    m = re.match(r"pas (\d+) : (\w+) (\[.*\]) -> kappa (\d+)", line)
    if not m: continue
    kind, data = m.group(2), ast.literal_eval(m.group(3))
    data = tuple(frozenset(x) if isinstance(x, list) else x for x in data)
    o = {tuple(img(g, x) for x in data) for g in grp}
    if kind == 'suppr':
        rem = {d[0] for d in o}; nE = maximal([f for f in cur if f not in rem]); small, big = nE, cur
    elif kind == 'reduc':
        rem = {d[0] for d in o}; nE = maximal([f for f in cur if f not in rem] + [d[0] - {d[1]} for d in o]); small, big = nE, cur
    else:
        nE = maximal(list(set(cur) | {d[0] | {d[1]} for d in o})); small, big = cur, nE
    assert len(nE) == int(m.group(4))
    print(f"pas {m.group(1)} ({kind}) : kappa {len(nE)} ; H~(cone de l'inclusion) = {mapping_cone_betti(small, big)}", flush=True)
    cur = nE
