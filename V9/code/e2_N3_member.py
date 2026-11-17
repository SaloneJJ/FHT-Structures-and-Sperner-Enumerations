"""V34 / E2 : K26 et K26* sont-ils dans la meme classe N3 ? Pour chaque pas de la descente (h2_descent.log)
et chaque sous-groupe H de la liste (264), on teste l'acyclicite (Z/2) du cone de l'inclusion K'_H c K_H."""
import sys, re, ast
sys.path.insert(0, '.')
from fht_tools import *
from chain_inv import subgroup_list, factors_of, orbits
K, V = build_K3('V25_K3_pieces.json', with_pendants=False)
grp = group_perms_K3(V)
facs = factors_of(V); SL = subgroup_list(facs)
def img(g, x): return frozenset(g[v] for v in x) if isinstance(x, frozenset) else g[x]
def KH(F, O):
    out = [frozenset(i for i, o in enumerate(O) if o <= f) for f in F]
    return maximal([s for s in out if s])
cur = K; bad = 0; tested = 0
for line in open('h2_descent.log'):
    m = re.match(r"pas (\d+) : (\w+) (\[.*\]) -> kappa (\d+)", line)
    if not m: continue
    kind, data = m.group(2), tuple(frozenset(x) if isinstance(x, list) else x for x in ast.literal_eval(m.group(3)))
    o = {tuple(img(g, x) for x in data) for g in grp}
    if kind == 'suppr': nE = maximal([f for f in cur if f not in {d[0] for d in o}]); small, big = nE, cur
    elif kind == 'reduc': nE = maximal([f for f in cur if f not in {d[0] for d in o}] + [d[0] - {d[1]} for d in o]); small, big = nE, cur
    else: nE = maximal(list(set(cur) | {d[0] | {d[1]} for d in o})); small, big = cur, nE
    nbad = 0
    for lab, gens in SL:
        O = orbits(V, gens)
        s, b = KH(small, O), KH(big, O)
        cone = betti_z2(list(b) + [f | {'*c*'} for f in s])
        tested += 1
        if cone: nbad += 1
    bad += nbad
    print(f"pas {m.group(1)} ({kind}) : sous-groupes ou le cone n'est pas acyclique : {nbad}/{len(SL)}", flush=True)
    cur = nE
print("total des echecs :", bad, "sur", tested)
