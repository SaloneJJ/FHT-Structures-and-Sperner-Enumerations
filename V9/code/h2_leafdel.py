"""Extension (hors V28) : suppression d'une orbite de feuilles (x retire de toutes les aretes), sur K26 et K26*.
Seul moyen de faire baisser v ; on teste l'homologie (Betti Z/2) et |Aut|."""
import sys, json
sys.path.insert(0, '.')
from fht_tools import *
K, V = build_K3('V25_K3_pieces.json', with_pendants=False)
Ks = [frozenset(f) for f in json.load(open('V32_FHT_K26etoile.json'))['facettes']]
grp = group_perms_K3(V)
orbs, seen = [], set()
for x in V:
    if x in seen: continue
    o = frozenset(g[x] for g in grp); seen |= o; orbs.append(o)
for name, F in (('K26', K), ('K26*', Ks)):
    b0 = betti_z2(F)
    print(name, ": Betti", b0)
    for o in orbs:
        nE = maximal([f - o for f in F if f - o])
        nV = [v for v in V if v not in o]
        ok = connected_cover(nE, nV)
        b = betti_z2(nE) if ok else None
        print(f"   - {sorted(o)} : v {len(nV)}, kappa {len(nE)}, Betti {b}{'  <== conserve' if b == b0 else ''}")
