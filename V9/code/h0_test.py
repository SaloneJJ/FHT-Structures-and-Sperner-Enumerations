"""H0.2 : lectures R0, R1, R2 sur tous les FHT a v <= 5 (axiome A5 par composant : 1, 1, 4, 23, 231)."""
import sys, collections
sys.path.insert(0, '.')
from fht_hier import *
from fht_enum_small import fhts, canon
from fht_tools import aut_order
stat = collections.Counter(); ex = []
for n in range(2, 6):
    reps = {}
    for E in fhts(tuple(range(n))): reps.setdefault(canon(n, E), E)
    for key, E in sorted(reps.items()):
        t = decompose(set(range(n)), E)
        b0, b1, b2, bh = betti_z2(R0(t)), betti_z2(R1(t)[:]), betti_z2(R2(t)[0]), betti_z2(contracted_host(t))
        flat = depth(t) == 1
        stat[(n, 'plats' if flat else 'hier')] += 1
        if flat: stat['plat : R1 = R0'] += (b1 == b0); stat['plats total'] += 1; continue
        stat['hier total'] += 1
        stat['hier : R0 = hote contracte'] += (b0 == bh)
        stat['hier : R2 = R0'] += (b2 == b0)
        if b1 != b0:
            stat['hier : R1 != R0'] += 1
            # ombre plate : R1 vu comme FHT plat, meme feuilles ; Aut(F) <= Aut(ombre)
            aF = aut_order(list(range(n)), E); aS = aut_order(list(range(n)), R1(t))
            ex.append((n, sorted(map(lambda e: ''.join(str(x + 1) for x in sorted(e)), E)), 'R0', b0, 'R1', b1, '|Aut F|', aF, '|Aut ombre|', aS))
print({k: v for k, v in stat.items()})
for e in ex: print("  ", e)
