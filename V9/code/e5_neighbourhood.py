"""V36 / E5 : voisinage morphique (mouvements d'orbite : suppr, reduc, agrand <= 2 absorptions) d'un representant,
classe par premier niveau de la chaine qui est brise :
  'N2' Betti change ; 'N2c' Betti egaux mais cone d'homologie non acyclique ; 'N3' un cone K_H non acyclique ;
  'N0' symetrie accidentelle (|Aut| != |G|) ; sinon 'classe' (reste dans la classe N3 rigide), avec dvol = variation de volume.
Usage : python3 e5_neighbourhood.py fichier.json|K26"""
import sys, itertools, time, json, collections
sys.path.insert(0, '.')
from fht_tools import *
from chain_inv import subgroup_list, factors_of, orbits
src = sys.argv[1]
K0, V = build_K3('V25_K3_pieces.json', with_pendants=False)
if src == 'K26': cur = K0
else:
    d = json.load(open(src)); cur = [frozenset(f) for f in (d['facettes'] if isinstance(d, dict) else d)]
grp = group_perms_K3(V); SL = subgroup_list(factors_of(V)); ORB = [orbits(V, g) for _, g in SL]
def img(g, x): return frozenset(g[v] for v in x) if isinstance(x, frozenset) else g[x]
def vol(F):
    S = set()
    for f in F:
        f = tuple(sorted(f))
        for k in range(1, len(f) + 1): S.update(itertools.combinations(f, k))
    return len(S)
def KH(F, O): return maximal([s for s in (frozenset(i for i, o in enumerate(O) if o <= f) for f in F) if s])
b0 = betti_z2(cur); v0 = vol(cur)
reps, s = [], set()
for e in cur:
    if e in s: continue
    o = {img(g, e) for g in grp}; s |= o; reps.append(e)
moves = []
for e in reps:
    moves.append(('suppr', (e,)))
    if len(e) >= 3: moves += [('reduc', (e, x)) for x in sorted(e)]
    moves += [('agrand', (e, x)) for x in V if x not in e]
seen = set(); tally = collections.Counter(); inclass = []
t0 = time.time()
for kind, data in moves:
    o = frozenset(tuple(img(g, x) for x in data) for g in grp)
    if (kind, o) in seen: continue
    seen.add((kind, o))
    if kind == 'suppr': nE = [f for f in cur if f not in {d[0] for d in o}]; sm = True
    elif kind == 'reduc': nE = [f for f in cur if f not in {d[0] for d in o}] + [d[0] - {d[1]} for d in o]; sm = True
    else:
        new = {d[0] | {d[1]} for d in o}
        if any(sum(1 for f in cur if f <= ne) > 2 for ne in new): tally['agrand non elementaire'] += 1; continue
        nE = list(set(cur) | new); sm = False
    nE = maximal(nE)
    if not connected_cover(nE, V): tally['non connexe'] += 1; continue
    if betti_z2(nE) != b0: tally['N2'] += 1; continue
    small, big = (nE, cur) if sm else (cur, nE)
    if betti_z2(list(big) + [f | {'*c*'} for f in small]): tally['N2c'] += 1; continue
    if any(betti_z2(list(KH(big, O)) + [f | {'*c*'} for f in KH(small, O)]) for O in ORB): tally['N3'] += 1; continue
    a = aut_order(V, nE)
    if a != 13824: tally['N0 (Aut %d)' % a] += 1; continue
    tally['classe'] += 1; inclass.append((kind, vol(nE) - v0, len(nE) - len(cur)))
print(f"{src} : kappa {len(cur)}, volume {v0} ; {len(seen)} classes de mouvements ({time.time()-t0:.0f} s)")
for k, n in sorted(tally.items()): print(f"   {k:28s} {n}")
print("   dans la classe : variations de volume", sorted(collections.Counter(d for _, d, _ in inclass).items())[:15])
print("   dans la classe, volume plus petit :", sum(1 for _, d, _ in inclass if d < 0))
