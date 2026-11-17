"""V36 / E4 : descente gloutonne par mouvements d'orbite qui minimise le VOLUME (nombre de simplexes)
en restant dans la classe N3 : Betti Z/2 conserves, puis pour le mouvement retenu, cones acycliques
pour l'homologie et pour les 264 K_H, et |Aut| = |G| (rigidite, niveau N0). Usage : python3 e4_volume_descent.py depart (K26 | K26etoile)"""
import sys, itertools, time, json
sys.path.insert(0, '.')
from fht_tools import *
from chain_inv import subgroup_list, factors_of, orbits
start = sys.argv[1]
K0, V = build_K3('V25_K3_pieces.json', with_pendants=False)
K = K0 if start == 'K26' else [frozenset(f) for f in json.load(open('V32_FHT_K26etoile.json'))['facettes']]
grp = group_perms_K3(V)
SL = subgroup_list(factors_of(V)); ORB = [orbits(V, g) for _, g in SL]
def img(g, x): return frozenset(g[v] for v in x) if isinstance(x, frozenset) else g[x]
def vol(F):
    S = set()
    for f in F:
        f = tuple(sorted(f))
        for k in range(1, len(f) + 1): S.update(itertools.combinations(f, k))
    return len(S)
def KH(F, O): return maximal([s for s in (frozenset(i for i, o in enumerate(O) if o <= f) for f in F) if s])
def cone_ok(small, big):
    if betti_z2(list(big) + [f | {'*c*'} for f in small]): return False
    for O in ORB:
        s, b = KH(small, O), KH(big, O)
        if betti_z2(list(b) + [f | {'*c*'} for f in s]): return False
    return True
b0 = betti_z2(K); log = open(f'e4_{start}.log', 'w')
def P(*a): print(*a, flush=True); print(*a, file=log, flush=True)
cur = K; P(f"depart {start} : kappa {len(cur)}, volume {vol(cur)}")
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
    v0 = vol(cur)
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
        vn = vol(nE)
        if vn >= v0: continue
        cands.append((vn, kind, data, nE, sm))
    cands.sort(key=lambda c: c[0])
    chosen = None
    for vn, kind, data, nE, sm in cands:
        if not connected_cover(nE, V) or betti_z2(nE) != b0: continue
        small, big = (nE, cur) if sm else (cur, nE)
        if not cone_ok(small, big): continue
        if aut_order(V, nE) != 13824: continue      # N0 : l'algebre se lit sur Aut ; symetrie accidentelle interdite (rigidite)
        chosen = (vn, kind, data, nE); break
    if chosen is None: P(f"pas {step} : aucun mouvement ne baisse le volume dans la classe N3 ({len(cands)} candidats ; {time.time()-t0:.0f} s)"); break
    vn, kind, data, nE = chosen
    P(f"pas {step} : {kind} -> volume {vn}, kappa {len(nE)}, |Aut| {aut_order(V, nE)} ({time.time()-t0:.0f} s)")
    cur = nE
    json.dump([sorted(f) for f in cur], open(f'e4_{start}_final.json', 'w'))
P(f"final {start} : volume {vol(cur)}, kappa {len(cur)}, |Aut| {aut_order(V, cur)}, facettes {sorted(set(len(f) for f in cur))}")
