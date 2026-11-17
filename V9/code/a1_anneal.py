"""V37 / A1 : recuit simule sur le volume (nombre de simplexes) dans la classe rigide N3 de K26.
Etat : FHT plat sur les 26 feuilles ; proposition : un mouvement d'orbite tire au hasard (suppr, reduc, agrand
elementaire) ; contraintes : connexe, Betti Z/2, cone d'homologie acyclique, |Aut| = 13824, cones K_H acycliques (264 H).
Acceptation de Metropolis sur le volume. Usage : python3 a1_anneal.py depart(K26|fichier.json) graine pas T0 T1"""
import sys, itertools, time, json, random, math
sys.path.insert(0, '.')
from fht_tools import *
from chain_inv import subgroup_list, factors_of, orbits
src, seed, nsteps, T0, T1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
random.seed(seed)
K0, V = build_K3('V25_K3_pieces.json', with_pendants=False)
if src == 'K26': cur = K0
else:
    d = json.load(open(src)); cur = [frozenset(f) for f in (d['facettes'] if isinstance(d, dict) else d)]
grp = group_perms_K3(V); ORB = [orbits(V, g) for _, g in subgroup_list(factors_of(V))]
def img(g, x): return frozenset(g[v] for v in x) if isinstance(x, frozenset) else g[x]
def vol(F):
    S = set()
    for f in F:
        f = tuple(sorted(f))
        for k in range(1, len(f) + 1): S.update(itertools.combinations(f, k))
    return len(S)
def KH(F, O): return maximal([s for s in (frozenset(i for i, o in enumerate(O) if o <= f) for f in F) if s])
def acyc(small, big): return not betti_z2(list(big) + [f | {'*c*'} for f in small])
b0 = betti_z2(cur); vcur = vol(cur); best = (vcur, cur)
tag = f"{src.split('/')[-1].replace('.json','')}_s{seed}"
log = open(f'a1_{tag}.log', 'w')
def P(*a): print(*a, flush=True); print(*a, file=log, flush=True)
P(f"depart {src} graine {seed} : volume {vcur}, kappa {len(cur)}")
t0 = time.time(); acc = 0
for step in range(nsteps):
    T = T0 * (T1 / T0) ** (step / max(1, nsteps - 1))
    e = random.choice(cur); kind = random.choice(['suppr', 'reduc', 'agrand', 'agrand'])
    if kind == 'reduc':
        if len(e) < 3: continue
        data = (e, random.choice(sorted(e)))
    elif kind == 'agrand': data = (e, random.choice([v for v in V if v not in e]))
    else: data = (e,)
    o = {tuple(img(g, x) for x in data) for g in grp}
    if kind == 'suppr': nE = [f for f in cur if f not in {d[0] for d in o}]; sm = True
    elif kind == 'reduc': nE = [f for f in cur if f not in {d[0] for d in o}] + [d[0] - {d[1]} for d in o]; sm = True
    else:
        new = {d[0] | {d[1]} for d in o}
        if any(sum(1 for f in cur if f <= ne) > 2 for ne in new): continue
        nE = list(set(cur) | new); sm = False
    nE = maximal(nE)
    vn = vol(nE); dv = vn - vcur
    if dv > 0 and random.random() > math.exp(-dv / T): continue
    if not connected_cover(nE, V) or betti_z2(nE) != b0: continue
    small, big = (nE, cur) if sm else (cur, nE)
    if not acyc(small, big): continue
    if aut_order(V, nE) != 13824: continue
    if not all(acyc(KH(small, O), KH(big, O)) for O in ORB): continue
    cur, vcur = nE, vn; acc += 1
    if vcur < best[0]:
        best = (vcur, cur); json.dump([sorted(f) for f in cur], open(f'a1_{tag}_best.json', 'w'))
    if acc % 10 == 0: P(f"pas {step} T={T:.0f} acceptes {acc} volume {vcur} meilleur {best[0]} ({time.time()-t0:.0f} s)")
P(f"fin : acceptes {acc}, volume courant {vcur}, meilleur {best[0]} ({time.time()-t0:.0f} s)")
json.dump([sorted(f) for f in best[1]], open(f'a1_{tag}_best.json', 'w'))
