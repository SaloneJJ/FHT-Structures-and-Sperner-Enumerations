"""V148 (c) : structures de d-cube natives dans K26 (V148 §1, criteres ecrits avant calcul).
Une structure de d-cube sur S (|S| = 2^d) = d bissections S = A u B en deux simplexes de K de taille 2^(d-1),
dont les d vecteurs de signes sont deux a deux distincts. Orbites sous Aut(K26) = S4^3 (13 824).
Version 2 (sortie compacte ; la version 1 imprimait toutes les orbites de d = 2).
Usage : python3 V148_FHT_c_structures_cube.py [d ...]      JJ SALONE, FHT 2026"""
import itertools, sys, time
import numpy as np
from V86_FHT_k26 import build_K
from fht_tools import group_perms_K3

t00 = time.time()
K, I, II = build_K()
V = sorted(set().union(*K)); n = len(V); idx = {v: i for i, v in enumerate(V)}
Kmask = [sum(1 << idx[v] for v in f) for f in K]; Kset = set(Kmask)
print(f'K26 : {n} sommets, {len(K)} facettes', flush=True)
perms = [[idx[p[v]] for v in V] for p in group_perms_K3(V)]
def img_mask(m, p):
    r = 0
    for i in range(n):
        if m >> i & 1: r |= 1 << p[i]
    return r
perms = [p for p in perms if set(img_mask(m, p) for m in Kmask) == Kset]
print(f'|Aut| (S4^3 qui preserve K) = {len(perms)}', flush=True)
G = np.array(perms, dtype=np.int64)
def pc(x): return bin(x).count('1')

def simplices(k):
    s = set()
    for f in Kmask:
        l = [i for i in range(n) if f >> i & 1]
        for c in itertools.combinations(l, k):
            s.add(sum(1 << i for i in c))
    return np.array(sorted(s), dtype=np.int64)

def bis_pairs(M):
    ia, ib = [], []
    for i in range(len(M)):
        sel = np.nonzero(((M[i + 1:] & M[i]) == 0))[0] + i + 1
        if len(sel):
            ia.append(np.full(len(sel), i, dtype=np.int32)); ib.append(sel.astype(np.int32))
    ia = np.concatenate(ia); ib = np.concatenate(ib)
    return M[ia] | M[ib], ia, ib

def is_cube(S, Alist):
    """les 2^d atomes (intersections des A_i ou de leurs complements dans S) ont tous exactement un element"""
    atoms = [S]
    for A in Alist:
        atoms = [a & A for a in atoms] + [a & ~A for a in atoms]
    return all(pc(a) == 1 for a in atoms)

def orbit_image(S):
    res = np.zeros(len(G), dtype=np.int64)
    for i in range(n):
        if S >> i & 1: res |= (1 << G[:, i])
    return res

dims = [int(a) for a in sys.argv[1:]] or [2, 3, 4]
XOCT = {f'{a}_x{c}' for a in 'LR' for c in range(4)}
for d in dims:
    k = 2 ** (d - 1)
    print(f'\n=== d = {d} : bissections en simplexes de taille {k} ===', flush=True)
    M = simplices(k)
    print(f'  simplexes de taille {k} : {len(M)}', flush=True)
    SA, ia, ib = bis_pairs(M)
    print(f'  paires disjointes : {len(SA)}  ({time.time()-t00:.0f} s)', flush=True)
    order = np.argsort(SA, kind='stable'); SA = SA[order]; ia = ia[order]
    uS, start, cnt = np.unique(SA, return_index=True, return_counts=True)
    cand = np.nonzero(cnt >= d)[0]
    print(f'  ensembles S avec >= {d} bissections : {len(cand)}', flush=True)
    structs = {}
    for ci, c in enumerate(cand):
        S = int(uS[c]); As = [int(M[ia[j]]) for j in range(start[c], start[c] + cnt[c])]
        found = [comb for comb in itertools.combinations(As, d) if is_cube(S, comb)]
        if found: structs[S] = found
        if ci % 50000 == 0 and ci: print(f'    {ci}/{len(cand)} ({time.time()-t00:.0f} s)', flush=True)
    tot = sum(len(v) for v in structs.values())
    print(f'  ensembles S portant une structure de {d}-cube : {len(structs)} ; structures au total : {tot}', flush=True)
    if not structs: continue
    seen = set(); rows = []
    for S in sorted(structs):
        if S in seen: continue
        res = orbit_image(S)
        orbit = set(int(x) for x in np.unique(res)); seen |= orbit
        stabidx = np.nonzero(res == S)[0]
        stab = [perms[j] for j in stabidx]
        sts = [frozenset(A for A in comb) for comb in structs[S]]      # structure = ensemble des A_i (bissection = {A, S^A})
        def canon(A): return min(A, S ^ A)
        stset = set(frozenset(canon(A) for A in st) for st in sts)
        left = set(stset); sorbs = []
        while left:
            s0 = next(iter(left)); st_sz = 0; orb = set()
            for p in stab:
                im = frozenset(canon(img_mask(A, p)) for A in s0); orb.add(im)
                if im == s0: st_sz += 1
            left -= orb; sorbs.append((len(orb), st_sz))
        rows.append((S, [V[i] for i in range(n) if S >> i & 1], len(orbit), len(stab), len(stset), sorted(sorbs, key=lambda z: -z[1])))
    print(f'  orbites de S portant une structure : {len(rows)}')
    maxstab = max(s for r in rows for (_, s) in r[5])
    print(f'  plus grand stabilisateur de structure : {maxstab}')
    rows.sort(key=lambda r: (-max(s for _, s in r[5]), r[2]))
    xrow = [r for r in rows if set(r[1]) == XOCT]
    shown = rows[:12] + [r for r in xrow if r not in rows[:12]]
    for S, names, osz, stabS, nst, sorbs in shown:
        tag = '   <-- octuplet des x' if set(names) == XOCT else ''
        print(f'  S = {" ".join(names)}{tag}\n     |orbite de S| = {osz} ; |Stab(S)| = {stabS} ; structures sur S : {nst} ; (taille d\'orbite, |Stab|) : {sorbs}')
    from collections import Counter
    print('  histogramme (|Stab| max de la structure) -> nombre d\'orbites de S :',
          dict(sorted(Counter(max(s for _, s in r[5]) for r in rows).items(), reverse=True)))
    # criteres K2 : combien d'orbites ont un stabilisateur de structure >= celui de l'octuplet x ?
    if xrow:
        sx = max(s for _, s in xrow[0][5])
        print(f'  stabilisateur de structure (x) = {sx} ; orbites de S avec |Stab| de structure >= {sx} : {sum(1 for r in rows if max(s for _, s in r[5]) >= sx)} ; strictement > : {sum(1 for r in rows if max(s for _, s in r[5]) > sx)}')
print(f'\ntemps total {time.time()-t00:.0f} s')
