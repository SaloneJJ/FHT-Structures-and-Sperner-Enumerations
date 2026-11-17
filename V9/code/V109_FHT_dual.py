"""V109 - Modele dual exact (MacWilliams) de la jauge Z/2 hierarchique : surfaces fermees.
Z_jauge(r) = Sum_{F in im A} Prod r_p^F_p  =  (|im A| / 2^N) Prod (1 + r_p)  Sum_{y in ker A^T} Prod s_p^y_p,
s_p = (1 - r_p) / (1 + r_p). ker A^T = im G^T (bords des cellules de dimension 3) + H2 (3 surfaces globales
du tore du niveau superieur). Variables duales : un bit par cellule de dimension 3 (+ 3 bits globaux) ;
y = G^T sigma + surfaces globales. Metropolis par classes de couleur.
Identite exacte, plaquette par plaquette :  <F_p>_r = r_p/(1+r_p) - 2 r_p / (1 - r_p^2) <y_p>_s.
Usage : python3 V109_FHT_dual.py r L0 balayages [rB rT rS]  (r : poids primal ; le dual est simule a s(r))
JJ SALONE, FHT 2026"""
import sys, math, json
import numpy as np, scipy.sparse as sp
from V109_FHT_complexe import Complexe, PAIRS

L0 = int(sys.argv[2]); NSW = int(sys.argv[3])
if len(sys.argv) > 6: rB, rT, rS = map(float, sys.argv[4:7])
else: rB = rT = rS = float(sys.argv[1])
rng = np.random.default_rng()
K = Complexe(L0); NLEV = K.NLEV
r = dict(B=rB, T=rT, S=rS); s = {t: (1 - v) / (1 + v) for t, v in r.items()}
logs = np.array([math.log(s[t]) for t in K.typ])
# generateurs : cellules 3 + surfaces globales (plans du tore superieur)
rows, cols = [], []
Gc = K.G.tocsr()
for c in range(K.n3):
    for p in Gc.indices[Gc.indptr[c]:Gc.indptr[c + 1]]: rows.append(c); cols.append(p)
ng = K.n3; Ltop = K.Lv[NLEV]
for (i, j) in PAIRS:
    k = 3 - i - j
    for x in K.sites(NLEV):
        if x[k] == 0: rows.append(ng); cols.append(K.B(NLEV, x, i, j))
    ng += 1
Gen = sp.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=(ng, K.np_))
_P = (Gen @ K.A).tocsr(); assert (_P.data % 2 == 0).all()
GT = Gen.T.tocsc()
conf = (Gen @ Gen.T).tocsr()
color = -np.ones(ng, dtype=np.int64)
for a in rng.permutation(ng):
    used = set(color[conf.indices[conf.indptr[a]:conf.indptr[a + 1]]].tolist()); c = 0
    while c in used: c += 1
    color[a] = c
classes = [np.where(color == c)[0] for c in range(color.max() + 1)]
Ms = [GT[:, C].tocsc() for C in classes]; MTs = [M.T.tocsr() for M in Ms]
y = np.zeros(K.np_, dtype=np.int64)
def sweep():
    global y
    for C, M, MT in zip(classes, Ms, MTs):
        d = MT @ ((1 - 2 * y) * logs)
        acc = (d >= 0) | (rng.random(len(C)) < np.exp(np.minimum(d, 0)))
        if acc.any(): y = (y + M[:, acc] @ np.ones(acc.sum(), dtype=np.int64)) % 2
for _ in range(NSW // 4): sweep()
fl = {(l, t): [0.0, 0] for l in range(NLEV + 1) for t in 'BTS'}
nm = 0
for it in range(NSW):
    sweep()
    if it % 2: continue
    nm += 1
    for l in range(NLEV + 1):
        ml = K.lev == l
        for t in 'BTS':
            mk = ml & (K.typ == t)
            if mk.any(): fl[(l, t)][0] += y[mk].sum(); fl[(l, t)][1] += mk.sum()
out = dict(r=dict(B=rB, T=rT, S=rS), s={t: round(v, 5) for t, v in s.items()}, L0=L0, balayages=NSW, generateurs=ng,
           classes=len(classes), surfaces_par_niveau={}, flux_predit_par_niveau={})
for (l, t), (a, n) in fl.items():
    if n == 0: continue
    yy = a / n; rp = r[t]
    out['surfaces_par_niveau']['%d%s' % (l, t)] = round(yy, 5)
    out['flux_predit_par_niveau']['%d%s' % (l, t)] = round(rp / (1 + rp) - 2 * rp / (1 - rp * rp) * yy, 5)
print(json.dumps(out, ensure_ascii=False))
