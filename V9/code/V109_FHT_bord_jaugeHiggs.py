"""V109 - Le niveau de la tuile seul (niveaux superieurs plats) : modele de jauge-Higgs Z/2 en dimension 3
(liens horizontaux = jauge, liens verticaux = champ de Higgs). Test de l'auto-dualite E <-> B :
au point auto-dual r = sqrt(2) - 1, l'identite de dualite impose <F_B> + <F_E> = r/(1+r) = 1 - 1/sqrt(2).
Usage : python3 V109_FHT_bord_jaugeHiggs.py rB rE L balayages   -- JJ SALONE, FHT 2026"""
import sys, math, json, itertools
import numpy as np, scipy.sparse as sp
rB, rE = float(sys.argv[1]), float(sys.argv[2]); L = int(sys.argv[3]); NSW = int(sys.argv[4])
rng = np.random.default_rng()
N = L ** 3
def s(x): return ((x[0] % L) * L + (x[1] % L)) * L + (x[2] % L)
def sh(x, k): y = list(x); y[k] += 1; return tuple(y)
h = lambda x, k: 3 * s(x) + k
v = lambda x: 3 * N + s(x)
nl = 4 * N; rows, cols, typ = [], [], []
for x in itertools.product(range(L), repeat=3):
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        p = len(typ); typ.append('B')
        for a in (h(x, i), h(sh(x, i), j), h(sh(x, j), i), h(x, j)): rows.append(p); cols.append(a)
    for k in range(3):
        p = len(typ); typ.append('E')
        for a in (h(x, k), v(x), v(sh(x, k))): rows.append(p); cols.append(a)
A = sp.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=(len(typ), nl))
typ = np.array(typ); logr = np.where(typ == 'B', math.log(rB), math.log(rE))
At = A.T.tocsr(); conf = (At @ A).tocsr()
color = -np.ones(nl, dtype=np.int64)
for a in rng.permutation(nl):
    used = set(color[conf.indices[conf.indptr[a]:conf.indptr[a + 1]]].tolist()); c = 0
    while c in used: c += 1
    color[a] = c
classes = [np.where(color == c)[0] for c in range(color.max() + 1)]
Ms = [A[:, C].tocsc() for C in classes]; MTs = [M.T.tocsr() for M in Ms]
F = np.zeros(len(typ), dtype=np.int64)
def sweep():
    global F
    for C, M, MT in zip(classes, Ms, MTs):
        d = MT @ ((1 - 2 * F) * logr)
        acc = (d >= 0) | (rng.random(len(C)) < np.exp(np.minimum(d, 0)))
        if acc.any(): F = (F + M[:, acc] @ np.ones(acc.sum(), dtype=np.int64)) % 2
for _ in range(NSW // 4): sweep()
xs = []
for it in range(NSW):
    sweep()
    if it % 2 == 0: xs.append((F[typ == 'B'].mean(), F[typ == 'E'].mean()))
xs = np.array(xs); xB, xE = xs.mean(0); eB, eE = xs.std(0) / math.sqrt(len(xs) / 10)
print(json.dumps(dict(rB=rB, rE=rE, L=L, balayages=NSW, xB=round(xB, 5), xE=round(xE, 5), somme=round(xB + xE, 5),
                      erreur=round(math.hypot(eB, eE), 5),
                      prediction_autoduale=round(rB / (1 + rB), 5) if abs(rE - (1 - rB) / (1 + rB)) < 1e-9 else None)))
