"""V109 - Jauge Z/2 de Wilson sur le complexe hierarchique (bibliotheque V109_FHT_complexe), poids par type
(rB, rT, rS), constants d'un niveau a l'autre. Suit les bits de lien b (et non seulement F) pour pouvoir
evaluer le chiffre de temps avec les trois cups : 'v103' (V101-V106), 'serre' (cubique tire en arriere),
'aw' (Alexander-Whitney, Freudenthal). Mesures : flux par niveau et type ; pour chaque cup et chaque niveau,
p(q=1), entropie marginale, entropie sachant le parent, entropie de bloc par site, information mutuelle.
Usage : python3 V109_FHT_jauge_cups.py rB rT rS L0 balayages depart   -- JJ SALONE, FHT 2026"""
import sys, math, json
import numpy as np
from V109_FHT_complexe import Complexe, Cups, qpairs_v103

rB, rT, rS = map(float, sys.argv[1:4]); L0 = int(sys.argv[4]); NSW = int(sys.argv[5])
START = sys.argv[6] if len(sys.argv) > 6 else 'chaud'
rng = np.random.default_rng()
K = Complexe(L0); NLEV = K.NLEV; A = K.A.tocsr()
w = dict(B=rB, T=rT, S=rS)
logr = np.array([math.log(w[t]) for t in K.typ])
At = A.T.tocsr(); conf = (At @ A).tocsr()
color = -np.ones(K.nl, dtype=np.int64)
for a in rng.permutation(K.nl):
    used = set(color[conf.indices[conf.indptr[a]:conf.indptr[a + 1]]].tolist()); c = 0
    while c in used: c += 1
    color[a] = c
classes = [np.where(color == c)[0] for c in range(color.max() + 1)]
Acs = [A[:, C].tocsc() for C in classes]; AcsT = [M.T.tocsr() for M in Acs]
b = rng.integers(0, 2, K.nl).astype(np.int64) if START == 'chaud' else np.zeros(K.nl, dtype=np.int64)
F = (A @ b) % 2
def sweep():
    global F
    for C, M, MT in zip(classes, Acs, AcsT):
        d = MT @ ((1 - 2 * F) * logr)
        acc = (d >= 0) | (rng.random(len(C)) < np.exp(np.minimum(d, 0)))
        if acc.any():
            F = (F + M[:, acc] @ np.ones(acc.sum(), dtype=np.int64)) % 2
            b[C[acc]] ^= 1
CU = Cups(K); QP = qpairs_v103(K)
import itertools
parent_of = [np.array([K.sidx(l + 1, K.par(l, x)) for x in K.sites(l)]) for l in range(NLEV)]
blocks = []
for l in range(NLEV):
    L = K.Lv[l]
    blocks.append(np.array([[K.sidx(l, tuple(2 * c + e for c, e in zip(X, d))) for d in itertools.product((0, 1), repeat=3)]
                            for X in itertools.product(range(L // 2), repeat=3)]))
CUPS = ['v103', 'serre', 'aw']
def qall(l):
    P1, P2 = QP[l]
    return dict(v103=((F[P1] * F[P2]).sum(1)) % 2, serre=CU.q_serre(b, l), aw=CU.q_aw(b, l))
for _ in range(NSW // 4): sweep()
acc = {c: [np.zeros(2) for _ in range(NLEV)] for c in CUPS}
joint = {c: [np.zeros(4) for _ in range(NLEV)] for c in CUPS}
blk = {c: [np.zeros(256) for _ in range(NLEV)] for c in CUPS}
agree = np.zeros((NLEV, 3))
fluxl = {(l, t): [0, 0] for l in range(NLEV + 1) for t in 'BTS'}
nmeas = 0
for s in range(NSW):
    sweep()
    if s % 2: continue
    nmeas += 1
    assert ((A @ b) % 2 == F).all()
    Q = [qall(l) for l in range(NLEV)]
    for l in range(NLEV):
        agree[l] += [(Q[l]['v103'] == Q[l]['serre']).mean(), (Q[l]['serre'] == Q[l]['aw']).mean(), (Q[l]['v103'] == Q[l]['aw']).mean()]
        for c in CUPS:
            q = Q[l][c]; acc[c][l] += [q.sum(), q.size]
            code = (q[blocks[l]] * (1 << np.arange(8))).sum(1); blk[c][l] += np.bincount(code, minlength=256)
            if l < NLEV - 1: joint[c][l] += np.bincount(2 * q + Q[l + 1][c][parent_of[l]], minlength=4)
    for l in range(NLEV + 1):
        ml = K.lev == l
        for t in 'BTS':
            mk = ml & (K.typ == t)
            if mk.any(): fluxl[(l, t)][0] += F[mk].sum(); fluxl[(l, t)][1] += mk.sum()
def Hn(c):
    c = c[c > 0] / c.sum(); return float(-(c * np.log2(c)).sum())
Hb = lambda p: 0.0 if p in (0, 1) else -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
out = dict(rB=rB, rT=rT, rS=rS, L0=L0, balayages=NSW, mesures=nmeas, depart=START,
           flux_par_niveau={'%d%s' % k: round(f[0] / f[1], 5) for k, f in fluxl.items() if f[1]},
           accord_v103_serre_aw=[[round(v / nmeas, 4) for v in row] for row in agree], cups={})
for c in CUPS:
    lv = []
    for l in range(NLEV):
        p = acc[c][l][0] / acc[c][l][1]
        rec = dict(niveau=l, p_q1=round(p, 4), H_marg=round(Hb(p), 4), H_bloc_par_site=round(Hn(blk[c][l]) / 8, 4))
        if l < NLEV - 1:
            j = joint[c][l]; hc = Hn(j) - Hn(np.array([j[0] + j[2], j[1] + j[3]]))
            rec['H_cond_parent'] = round(hc, 4); rec['I_parent'] = round(Hb(p) - hc, 5)
        lv.append(rec)
    out['cups'][c] = lv
print(json.dumps(out, ensure_ascii=False))
