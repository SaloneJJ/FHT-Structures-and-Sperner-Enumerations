"""V106 - Derive de V103 : flux moyens par niveau et par type (pour la recursion empirique d'un niveau au suivant).
V103 - Version rapide (numpy / scipy) de V101-V102 : jauge Z/2 de Wilson sur le complexe hierarchique
(stella, blocs 2^3, tore), poids par type de plaquette (B magnetique, T electrique triangle, S electrique carre)
et par niveau ; chiffre de temps q = F u F (cup cubique mod 2). Mise a jour par classes de couleur
(liens sans plaquette commune), Metropolis exact classe par classe.
Usage : python3 V103_FHT_jauge_rapide.py r0 modele L0 balayages depart(froid|chaud) [a]   -- JJ SALONE, FHT 2026"""
import sys, math, itertools, json
import numpy as np, scipy.sparse as sp
r0 = float(sys.argv[1]); model = sys.argv[2]; L0 = int(sys.argv[3]); NSW = int(sys.argv[4])
START = sys.argv[5] if len(sys.argv) > 5 else 'chaud'; AN = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
rng = np.random.default_rng()
NLEV = int(math.log2(L0)) - 1
S = lambda n: 10 * 4**n + 4 * 2**n + 1
H = lambda n: 2**(n + 1) + 1
def weights(l):
    if model == 'A': return dict(B=r0 * H(0) / H(l), T=r0 * S(0) / S(l), S=r0 * S(0) / S(l + 1))
    if model == 'B': return dict(B=r0, T=r0, S=r0 * S(l) / S(l + 1))
    if model == 'U': return dict(B=r0, T=r0, S=r0)
    if model == 'X': return dict(B=r0, T=AN * r0, S=AN * r0)
Lv = [L0 >> l for l in range(NLEV + 1)]
# indices de liens : h(l,x,k) et v(l,x)
off = {}; n = 0
for l in range(NLEV + 1):
    off[('h', l)] = n; n += 3 * Lv[l]**3
    if l < NLEV: off[('v', l)] = n; n += Lv[l]**3
nl = n
def sidx(l, x): L = Lv[l]; return ((x[0] % L) * L + (x[1] % L)) * L + (x[2] % L)
def h(l, x, k): return off[('h', l)] + 3 * sidx(l, x) + k
def v(l, x): return off[('v', l)] + sidx(l, x)
def sh(x, k, s=1): y = list(x); y[k] += s; return tuple(y)
par = lambda l, x: tuple((c % Lv[l]) // 2 for c in x)
rows, cols, typ, lev = [], [], [], []
pB, pE = {}, {}
for l in range(NLEV + 1):
    for x in itertools.product(range(Lv[l]), repeat=3):
        for i, j in [(0, 1), (0, 2), (1, 2)]:
            p = len(typ); pB[(l, sidx(l, x), i, j)] = p; typ.append('B'); lev.append(l)
            for a in (h(l, x, i), h(l, sh(x, i), j), h(l, sh(x, j), i), h(l, x, j)): rows.append(p); cols.append(a)
        if l < NLEV:
            for k in range(3):
                y = sh(x, k); P, Pp = par(l, x), par(l, y)
                ls = [h(l, x, k), v(l, y), v(l, x)]
                t = 'T' if P == Pp else 'S'
                if t == 'S': ls.append(h(l + 1, P, k))
                p = len(typ); pE[(l, sidx(l, x), k)] = p; typ.append(t); lev.append(l)
                for a in ls: rows.append(p); cols.append(a)
npq = len(typ)
A = sp.csr_matrix((np.ones(len(rows), dtype=np.int32), (rows, cols)), shape=(npq, nl))
A.data[:] = 1
W = [weights(l) for l in range(NLEV + 1)]
logr = np.array([math.log(W[lev[p]][typ[p]]) for p in range(npq)])
# coloration gloutonne des liens (aucune plaquette commune dans une classe)
At = A.T.tocsr(); conf = (At @ A).tocsr()
color = -np.ones(nl, dtype=np.int64)
for a in rng.permutation(nl):
    used = set(color[conf.indices[conf.indptr[a]:conf.indptr[a + 1]]].tolist())
    c = 0
    while c in used: c += 1
    color[a] = c
classes = [np.where(color == c)[0] for c in range(color.max() + 1)]
Acs = [A[:, C].tocsc() for C in classes]; AcsT = [M.T.tocsr() for M in Acs]
b = rng.integers(0, 2, nl).astype(np.int8) if START == 'chaud' else np.zeros(nl, dtype=np.int8)
F = (A @ b.astype(np.int64)) % 2
def sweep():
    global F
    for C, M, MT in zip(classes, Acs, AcsT):
        d = MT @ ((1 - 2 * F) * logr)
        acc = (d >= 0) | (rng.random(len(C)) < np.exp(np.minimum(d, 0)))
        if acc.any():
            F = (F + M[:, acc] @ np.ones(acc.sum(), dtype=np.int64)) % 2
# chiffre de temps : paires (plaquette, plaquette) par site
qpairs = []; parent_of = []; tau = []
for l in range(NLEV):
    P1, P2 = [], []
    for x in itertools.product(range(Lv[l]), repeat=3):
        a1, a2 = [], []
        for k in range(3):
            i, j = [t for t in range(3) if t != k]
            Pp = par(l, sh(x, k))
            a1.append(pE[(l, sidx(l, x), k)]); a2.append(pB[(l + 1, sidx(l + 1, Pp), i, j)])
            a1.append(pB[(l, sidx(l, x), i, j)]); a2.append(pE[(l, sidx(l, sh(sh(x, i), j)), k)])
        P1.append(a1); P2.append(a2)
    qpairs.append((np.array(P1), np.array(P2)))
    parent_of.append(np.array([sidx(l + 1, par(l, x)) for x in itertools.product(range(Lv[l]), repeat=3)]))
    tau.append(np.array([sum(1 for c in x if c % 2 == 0) % 2 for x in itertools.product(range(Lv[l]), repeat=3)]))
blocks = []
for l in range(NLEV):
    L = Lv[l]; B = []
    for X in itertools.product(range(L // 2), repeat=3):
        B.append([sidx(l, tuple(2 * c + e for c, e in zip(X, d))) for d in itertools.product((0, 1), repeat=3)])
    blocks.append(np.array(B))
for _ in range(NSW // 4): sweep()
acc_q = [np.zeros(2) for _ in range(NLEV)]; acc_L = [np.zeros(2) for _ in range(NLEV)]; acc_R = [np.zeros(2) for _ in range(NLEV)]
joint = [np.zeros(4) for _ in range(NLEV)]; blk = [np.zeros(256) for _ in range(NLEV)]
flux = {t: [0, 0] for t in 'BTS'}; fluxl = {(l, t): [0, 0] for l in range(NLEV + 1) for t in 'BTS'}; typ_arr = np.array(typ); lev_arr = np.array(lev)
nmeas = 0
for s in range(NSW):
    sweep()
    if s % 2: continue
    nmeas += 1
    Q = [((F[P1] * F[P2]).sum(1)) % 2 for P1, P2 in qpairs]
    for l in range(NLEV):
        q = Q[l]; acc_q[l] += [q.sum(), q.size]
        acc_L[l] += [q[tau[l] == 0].sum(), (tau[l] == 0).sum()]; acc_R[l] += [q[tau[l] == 1].sum(), (tau[l] == 1).sum()]
        code = (q[blocks[l]] * (1 << np.arange(8))).sum(1); blk[l] += np.bincount(code, minlength=256)
        if l < NLEV - 1:
            joint[l] += np.bincount(2 * q + Q[l + 1][parent_of[l]], minlength=4)
    for l in range(NLEV + 1):
        ml = lev_arr == l
        for t in 'BTS':
            mk = ml & (typ_arr == t)
            if mk.any(): fluxl[(l, t)][0] += F[mk].sum(); fluxl[(l, t)][1] += mk.sum()
    m0 = lev_arr == 0
    for t in 'BTS':
        mk = m0 & (typ_arr == t); flux[t][0] += F[mk].sum(); flux[t][1] += mk.sum()
def Hn(c):
    c = c[c > 0] / c.sum(); return float(-(c * np.log2(c)).sum())
Hb = lambda p: 0.0 if p in (0, 1) else -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
out = dict(r0=r0, modele=model, depart=START, a=AN, L0=L0, balayages=NSW, mesures=nmeas, classes=len(classes),
           flux_niveau0={t: round(f[0] / f[1], 4) for t, f in flux.items() if f[1]}, flux_par_niveau={'%d%s' % k: round(f[0] / f[1], 5) for k, f in fluxl.items() if f[1]}, niveaux=[])
for l in range(NLEV):
    p = acc_q[l][0] / acc_q[l][1]
    rec = dict(niveau=l, poids={k: round(vv, 5) for k, vv in W[l].items()}, p_q1=round(p, 4),
               p_L=round(acc_L[l][0] / acc_L[l][1], 4), p_R=round(acc_R[l][0] / acc_R[l][1], 4),
               H_marg=round(Hb(p), 4), H_bloc_par_site=round(Hn(blk[l]) / 8, 4))
    if l < NLEV - 1:
        j = joint[l]; rec['H_cond_parent'] = round(Hn(j) - Hn(np.array([j[0] + j[2], j[1] + j[3]])), 4)
    out['niveaux'].append(rec)
print(json.dumps(out, ensure_ascii=False))
