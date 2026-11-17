"""V109 - Tests de structure : bords (G A = 0, H4 G = 0), rangs sur F2, coordinations,
propriete de bord des trois cups (sum q invariante par retournement d'un lien interieur).
Usage : python3 V109_FHT_tests_structure.py L0   -- JJ SALONE, FHT 2026"""
import sys, collections
import numpy as np
from V109_FHT_complexe import Complexe, Cups, qpairs_v103

def rank_f2(M):
    M = (M.toarray() % 2).astype(np.uint8) if hasattr(M, 'toarray') else (M % 2).astype(np.uint8)
    M = M.copy(); r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0: continue
        p = r + piv[0]
        if p != r: M[[r, p]] = M[[p, r]]
        nz = np.nonzero(M[:, c])[0]; nz = nz[nz != r]
        M[nz] ^= M[r]
        r += 1
        if r == rows: break
    return r

L0 = int(sys.argv[1]); K = Complexe(L0)
print('L0', L0, 'niveaux', K.NLEV, 'liens', K.nl, 'plaquettes', K.np_, 'cellules3', K.n3, 'cellules4', K.n4)
print('G A = 0 :', ((K.G @ K.A).toarray() % 2 == 0).all(), ' H4 G = 0 :', ((K.H4 @ K.G).toarray() % 2 == 0).all())
if L0 <= 8:
    rA = rank_f2(K.A); rG = rank_f2(K.G); rH = rank_f2(K.H4)
    nsites = sum(K.Lv[l] ** 3 for l in range(K.NLEV + 1))
    print('rang A', rA, ' dim ker A^T', K.np_ - rA, ' rang G', rG, ' H2 (ker A^T / im G^T)', K.np_ - rA - rG)
    print('sites', nsites, ' rang de jauge (liens - rang A)', K.nl - rA, ' rang H4', rH, ' n3 - rang G', K.n3 - rG)
# coordinations par type et niveau
A = K.A.tocsc(); At = K.A.tocsr(); G = K.G.tocsc(); Gr = K.G.tocsr()
cnt = collections.Counter()
for p in range(K.np_):
    nlinks = At.indptr[p + 1] - At.indptr[p]; ncells = G.indptr[p + 1] - G.indptr[p]
    l = K.lev[p]; tag = 'bas' if l == 0 else ('haut' if l == K.NLEV else 'milieu')
    cnt[(K.typ[p], tag, nlinks, ncells)] += 1
print('plaquettes : (type, niveau, nb liens, nb cellules3) -> nombre')
for k in sorted(cnt): print('  ', k, cnt[k])
cl = collections.Counter()
for a in range(K.nl):
    kind, l = K.link_type(a); tag = 'bas' if l == 0 else ('haut' if l == K.NLEV else 'milieu')
    cl[(kind, tag, A.indptr[a + 1] - A.indptr[a])] += 1
print('liens : (type, niveau, nb plaquettes) -> nombre')
for k in sorted(cl): print('  ', k, cl[k])
c3 = collections.Counter()
for c in range(K.n3):
    l = K.lev3[c]; tag = 'bas' if l == 0 else ('haut' if l == K.NLEV else 'milieu')
    c3[(K.typ3[c], tag, Gr.indptr[c + 1] - Gr.indptr[c])] += 1
print('cellules3 : (type, niveau, nb plaquettes) -> nombre')
for k in sorted(c3): print('  ', k, c3[k])
# propriete de bord des cups
rng = np.random.default_rng(1)
C = Cups(K); QP = qpairs_v103(K)
def totals(b):
    F = (K.A @ b) % 2
    tv = sum(int(((F[P1] * F[P2]).sum(1) % 2).sum()) for P1, P2 in QP) % 2
    ts = sum(int(C.q_serre(b, l).sum()) for l in range(K.NLEV)) % 2
    ta = sum(int(C.q_aw(b, l).sum()) for l in range(K.NLEV)) % 2
    return tv, ts, ta
interior = [a for a in range(K.nl) if (K.link_type(a)[0] == 'v') or (0 < K.link_type(a)[1] < K.NLEV)]
chg = np.zeros(3); ntest = 200
for _ in range(ntest):
    b = rng.integers(0, 2, K.nl)
    t0 = totals(b); a = rng.choice(interior); b[a] ^= 1; t1 = totals(b)
    chg += [t0[m] != t1[m] for m in range(3)]
print('retournement d un lien interieur : fraction de cas ou sum q change (v103, serre, aw) :', (chg / ntest).round(3))
# egalite serre / aw a cobord pres : sum sur tout le complexe
same = 0
for _ in range(100):
    b = rng.integers(0, 2, K.nl); t = totals(b); same += (t[1] == t[2])
print('sum q serre = sum q aw (configurations aleatoires) :', same, '/ 100')
