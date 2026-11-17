"""Relecture V167 : socle. Charge les definitions de V162_FHT_R3_jauge.py (modele de reference, non modifie),
reconstruit H(theta) par composantes de Fourier de facon independante (moindres carres sur theta aleatoires,
frequences candidates k/(2L), |k| <= KMAX choisi plus grand que celui de V166), et fournit H, dH/dtheta.
Usage : import rv3_core as C ; C.setup(N, q)"""
import sys, os, numpy as np, scipy.sparse as sp
CODE = '/home/claude/fht/05_code/code_v30/V162_FHT_R3_jauge.py'
G = {}
def setup(N, q, m=-1.0, TR='b', ROOT='r1'):
    sys.argv = ['x', str(N), str(m), TR, str(q), 'none', '1', '64', ROOT, '0/1']
    src = open(CODE).read().split("print(f'=== V162")[0]
    g = {'__name__': 'v162'}; exec(src, g); G.clear(); G.update(g)
    return g
def Href(th):
    return G['hamiltonian'](th, 0.0)[0]
def fourier_fit(kmax_factor=2, nsamp=None, seed=7):
    """Moindres carres sur des theta aleatoires, avec un ensemble de frequences PLUS LARGE que celui de V166
    (|k| <= kmax_factor * 2 q L). Si V166 est complet, les composantes hors {0, +-2q 2^j} doivent etre nulles."""
    L, q, Nc = G['L'], G['q'], G['Nc']; n4 = 4 * Nc
    kmax = int(round(kmax_factor * 2 * L * abs(q))); ks = np.arange(-kmax, kmax + 1); K = len(ks)
    ns = nsamp or (K + 40)
    rng = np.random.default_rng(seed); ths = rng.random(ns) * 4 * np.pi * L * kmax_factor
    h0 = Href(0.777).tocoo(); keys = np.unique(h0.row.astype(np.int64) * n4 + h0.col)
    D = np.zeros((ns, len(keys)), complex)
    for i, t in enumerate(ths):
        h = Href(t).tocoo(); kk = h.row.astype(np.int64) * n4 + h.col; pos = np.searchsorted(keys, kk)
        assert np.all(keys[np.minimum(pos, len(keys) - 1)] == kk)
        np.add.at(D[i], pos, h.data)
    A = np.exp(1j * np.outer(ths, ks) / (2 * L))
    Cf, res, rk, sv = np.linalg.lstsq(A, D, rcond=None)
    return ks, Cf, keys, n4, ths, D, A
class Fast:
    def __init__(self, ks, Cf, keys, n4, tol=1e-12):
        keep = np.abs(Cf).max(axis=1) > tol
        self.ks, self.F = ks[keep], Cf[keep]; self.rows = (keys // n4).astype(int); self.cols = (keys % n4).astype(int); self.n4 = n4
        self.L = G['L']
    def H(self, th, deriv=0):
        ph = np.exp(1j * th * self.ks / (2 * self.L)) * (1j * self.ks / (2 * self.L)) ** deriv
        return sp.csr_matrix((ph @ self.F, (self.rows, self.cols)), shape=(self.n4, self.n4))
def prof(v):
    Nc, N, depth = G['Nc'], G['N'], G['depth']
    return np.bincount(depth, weights=np.sum(np.abs(v.reshape(Nc, 4)) ** 2, axis=1), minlength=N + 1)
def profs(V):
    Nc, N, depth = G['Nc'], G['N'], G['depth']
    w = np.sum(np.abs(V.reshape(Nc, 4, -1)) ** 2, axis=1)   # Nc x k
    out = np.zeros((N + 1, V.shape[1]))
    for n in range(N + 1): out[n] = w[depth == n].sum(axis=0)
    return out
