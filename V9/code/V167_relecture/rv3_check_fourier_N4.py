"""Point 1 a N = 4 : modele a 11 frequences {0, +-2q 2^j} ajuste par moindres carres sur 30 theta aleatoires dans [0, 4 pi L) ;
le residu sur 30 points (pour 11 inconnues) teste la completude ; puis comparaison avec la DFT de V166 et derivee.
Usage : python3 rv3_check_fourier_N4.py q"""
import sys, time, numpy as np
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
q = float(sys.argv[1]); N = 4; t0 = time.time()
g = C.setup(N, q); L = g['L']; Nc = g['Nc']; n4 = 4 * Nc
ks = np.array(sorted(set([0] + [s * int(2*q*2**j) for j in range(N+1) for s in (1,-1)])))
h0 = C.Href(0.777).tocoo(); keys = np.unique(h0.row.astype(np.int64) * n4 + h0.col)
print(f'=== N = 4, q = {q}, dim {n4}, nnz {len(keys)} ; {len(ks)} frequences : {list(ks)} ===', flush=True)
def vec(t):
    h = C.Href(t).tocoo(); kk = h.row.astype(np.int64) * n4 + h.col; pos = np.searchsorted(keys, kk)
    assert np.all(keys[np.minimum(pos, len(keys) - 1)] == kk); d = np.zeros(len(keys), complex); np.add.at(d, pos, h.data); return d
rng = np.random.default_rng(5); ths = rng.random(30) * 4 * np.pi * L
D = np.array([vec(t) for t in ths]); A = np.exp(1j * np.outer(ths, ks) / (2 * L))
Cf = np.linalg.lstsq(A, D, rcond=None)[0]
print(f'  residu d ajustement (30 theta, 11 frequences) : {np.abs(A @ Cf - D).max():.1e} ; temps {time.time()-t0:.0f} s', flush=True)
del D
fast = C.Fast(ks, Cf, keys, n4)
# DFT V166 (4qL+1 points), on ne garde que les composantes non nulles
kmax = int(round(2 * L * q)); M = 2 * kmax + 1; kv = np.arange(-kmax, kmax + 1)
E = np.exp(-2j * np.pi * np.outer(kv, np.arange(M)) / M) / M
Fv = np.zeros((len(kv), len(keys)), complex)
for i in range(M):
    d = vec(2 * np.pi * 2 * L * i / M); Fv += np.outer(E[:, i], d)
on = np.isin(kv, ks)
print(f'  DFT V166 : composantes hors des 11 frequences : max {np.abs(Fv[~on]).max():.1e} ; ecart aux miennes {max(np.abs(Fv[kv == k][0] - Cf[ks == k][0]).max() for k in ks):.1e}', flush=True)
del Fv
e0 = e1 = 0
for t in [0.0, np.pi, 0.05, 1.234, 4.5]:
    H = C.Href(t); e0 = max(e0, abs(fast.H(t) - H).max()); hh = 1e-3
    Dn = (-C.Href(t + 2*hh) + 8 * C.Href(t + hh) - 8 * C.Href(t - hh) + C.Href(t - 2*hh)) / (12 * hh); e1 = max(e1, abs(fast.H(t, 1) - Dn).max())
print(f'  theta = 0, pi, 0.05, 1.234, 4.5 : |H_Fourier - H_V162| max {e0:.1e} ; dH/dtheta vs differences a 5 points {e1:.1e}')
np.savez(f'/home/claude/fht/relecture3/fourier_N4_q{int(q)}.npz', ks=ks, Cf=Cf, keys=keys, n4=n4)
print(f'  temps {time.time() - t0:.0f} s')
