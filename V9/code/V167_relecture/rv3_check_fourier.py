"""Point 1 : exactitude et completude de la reconstruction de Fourier de V166, derivee analytique.
Usage : python3 rv3_check_fourier.py N q"""
import sys, time, numpy as np, scipy.sparse as sp
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
N, q = int(sys.argv[1]), float(sys.argv[2]); t0 = time.time()
g = C.setup(N, q); L = g['L']
ks, Cf, keys, n4, ths, D, A = C.fourier_fit(kmax_factor=2 if N <= 3 else 1.25)
amp = np.abs(Cf).max(axis=1); nz = ks[amp > 1e-10]
print(f'=== N = {N}, q = {q}, L = {L}, dim {n4} ; moindres carres sur {len(ths)} theta aleatoires, |k| <= {ks.max()} (V166 : {int(2*L*q)}) ===')
print(f'  residu d ajustement max {np.abs(A @ Cf - D).max():.1e}')
print(f'  frequences k (unite 1/2L) avec |C_k| > 1e-10 : {list(nz)}')
print(f'  attendu [T] : 0 et +-2q*2^j, j = 0..N : {sorted(set([0] + [s * int(2*q*2**j) for j in range(N+1) for s in (1,-1)]))}')
print(f'  plus grande composante hors de cet ensemble : {amp[~np.isin(ks, nz)].max() if (~np.isin(ks, nz)).any() else 0:.1e}')
# reconstruction V166 (DFT sur 4qL+1 points), re-implementee a l'identique pour la comparer a la mienne
kmax = int(round(2 * L * abs(q))); M = 2 * kmax + 1; thm = 2 * np.pi * 2 * L * np.arange(M) / M; kv = np.arange(-kmax, kmax + 1)
Fv = np.zeros((len(kv), len(keys)), complex)
for i, t in enumerate(thm):
    h = C.Href(t).tocoo(); kk = h.row.astype(np.int64) * n4 + h.col; pos = np.searchsorted(keys, kk); d = np.zeros(len(keys), complex); np.add.at(d, pos, h.data)
    Fv += np.outer(np.exp(-2j * np.pi * kv * i / M), d) / M
dif = 0
for kk_, row in zip(kv, Fv):
    j = np.where(ks == kk_)[0][0]; dif = max(dif, np.abs(row - Cf[j]).max())
print(f'  |C_k(DFT V166) - C_k(moindres carres)| max {dif:.1e}')
fast = C.Fast(ks, Cf, keys, n4)
rng = np.random.default_rng(11); e0 = e1 = e2 = 0
for t in list(rng.random(6) * 2 * np.pi) + [0.0, np.pi, np.pi / q, 1e3 * np.pi + 0.3]:
    H = C.Href(t); e0 = max(e0, abs(fast.H(t) - H).max())
    hh = 1e-3; Dn = (-C.Href(t + 2*hh) + 8 * C.Href(t + hh) - 8 * C.Href(t - hh) + C.Href(t - 2*hh)) / (12 * hh)   # 5 points, erreur O(h^4)
    e1 = max(e1, abs(fast.H(t, 1) - Dn).max())
    D2 = (C.Href(t + hh) - 2 * H + C.Href(t - hh)) / hh**2; e2 = max(e2, abs(fast.H(t, 2) - D2).max())
print(f'  10 theta (dont 0, pi, pi/q, 1000 pi + 0.3) : |H_Fourier - H_V162| max {e0:.1e} ; dH/dtheta vs differences a 5 points (h = 1e-3) {e1:.1e} ; d2H vs differences {e2:.1e}')
print(f'  norme de dH/dtheta : max |element| {abs(fast.H(0.3, 1)).max():.3e}')
# periodicite H(theta + 2 pi / q) ~ H(theta) par jauge diagonale g = exp(i 2 pi Z / L) (Z coordonnee fine du representant)
print(f'  temps {time.time() - t0:.0f} s')
