# V198 — renormalisation avec amplitudes par profondeur de fissure : t(1)=ε, t(2)=1, t(3)=τ. Critères : V197 §3.
import numpy as np, itertools as it
from scipy.optimize import least_squares
exec(open("V195_FHT_modes_coin.py").read().split('print("=== E1')[0])     # eta, cluster(n, tf, tc)
def hier(eps, tau, N=8):
    t = {1: eps, 2: 1.0, 3: tau}
    sites = list(it.product(range(N), repeat=3)); idx = {s: i for i, s in enumerate(sites)}
    H = np.zeros((N**3,)*2)
    for s in sites:
        for d in range(3):
            if s[d] + 1 < N:
                x = s[d]; h = 1
                while x % 2 == 1: h += 1; x //= 2
                u = list(s); u[d] += 1; H[idx[s], idx[tuple(u)]] = H[idx[tuple(u)], idx[s]] = t[h]*eta(d, s)
    return H
def eff_spec(p):
    H, _ = cluster(2, p[0], p[1]); return np.sort(np.linalg.eigvalsh(H))
print(" ε    τ    | 64e     65e    rapport | t_f'     t_c'     ε'      | résidu | E_bloc/√3  t_c'/τ")
for eps in [0.2, 0.3, 0.5]:
    Hb, _ = cluster(2, eps, 1.0); wb = np.sort(np.abs(np.linalg.eigvalsh(Hb))); Eb = wb[7]
    for tau in [0.05, 0.1, 0.2, 0.3, 0.5]:
        H = hier(eps, tau); w = np.linalg.eigvalsh(H); a = np.sort(np.abs(w)); o = np.argsort(np.abs(w))
        lowE = np.sort(w[o[:64]])
        best = None
        for g in [(Eb/np.sqrt(3), tau*0.5), (Eb/np.sqrt(3), tau), (Eb/np.sqrt(3), tau*0.2)]:
            r = least_squares(lambda p: eff_spec(p) - lowE, g, bounds=([0, 0], [10, 10]))
            if best is None or r.cost < best.cost: best = r
        tf2, tc2 = best.x; res = np.linalg.norm(eff_spec(best.x) - lowE)/np.linalg.norm(lowE)
        print(" %.1f  %.2f | %.4f  %.4f  %5.2f   | %.5f  %.5f  %.4f  | %.3f  | %.5f    %.4f" % (eps, tau, a[63], a[64], a[64]/a[63], tf2, tc2, tf2/tc2, res, Eb/np.sqrt(3), tc2/tau))
