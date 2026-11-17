"""V110 - Test : la somme de q (cup de Serre) sur les sites R seuls est-elle un terme de bord ?
(retournement de liens interieurs). Usage : python3 V110_FHT_theta_test.py L0   -- JJ SALONE, FHT 2026"""
import sys
import numpy as np
from V109_FHT_complexe import Complexe, Cups
L0 = int(sys.argv[1]); K = Complexe(L0); C = Cups(K); rng = np.random.default_rng(2)
tau = [np.array([sum(1 for c in x if c % 2 == 0) % 2 for x in K.sites(l)]) for l in range(K.NLEV)]
def sums(b):
    tot = R = L = 0
    for l in range(K.NLEV):
        q = C.q_serre(b, l); tot += q.sum(); R += q[tau[l] == 1].sum(); L += q[tau[l] == 0].sum()
    return tot % 2, R % 2, L % 2
interior = [a for a in range(K.nl) if (K.link_type(a)[0] == 'v') or (0 < K.link_type(a)[1] < K.NLEV)]
ch = np.zeros(3); n = 300
for _ in range(n):
    b = rng.integers(0, 2, K.nl); s0 = sums(b); b[rng.choice(interior)] ^= 1; s1 = sums(b)
    ch += [s0[m] != s1[m] for m in range(3)]
print('fraction de changements (total, R seuls, L seuls) :', (ch / n).round(3))
