"""Controle plat (paroi de Kaplan) recalcule par impulsion laterale, representation des Gamma propre.
Usage : python3 rel_flat.py N m"""
import sys, itertools, numpy as np
N = int(sys.argv[1]); m = float(sys.argv[2]); r = 1.0
s0 = np.eye(2); s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1., -1.]).astype(complex)
G = [np.kron(s1, s3), np.kron(s2, s3), np.kron(s3, s3), np.kron(s0, s1), np.kron(s0, s2)]
CHI = np.kron(s0, s3); L = 2 ** N; nl = N + 1
win = 1.1 * np.sin(2 * np.pi * 2 / L)
Sv = np.zeros((nl, nl)); Sv[np.arange(1, nl), np.arange(nl - 1)] = 1            # (S psi)(n) = psi(n-1)
Hv = np.kron(Sv, G[3] / 2j - r / 2 * G[4]); Hv = Hv + Hv.conj().T
cnt = {}; E1 = []
for n in itertools.product(range(L), repeat=3):
    k = 2 * np.pi * np.array(n) / L
    h = sum(np.sin(k[j]) * G[j] for j in range(3)) + (m + 4 * r - r * np.cos(k).sum()) * G[4]
    Hk = np.kron(np.eye(nl), h) + Hv
    E, V = np.linalg.eigh(Hk)
    # rotation dans les sous-espaces degeneres (paires +/-k ne se melangent pas ici ; degenerescences bord/couche 0 a k=0)
    for a in range(len(E)):
        if abs(E[a]) >= win: continue
        p = V[:, a].reshape(nl, 4); w = (np.abs(p) ** 2).sum(1); nb = np.arange(nl) @ w
        chi = np.real(np.einsum('ia,ab,ib->', p.conj(), CHI, p))
        key = ('couche N' if nb > N - 1.5 else ('couche 0' if nb < 1 else 'milieu'), '+' if chi > 0.3 else ('-' if chi < -0.3 else '0'))
        cnt[key] = cnt.get(key, 0) + 1
    if sum(min(x, L - x) ** 2 for x in n) == 1: E1.append(np.sort(np.abs(E))[0])
print(f'plat N = {N}, m = {m} : fenetre {win:.4f} :', cnt, f'; E(|n|^2=1) = {min(E1):.5f}')
