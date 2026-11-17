#!/usr/bin/env python3
"""alpha-1a : multiplicite N_eff du terme de jauge (reponse de l'action spectrale a un champ U(1) lent). Criteres : CRITERES_ALPHA1.md.
Usage : python3 alpha1.py modele L   (modele : wilson | natif0 | natif1)"""
import sys, math, time
import numpy as np
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
Xg = [k3(sx, I2, I2), k3(sz, sx, I2), k3(sz, sz, sx)]; Yg = [k3(sy, I2, I2), k3(sz, sy, I2), k3(sz, sz, sy)]     # comme V352_arbre
G5 = sum(Xg) / np.sqrt(3); Am = (Xg[0] - Xg[1]) / np.sqrt(2); Bm = (Xg[0] + Xg[1] - 2 * Xg[2]) / np.sqrt(6); FG = 1j * Am @ Bm @ G5
al = [np.kron(sx, s) for s in (sx, sy, sz)]; be = np.kron(sz, I2)
def terms(model):
    """renvoie (dim, onsite(k_dir) pour une direction en impulsion, hop_fwd pour la direction y, onsite_y, constant)"""
    if model == 'wilson':
        d = 4; r = 1.0
        mom = lambda dd, k: np.sin(k) * al[dd] + r * (1 - np.cos(k)) * be
        hop = (al[1] / 2j - r / 2 * be); ons_y = r * be; const = np.zeros((4, 4))
    else:
        d = 8; mu = 0.0 if model == 'natif0' else 1.0
        mom = lambda dd, k: (1 - np.cos(k)) * Xg[dd] - np.sin(k) * Yg[dd]
        hop = (-1 / 2j) * Yg[1] - 0.5 * Xg[1]; ons_y = Xg[1]; const = mu * FG
    return d, mom, hop, ons_y, const
def chi_heat(u): return np.exp(-u)
def chi_h(u):
    x = np.sqrt(u); return np.log1p(np.exp(-x)) + x / (1 + np.exp(x))
LAMS = np.array([0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0])
def spectral_sums(model, L, a, n):
    d, mom, hop, ons_y, const = terms(model); ks = 2 * np.pi * np.arange(L) / L; ys = np.arange(L)
    A = a * np.cos(2 * np.pi * n * ys / L)
    Sh = np.zeros(len(LAMS)); Se = np.zeros(len(LAMS))
    for kx in ks:
        Hs = np.zeros((L, d * L, d * L), complex)
        for iz, kz in enumerate(ks):
            H = Hs[iz]
            for y in range(L):
                blk = mom(0, kx + A[y]) + mom(2, kz) + ons_y + const
                H[d * y:d * y + d, d * y:d * y + d] += blk
                y2 = (y + 1) % L; H[d * y:d * y + d, d * y2:d * y2 + d] += hop; H[d * y2:d * y2 + d, d * y:d * y + d] += hop.conj().T
        E2 = np.linalg.eigvalsh(Hs) ** 2
        for i, Lm in enumerate(LAMS):
            Sh[i] += chi_heat(E2 / Lm ** 2).sum(); Se[i] += chi_h(E2 / Lm ** 2).sum()
    return Sh, Se
if __name__ == '__main__':
    model, L = sys.argv[1], int(sys.argv[2]); t0 = time.time(); V = L ** 3
    S0 = spectral_sums(model, L, 0.0, 1)
    res = {}
    for n in (1, 2):
        q = 2 * np.pi * n / L; out = []
        for a in (0.02, 0.01):
            Sp = spectral_sums(model, L, a, n); Sm = spectral_sums(model, L, -a, n)
            out.append([(Sp[j] + Sm[j] - 2 * S0[j]) / a ** 2 for j in (0, 1)])
        d2 = [(4 * out[1][j] - out[0][j]) / 3 for j in (0, 1)]          # Richardson en a
        kh = np.sqrt(np.pi) / (12 * np.pi ** 2 * LAMS); ke = 1 / (24 * np.pi ** 2 * LAMS)
        res[n] = (d2[0] / (kh * V * q ** 2), d2[1] / (ke * V * q ** 2), q)
    q1, q2 = res[1][2], res[2][2]
    print(f'modele {model}, L = {L} ({time.time() - t0:.0f} s) ; N_eff extrapole a q -> 0 (lineaire en q^2) :')
    for i, Lm in enumerate(LAMS):
        ex = [(res[1][j][i] * q2 ** 2 - res[2][j][i] * q1 ** 2) / (q2 ** 2 - q1 ** 2) for j in (0, 1)]
        print(f'  Lambda = {Lm:5.2f} : chaleur N_eff = {ex[0]:7.3f} (q1 {res[1][0][i]:7.3f}, q2 {res[2][0][i]:7.3f}) ; h : N_eff = {ex[1]:7.3f} (q1 {res[1][1][i]:7.3f})')
