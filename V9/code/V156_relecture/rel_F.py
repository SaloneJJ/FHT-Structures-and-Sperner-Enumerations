"""Modele F reconstruit independamment a partir de la DESCRIPTION de V146 §0 (pas du code des auteurs), graine propre.
H_F = I (+) II, I = A (x) B (x) 3 generations, A = (s_a, d1, d2, s_b), gamma_A = (+,-,-,+), B = (lepton, 3 couleurs).
pi(lambda, q, m) : sur I agit sur l'indice A par diag(lambda, q, conj(lambda)) ; sur II agit sur l'indice B par diag(lambda, m).
J_F = echange I <-> II avec conjugaison ; gamma_F(II) = -gamma_F(I) image.
D_F : Yukawa s_a <-> d1 (Y_up), s_b <-> d2 (Y_down) dans I, conjugue dans II, + Majorana symetrique sur nu_R (a = 0, b = 0)."""
import itertools, numpy as np
nA, nB, nG = 4, 4, 3; dI = nA * nB * nG; dF = 2 * dI
def idx(a, b, g, part=0): return part * dI + (a * nB + b) * nG + g
gA = np.array([1, -1, -1, 1])
gF = np.array([gA[a] * (1 if part == 0 else -1) for part in (0, 1) for a in range(nA) for b in range(nB) for g in range(nG)], float)
GF = np.diag(gF)
UF = np.zeros((dF, dF)); UF[:dI, dI:] = np.eye(dI); UF[dI:, :dI] = np.eye(dI)
def Jconj(U, X): return U @ X.conj() @ U.conj().T          # J X J^-1 pour J = U o conj
def pi(lam, q, m):
    P = np.zeros((dF, dF), complex)
    MA = np.zeros((4, 4), complex); MA[0, 0] = lam; MA[1:3, 1:3] = q; MA[3, 3] = np.conj(lam)
    P[:dI, :dI] = np.kron(np.kron(MA, np.eye(nB)), np.eye(nG))
    MB = np.zeros((4, 4), complex); MB[0, 0] = lam; MB[1:, 1:] = m
    P[dI:, dI:] = np.kron(np.kron(np.eye(nA), MB), np.eye(nG))
    return P
pauli = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1., -1.]).astype(complex)]
def real_basis():
    """base reelle de A_F = C (+) H (+) M3(C) : 2 + 4 + 18 = 24 elements"""
    B = [pi(1, np.zeros((2, 2)), np.zeros((3, 3))), pi(1j, np.zeros((2, 2)), np.zeros((3, 3)))]
    for q in [np.eye(2)] + [1j * p for p in pauli]: B.append(pi(0, q, np.zeros((3, 3))))
    for i, j in itertools.product(range(3), repeat=2):
        for c in (1, 1j):
            mm = np.zeros((3, 3), complex); mm[i, j] = c; B.append(pi(0, np.zeros((2, 2)), mm))
    return B
def random_alg(rng):
    lam = complex(*rng.normal(size=2)); a, b = rng.normal(size=2) + 1j * rng.normal(size=2)
    return pi(lam, np.array([[a, b], [-np.conj(b), np.conj(a)]]), rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
def build_DF(rng, majorana=True):
    Y = {k: rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)) for k in ('n', 'e', 'u', 'd')}
    S = np.zeros((dI, dI), complex)
    for b in range(nB):
        yu, yd = (Y['n'], Y['e']) if b == 0 else (Y['u'], Y['d'])
        for g, h in itertools.product(range(nG), repeat=2):
            S[idx(1, b, g), idx(0, b, h)] += yu[g, h]; S[idx(2, b, g), idx(3, b, h)] += yd[g, h]
    S = S + S.conj().T
    D = np.zeros((dF, dF), complex); D[:dI, :dI] = S; D[dI:, dI:] = S.conj()
    if majorana:
        YR = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); YR = YR + YR.T
        T = np.zeros((dI, dI), complex)
        for g, h in itertools.product(range(nG), repeat=2): T[idx(0, 0, g), idx(0, 0, h)] = YR[g, h]
        D[dI:, :dI] = T; D[:dI, dI:] = T.conj().T
    return D
