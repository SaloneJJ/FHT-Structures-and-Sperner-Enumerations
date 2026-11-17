"""V339 relecture -- outils communs (écrits indépendamment de V335-V338).
Octet BBH construit en espace réel (sauts t_f intra-cellule, t_c inter-cellule, signes de Jordan-Wigner eta_d(s) = (-1)^{s_1+..+s_{d-1}}),
matrices de Clifford obtenues par dérivation numérique de H(k) au point du cône, et non recopiées des scripts relus.
Seule la géométrie de l'arbre augmenté (niveaux, liens U(1), transport (b) r0) est reprise de V162 (exec), comme autorisé."""
import numpy as np, itertools as it

SITES = list(it.product((0, 1), repeat=3))
IDX = {s: i for i, s in enumerate(SITES)}
SUB = np.array([(-1) ** sum(s) for s in SITES], float)          # marque L/R de l'octet (sous-réseau)

def octet_mats(tf=1.0, tc=1.0):
    """H0 (intra), M[d] (saut de la cellule x+e_d vers x) : (Hψ)_x = H0 ψ_x + Σ_d M_d ψ_{x+e_d} + M_d^† ψ_{x-e_d}"""
    H0 = np.zeros((8, 8)); M = [np.zeros((8, 8)) for _ in range(3)]
    for s in SITES:
        for d in range(3):
            eta = (-1) ** sum(s[:d])
            if s[d] == 0:
                t = list(s); t[d] = 1; t = tuple(t)
                H0[IDX[s], IDX[t]] += tf * eta; H0[IDX[t], IDX[s]] += tf * eta
            else:                                             # site s (s_d = 1) en x  <->  site s - e_d en x + e_d
                t = list(s); t[d] = 0; t = tuple(t)
                M[d][IDX[s], IDX[t]] += tc * eta
    return H0, M

def Hk(k, tf=1.0, tc=1.0, alt=False):
    H0, M = octet_mats(tf, tc); sg = -1 if alt else 1
    return H0 + sum(sg * (M[d] * np.exp(1j * k[d]) + M[d].T * np.exp(-1j * k[d])) for d in range(3))

def clifford_own(alt=False):
    """Γ5 (masse au cône), α_d = ∂H/∂q_d au cône, Γ4 orienté comme le modèle emprunté, A, B, F = iAB."""
    K0 = np.zeros(3) if alt else np.full(3, np.pi)
    Hm = Hk(K0, 0.4, 1.0, alt)                                 # masse au cône pour t_f != t_c
    G5 = Hm / np.sqrt(np.real(np.trace(Hm @ Hm)) / 8)
    e = 1e-5; al = []; cs = []
    for d in range(3):
        dq = np.zeros(3); dq[d] = e
        al.append((Hk(K0 + dq, 1, 1, alt) - Hk(K0 - dq, 1, 1, alt)) / (2 * e))
        cs.append((Hk(K0 + dq, 1, 1, alt) + Hk(K0 - dq, 1, 1, alt) - 2 * Hk(K0, 1, 1, alt)) / e ** 2)
    cs = [c / np.sqrt(np.real(np.trace(c @ c)) / 8) for c in cs]
    A = (cs[0] - cs[1]) / np.sqrt(2); B = (cs[0] + cs[1] - 2 * cs[2]) / np.sqrt(6)
    P = G5 @ al[0] @ al[1] @ al[2]
    return dict(G5=G5, al=al, cs=cs, A=A, B=B, P=P)

def borrowed_gammas():
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1]).astype(complex); I2 = np.eye(2)
    return [np.kron(sx, sx), np.kron(sx, sy), np.kron(sx, sz), np.kron(sy, I2), np.kron(sz, I2)]
