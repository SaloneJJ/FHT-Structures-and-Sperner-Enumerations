"""V352 : action spectrale (noyau de la chaleur) du Dirac natif de la piste B sur l'arbre augmenté.
Constructeur numpy seul (sans scipy) de Hnat, mode « unif » de V336 (t_f = t_c = r = μ = 1, m_tot = 0, jauge alternée),
sans flux, avec phase de Bloch k0 de la racine (forêt d'arbres : extension latérale infinie, Bloch exact).
Critères : 01_synthese/V352 §2 (écrits avant calcul).                                    JJ SALONE, FHT 2026"""
import itertools, numpy as np

sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
Xg = [k3(sx, I2, I2), k3(sz, sx, I2), k3(sz, sz, sx)]; Yg = [k3(sy, I2, I2), k3(sz, sy, I2), k3(sz, sz, sy)]
G5 = sum(Xg) / np.sqrt(3); Am = (Xg[0] - Xg[1]) / np.sqrt(2); Bm = (Xg[0] + Xg[1] - 2 * Xg[2]) / np.sqrt(6)
G4 = G5 @ Yg[0] @ Yg[1] @ Yg[2]; FG = 1j * Am @ Bm @ G5; Y123 = Yg[0] @ Yg[1] @ Yg[2]
OPS = {'G5': G5, 'AY': Am @ Y123, 'BY': Bm @ Y123}

def geometry(N):
    sides = [2 ** n for n in range(N + 1)]
    offs = np.concatenate([[0], np.cumsum([s ** 3 for s in sides])]).astype(int)
    return sides, offs, int(offs[-1])

def cell_ops(N, k0=(0., 0., 0.), levels=None):
    """T[d] (sauts latéraux avant, par niveau) et S (enfant <- parent, 1/√8), denses, avec phase de Bloch :
    lien de déplacement Δx (unités de cellule fine) -> exp(i k0·Δx / 2^N) (jauge des centres ; équivalente à V162 sans flux)."""
    sides, offs, Nc = geometry(N); L = 2 ** N; k0 = np.asarray(k0, float)
    T = [np.zeros((Nc, Nc), complex) for _ in range(3)]; S = np.zeros((Nc, Nc), complex)
    lv = range(N + 1) if levels is None else levels
    for n in lv:
        s = sides[n]; d = L // s
        for I, J, K in itertools.product(range(s), repeat=3):
            i = offs[n] + (I * s + J) * s + K
            for jd, sh in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
                nb = offs[n] + (((I + sh[0]) % s) * s + (J + sh[1]) % s) * s + (K + sh[2]) % s
                T[jd][i, nb] += np.exp(-1j * k0[jd] * d / L)          # T[i, nb] = u(nb -> i), comme V162 (conj du lien avant)
    for n in range(1, N + 1):
        if levels is not None and (n not in levels or n - 1 not in levels): continue
        s = sides[n]; d = L // s
        for I, J, K in itertools.product(range(s), repeat=3):
            ch = offs[n] + (I * s + J) * s + K; pa = offs[n - 1] + ((I // 2) * (s // 2) + J // 2) * (s // 2) + K // 2
            dx = np.array([(I % 2) - 0.5, (J % 2) - 0.5, (K % 2) - 0.5]) * d       # centre enfant - centre parent (unités fines)
            S[ch, pa] = np.exp(1j * np.dot(k0, dx) / L) / np.sqrt(8)
    return T, S

def Hnat(N, k0=(0., 0., 0.), s=+1, mu=1.0, r=1.0, tf=1.0, tc=1.0, vertical=True, levels=None):
    T, S = cell_ops(N, k0, levels)
    if not vertical: S = 0 * S
    Nc = T[0].shape[0]; Id = np.eye(Nc); TCs = -tc
    H = np.zeros((8 * Nc, 8 * Nc), complex)
    for d in range(3):
        H += TCs * np.kron((T[d] - T[d].conj().T) / 2j, Yg[d])
        H += s * (tf * np.kron(Id, Xg[d]) + TCs * np.kron((T[d] + T[d].conj().T) / 2, Xg[d]))
    if vertical:
        H += np.kron((S - S.conj().T) / 2j, G4) + s * (r / 2) * np.kron(2 * Id - S - S.conj().T, G5)
    H += mu * np.kron(Id, FG)
    return H

def Hins(N, k0, epsD, O='AY', **kw):
    """D = D_+ ⊕ D_- + ε_D (O ⊗ 1) ⊗ σ_x (fibre jouet γ_F = σ_z, D_F = σ_x), comme V342/V343."""
    Hp, Hm = Hnat(N, k0, +1, **kw), Hnat(N, k0, -1, **kw)
    Nc = Hp.shape[0] // 8; X = epsD * np.kron(np.eye(Nc), OPS[O])
    return np.block([[Hp, X], [X, Hm]])

def depth_of(N):
    sides, offs, Nc = geometry(N)
    return np.concatenate([np.full(sd ** 3, n) for n, sd in enumerate(sides)])
