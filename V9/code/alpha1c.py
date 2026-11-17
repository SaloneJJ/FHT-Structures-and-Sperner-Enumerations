#!/usr/bin/env python3
"""alpha-1c : reponse de l'action spectrale de l'arbre augmente complet (Hnat de V352, point R-unif) a un champ U(1) lent.
Criteres : CRITERES_ALPHA1c.md (ecrits avant).  Derivee seconde exacte (Kubo), une diagonalisation par point de Bloch.
Usage :
  python3 alpha1c.py valid                      # V1 (Kubo = differences finies) et V2 (niveau fin = alpha1.py natif1 16)
  python3 alpha1c.py run N M B BY variante      # variante : couple | decouple | fin      (ex. : 2 4 4 2 couple)
numpy seul.                                                                                  JJ SALONE, FHT 2026"""
import sys, math, time, itertools
import numpy as np
import warnings
warnings.filterwarnings('ignore', message='.*encountered in matmul')   # fausses alertes d'Accelerate (macOS, numpy 2) ; resultats identiques sous Linux sans alerte
sys.dont_write_bytecode = True

sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
Xg = [k3(sx, I2, I2), k3(sz, sx, I2), k3(sz, sz, sx)]; Yg = [k3(sy, I2, I2), k3(sz, sy, I2), k3(sz, sz, sy)]
G5 = sum(Xg) / np.sqrt(3); Am = (Xg[0] - Xg[1]) / np.sqrt(2); Bm = (Xg[0] + Xg[1] - 2 * Xg[2]) / np.sqrt(6)
G4 = G5 @ Yg[0] @ Yg[1] @ Yg[2]; FG = 1j * Am @ Bm @ G5            # identiques a V352_arbre.py
LAMS = np.array([0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0])
KH = lambda Lm: math.sqrt(math.pi) / (12 * math.pi ** 2 * Lm)      # kappa chaleur (alpha1)
KE = lambda Lm: 1 / (24 * math.pi ** 2 * Lm)                       # kappa h (alpha1)

# ---------------------------------------------------------------- geometrie : supercellule de M racines selon y
def geometry(N, M, levels):
    L = 2 ** N; lv = []; off = 0
    for n in levels:
        s = 2 ** n; dims = (s, M * s, s); lv.append((n, s, dims, off, L // s)); off += s * M * s * s
    return L, lv, off

def links(N, M, levels, vertical):
    """liste des liens : (i, j, Delta = r_i - r_j (unites fines), y_milieu, type) ; type 0..2 = lateral (T[d]), 3 = vertical (S, enfant <- parent)."""
    L, lv, Nc = geometry(N, M, levels); idx = {}; out = []
    for (n, s, dims, off, d) in lv:
        for I, J, K in itertools.product(range(dims[0]), range(dims[1]), range(dims[2])):
            idx[(n, I, J, K)] = off + (I * dims[1] + J) * dims[2] + K
    for (n, s, dims, off, d) in lv:
        for I, J, K in itertools.product(range(dims[0]), range(dims[1]), range(dims[2])):
            i = idx[(n, I, J, K)]; y = (J + 0.5) * d
            for jd, sh in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
                nb = idx[(n, (I + sh[0]) % dims[0], (J + sh[1]) % dims[1], (K + sh[2]) % dims[2])]
                out.append((i, nb, -d * np.array(sh, float), y, jd))       # r_i - r_nb (non replie) = -d e_jd
    if vertical:
        nl = [x[0] for x in lv]
        for (n, s, dims, off, d) in lv:
            if n - 1 not in nl: continue
            for I, J, K in itertools.product(range(dims[0]), range(dims[1]), range(dims[2])):
                ch = idx[(n, I, J, K)]; pa = idx[(n - 1, I // 2, J // 2, K // 2)]
                dx = np.array([(I % 2) - 0.5, (J % 2) - 0.5, (K % 2) - 0.5]) * d   # r_enfant - r_parent
                yc = (J + 0.5) * d; ymid = yc - dx[1] / 2
                out.append((ch, pa, dx, ymid, 3))
    return L, Nc, out

def mats(N, M, L, Nc, lk, kvec, a=0.0, nq=1, deriv=False):
    """T[d], S a amplitude de champ a ; si deriv : aussi dT/da, d2T/da2 (en a = 0 + a)."""
    T = [np.zeros((Nc, Nc), complex) for _ in range(3)]; S = np.zeros((Nc, Nc), complex)
    T1 = [np.zeros((Nc, Nc), complex) for _ in range(3)]; S1 = np.zeros((Nc, Nc), complex)
    T2 = [np.zeros((Nc, Nc), complex) for _ in range(3)]; S2 = np.zeros((Nc, Nc), complex)
    q = 2 * np.pi * nq / (M * L)
    for (i, j, D, ym, t) in lk:
        phi = math.cos(q * ym) * D[0]                                 # phase de Peierls par unite de a : A_x(y_milieu) * (x_i - x_j)
        e = np.exp(1j * (np.dot(kvec, D) + a * phi))
        amp = 1.0 if t < 3 else 1 / np.sqrt(8)
        if t < 3:
            T[t][i, j] += amp * e; T1[t][i, j] += 1j * phi * amp * e; T2[t][i, j] += -phi ** 2 * amp * e
        else:
            S[i, j] += amp * e; S1[i, j] += 1j * phi * amp * e; S2[i, j] += -phi ** 2 * amp * e
    return (T, S, T1, S1, T2, S2) if deriv else (T, S)

def Hhop(T, S, s=+1, tc=1.0, r=1.0, vertical=True):
    """partie lineaire en (T, S) de Hnat (V352) ; la partie constante est ajoutee par Hconst."""
    Nc = T[0].shape[0]; H = np.zeros((8 * Nc, 8 * Nc), complex)
    for d in range(3):
        H += -tc * np.kron((T[d] - T[d].conj().T) / 2j, Yg[d]) + s * (-tc) * np.kron((T[d] + T[d].conj().T) / 2, Xg[d])
    if vertical:
        H += np.kron((S - S.conj().T) / 2j, G4) + s * (r / 2) * np.kron(-S - S.conj().T, G5)
    return H

def Hconst(Nc, s=+1, mu=1.0, r=1.0, tf=1.0, vertical=True):
    Id = np.eye(Nc); H = np.zeros((8 * Nc, 8 * Nc), complex)
    for d in range(3): H += s * tf * np.kron(Id, Xg[d])
    if vertical: H += s * r * np.kron(Id, G5)
    return H + mu * np.kron(Id, FG)

# ---------------------------------------------------------------- fonctions spectrales F(E) = chi(E^2/Lambda^2), F', F''
def F(E, Lm, kind):
    x = E / Lm
    if kind == 'heat': return np.exp(-x ** 2)
    ax = np.abs(x); return np.log1p(np.exp(-ax)) + ax / (1 + np.exp(ax))
def Fp(E, Lm, kind):
    x = E / Lm
    if kind == 'heat': return -2 * x / Lm * np.exp(-x ** 2)
    return -(x / Lm) / (4 * np.cosh(x / 2) ** 2)
def Fpp(E, Lm, kind):
    x = E / Lm
    if kind == 'heat': return (-2 + 4 * x ** 2) / Lm ** 2 * np.exp(-x ** 2)
    c2 = np.cosh(x / 2) ** 2; return (-1 / (4 * c2) + x * np.tanh(x / 2) / (4 * c2)) / Lm ** 2

def kubo(E, A, H2d, Lm, kind, blk=512):
    """d2/da2 Tr F(H) = sum_i F'(E_i) (H2)_ii + sum_ij |H1_ij|^2 (F'(E_i) - F'(E_j)) / (E_i - E_j)   (A = |H1_ij|^2 en base propre)
    v2 (15/11, memoire) : somme par blocs de lignes, meme formule ; ecart attendu a v1 : arrondi seulement"""
    G = Fp(E, Lm, kind); tot = float(np.dot(G, H2d))
    for a in range(0, len(E), blk):
        b = min(a + blk, len(E))
        dE = E[a:b, None] - E[None, :]; dG = G[a:b, None] - G[None, :]
        small = np.abs(dE) < 1e-7
        Dd = np.where(small, Fpp(0.5 * (E[a:b, None] + E[None, :]), Lm, kind), dG / np.where(small, 1.0, dE))
        tot += float(np.sum(A[a:b] * Dd))
    return tot

# ---------------------------------------------------------------- calcul principal
def response(N, M, B, BY, variante, mu=1.0, nqs=(1, 2), verbose=True):
    vertical = variante == 'couple'
    levels = [N] if variante == 'fin' else list(range(N + 1))
    L, Nc, lk = links(N, M, levels, vertical)
    V = (B * L) * (BY * M * L) * (B * L)
    acc = {(nq, kind): np.zeros(len(LAMS)) for nq in nqs for kind in ('heat', 'h')}
    E_all = []; t0 = time.time(); npts = B * B * BY; ip = 0
    for mx, my, mz in itertools.product(range(B), range(BY), range(B)):
        kvec = np.array([2 * np.pi * mx / (B * L), 2 * np.pi * my / (BY * M * L), 2 * np.pi * mz / (B * L)])
        E = U = None
        for nq in nqs:
            T, S, T1, S1, T2, S2 = mats(N, M, L, Nc, lk, kvec, nq=nq, deriv=True)
            if E is None:
                H0 = Hhop(T, S, vertical=vertical) + Hconst(Nc, mu=mu, vertical=vertical)
                E, U = np.linalg.eigh(H0); E_all.append(E); del H0          # v2 : H0 libere (memoire)
            H1 = Hhop(T1, S1, vertical=vertical); H1U = H1 @ U; del H1
            A = np.abs(U.conj().T @ H1U) ** 2; del H1U
            H2 = Hhop(T2, S2, vertical=vertical); H2U = H2 @ U; del H2
            H2d = np.real(np.einsum('ij,ij->j', U.conj(), H2U)); del H2U
            for kind in ('heat', 'h'):
                for il, Lm in enumerate(LAMS):
                    acc[(nq, kind)][il] += kubo(E, A, H2d, Lm, kind)
            del A
        ip += 1
        if verbose: print(f'   point de Bloch {ip}/{npts}  ({time.time() - t0:.0f} s)', flush=True)
    res = {}
    for kind, kap in (('heat', KH), ('h', KE)):
        ne = {}
        for nq in nqs:
            q = 2 * np.pi * nq / (M * L); ne[nq] = np.array([acc[(nq, kind)][i] / (kap(Lm) * V * q ** 2) for i, Lm in enumerate(LAMS)])
        if len(nqs) == 2:
            q1, q2 = [2 * np.pi * n / (M * L) for n in nqs]
            ne['ext'] = (ne[nqs[0]] * q2 ** 2 - ne[nqs[1]] * q1 ** 2) / (q2 ** 2 - q1 ** 2)
        res[kind] = ne
    return res, np.sort(np.concatenate(E_all)), acc

def total_trace(N, M, B, BY, variante, a, nq, Lm, kind, mu=1.0):
    vertical = variante == 'couple'; levels = [N] if variante == 'fin' else list(range(N + 1))
    L, Nc, lk = links(N, M, levels, vertical); tot = 0.0
    for mx, my, mz in itertools.product(range(B), range(BY), range(B)):
        kvec = np.array([2 * np.pi * mx / (B * L), 2 * np.pi * my / (BY * M * L), 2 * np.pi * mz / (B * L)])
        T, S = mats(N, M, L, Nc, lk, kvec, a=a, nq=nq)
        E = np.linalg.eigvalsh(Hhop(T, S, vertical=vertical) + Hconst(Nc, mu=mu, vertical=vertical)); tot += F(E, Lm, kind).sum()
    return tot

def validate():
    print('V1 : Kubo = differences finies (Richardson en a), arbre couple N = 1, M = 2, B = 2, BY = 1', flush=True)
    res, _, acc = response(1, 2, 2, 1, 'couple', verbose=False); worst = 0.0
    for kind in ('heat', 'h'):
        for il in (0, 5, 8):
            Lm = LAMS[il]; nq = 1; S0 = total_trace(1, 2, 2, 1, 'couple', 0.0, nq, Lm, kind); d = []
            for a in (0.02, 0.01):
                d.append((total_trace(1, 2, 2, 1, 'couple', a, nq, Lm, kind) + total_trace(1, 2, 2, 1, 'couple', -a, nq, Lm, kind) - 2 * S0) / a ** 2)
            fd = (4 * d[1] - d[0]) / 3; kb = acc[(nq, kind)][il]; rel = abs(fd - kb) / abs(kb); worst = max(worst, rel)
            print(f'   {kind:4s} Lambda = {Lm:4.2f} : Kubo {kb:.8e}  DF {fd:.8e}  ecart rel. {rel:.1e}')
    print(f'   V1 : ecart max {worst:.1e} -> {"PASSE" if worst < 1e-5 else "ECHEC"}')
    print('V2 : niveau fin seul (N = 2, M = 4, B = 4, BY = 1) contre alpha1.py natif1 L = 16', flush=True)
    import alpha1 as A1
    res2, E2, _ = response(2, 4, 4, 1, 'fin', verbose=False)
    d, mom, hop, ons_y, const = A1.terms('natif1'); L = 16; ks = 2 * np.pi * np.arange(L) / L; Ea = []
    for kx in ks:
        for kz in ks:
            H = np.zeros((8 * L, 8 * L), complex)
            for y in range(L):
                H[8 * y:8 * y + 8, 8 * y:8 * y + 8] += mom(0, kx) + mom(2, kz) + ons_y + const
                y2 = (y + 1) % L; H[8 * y:8 * y + 8, 8 * y2:8 * y2 + 8] += hop; H[8 * y2:8 * y2 + 8, 8 * y:8 * y + 8] += hop.conj().T
            Ea.append(np.linalg.eigvalsh(H))
    Ea = np.sort(np.concatenate(Ea)); dspec = np.max(np.abs(Ea - E2))
    print(f'   spectres : ecart max {dspec:.1e} -> {"PASSE" if dspec < 1e-8 else "ECHEC"}')
    # N_eff alpha1 (differences finies, meme q que n_q = 1, 2 ici : q = 2 pi n / 16)
    S0 = A1.spectral_sums('natif1', 16, 0.0, 1); worst = 0.0
    for nq in (1, 2):
        q = 2 * np.pi * nq / 16; out = []
        for a in (0.02, 0.01):
            Sp = A1.spectral_sums('natif1', 16, a, nq); Sm = A1.spectral_sums('natif1', 16, -a, nq); out.append([(Sp[j] + Sm[j] - 2 * S0[j]) / a ** 2 for j in (0, 1)])
        d2 = [(4 * out[1][j] - out[0][j]) / 3 for j in (0, 1)]
        for j, kind, kap in ((0, 'heat', KH), (1, 'h', KE)):
            nea = d2[j] / (np.array([kap(Lm) for Lm in LAMS]) * 16 ** 3 * q ** 2); mine = res2[kind][nq]
            rel = np.max(np.abs(nea - mine) / np.maximum(np.abs(mine), 1e-3)); worst = max(worst, rel)
            print(f'   n_q = {nq}, {kind:4s} : N_eff alpha1 {np.round(nea[[0, 3, 5, 8]], 5)}  ici {np.round(mine[[0, 3, 5, 8]], 5)}  ecart rel. max {rel:.1e}')
    print(f'   V2 (N_eff, precision de Richardson) : ecart max {worst:.1e}')

if __name__ == '__main__':
    if sys.argv[1] == 'valid':
        validate()
    else:
        N, M, B, BY = map(int, sys.argv[2:6]); var = sys.argv[6]; mu = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
        print(f'alpha-1c : N = {N} ({N + 1} niveaux), M = {M} racines selon y, Bloch {B} x {BY} x {B}, variante {var}, mu = {mu}', flush=True)
        t0 = time.time(); res, E, _ = response(N, M, B, BY, var, mu=mu)
        print(f'termine en {time.time() - t0:.0f} s ; |E|min = {np.min(np.abs(E)):.4f}')
        print('N_eff (extrapole q -> 0 ; entre parentheses n_q = 1, 2)')
        for i, Lm in enumerate(LAMS):
            print(f'  Lambda = {Lm:5.2f} : chaleur {res["heat"]["ext"][i]:8.3f} ({res["heat"][1][i]:7.3f}, {res["heat"][2][i]:7.3f}) ; '
                  f'h {res["h"]["ext"][i]:8.3f} ({res["h"][1][i]:7.3f}, {res["h"][2][i]:7.3f})')
        np.savez(f'ALPHA1c_N{N}_M{M}_B{B}x{BY}_{var}_mu{mu:g}.npz', LAMS=LAMS, heat=res['heat']['ext'], h=res['h']['ext'],
                 heat1=res['heat'][1], heat2=res['heat'][2], h1=res['h'][1], h2=res['h'][2], E=E)
