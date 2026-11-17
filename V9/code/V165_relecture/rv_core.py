"""Relecture independante (V159/V161/V163) -- noyau ecrit sans reprendre le code de V158/V162.
Geometrie : arbre augmente (cote 2^n au niveau n) ou plat (cote L a tous les niveaux), n = 0..N, T^3 de cote L = 2^N.
Jauge : potentiel continu lineaire A(r) = A0(r) + (0, 0, theta/L), avec deux jauges au choix
   G1 : A0 = (0, B X, 0), section psi(r + L e_x) = exp(i B L Y) psi(r)
   G2 : A0 = (-B Y, 0, 0), section psi(r + L e_y) = exp(-i B L X) psi(r)
Lien t <- s : U = exp(i int_{s~}^{t} A) * g(s, m), s~ = s + L m image de s adjacente a t ; integrale exacte (A lineaire).
Convention : (T psi)(t) = sum U psi(s) ; terme psi^dag(t) U psi(s).
"""
import itertools
import numpy as np, scipy.sparse as sp

class Geo:
    def __init__(s, N, flat=False, gauge='G1', q=1.0, Phi=1.0):
        s.N, s.flat, s.gauge, s.q = N, flat, gauge, q
        s.L = L = 2 ** N
        s.B = 2 * np.pi * Phi / L ** 2
        s.sides = [L if flat else 2 ** n for n in range(N + 1)]
        cells = []
        for n, sd in enumerate(s.sides):
            for x in itertools.product(range(sd), repeat=3):
                cells.append((n, x))
        s.cells = cells
        s.index = {c: i for i, c in enumerate(cells)}
        s.Nc = len(cells)
        s.depth = np.array([c[0] for c in cells])
        s.coords = np.array([c[1] for c in cells])
        s.dsz = np.array([L / s.sides[c[0]] for c in cells])
        s.cen = (s.coords + 0.5) * s.dsz[:, None]
        s._build_links()

    # potentiel (sans theta) et section
    def A0(s, r):
        B = s.B
        if s.gauge == 'G1':
            return np.array([0.0, B * r[0], 0.0])
        return np.array([-B * r[1], 0.0, 0.0])

    def section_phase(s, r, m):
        """psi(r + L m) = exp(i phase) psi(r)"""
        B, L = s.B, s.L
        if s.gauge == 'G1':
            return m[0] * B * L * r[1]
        return -m[1] * B * L * r[0]

    def link_phase0(s, t, s_img, s_cell):
        """phase sans theta du lien t <- s ; s_img = position de l'image de s adjacente a t"""
        rt = s.cen[t]
        m = np.round((s_img - s.cen[s_cell]) / s.L).astype(int)
        dl = rt - s_img
        integ = 0.5 * (s.A0(s_img) + s.A0(rt)) @ dl          # exact pour A lineaire
        return integ + s.section_phase(s.cen[s_cell], m), dl[2] / s.L   # coefficient de theta : int A_z = theta dz / L

    def _build_links(s):
        s.lat = []    # par direction j : (t, s, phase0, cz)
        for j in range(3):
            R, C, P, Z = [], [], [], []
            for i, (n, x) in enumerate(s.cells):
                sd = s.sides[n]
                y = list(x); y[j] = (y[j] + 1) % sd
                k = s.index[(n, tuple(y))]
                img = s.cen[i].copy(); img[j] += s.dsz[i]
                ph, cz = s.link_phase0(i, img, k)
                R.append(i); C.append(k); P.append(ph); Z.append(cz)
            s.lat.append(tuple(map(np.array, (R, C, P, Z))))
        R, C, P, Z = [], [], [], []
        for i, (n, x) in enumerate(s.cells):
            if n == 0: continue
            px = x if s.flat else tuple(v // 2 for v in x)
            p = s.index[(n - 1, px)]
            ph, cz = s.link_phase0(i, s.cen[p].copy(), p)
            R.append(i); C.append(p); P.append(ph); Z.append(cz)
        s.ver = tuple(map(np.array, (R, C, P, Z)))
        s.wv = 1.0 if s.flat else 1 / np.sqrt(8)

    def ops(s, th, deriv=False):
        """T_j (j = 0,1,2) et S (enfant <- parent, poids wv) ; deriv : d/dtheta"""
        q = s.q; out = []
        for (R, C, P, Z) in s.lat + [s.ver]:
            v = np.exp(1j * q * (P + th * Z))
            if deriv: v = 1j * q * Z * v
            out.append(sp.csr_matrix((v, (R, C)), shape=(s.Nc, s.Nc)))
        out[3] = out[3] * s.wv
        return out[:3], out[3]

# ------------------------------------------------------------------ modeles
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1]); I2 = np.eye(2)
GAM = [np.kron(sx, sx), np.kron(sx, sy), np.kron(sx, sz), np.kron(sy, I2), np.kron(sz, I2)]

def H_wilson(g, th, m=-1.0, r=1.0, deriv=False, T=None, S=None):
    """modele V151 (Gamma empruntes), secteur gamma_F = +1"""
    if T is None: T, S = g.ops(th, deriv)
    kr = lambda A, G: sp.kron(A, sp.csr_matrix(G), format='csr')
    H = sp.csr_matrix((4 * g.Nc, 4 * g.Nc), dtype=complex)
    for j in range(3):
        H = H + kr((T[j] - T[j].getH()) / 2j, GAM[j]) - (r / 2) * kr(T[j] + T[j].getH(), GAM[4])
    H = H + kr((S - S.getH()) / 2j, GAM[3]) - (r / 2) * kr(S + S.getH(), GAM[4])
    if not deriv: H = H + (m + 4 * r) * kr(sp.identity(g.Nc), GAM[4])
    return H.tocsr()

def ks_phases(g):
    x1, x2, x3 = g.coords.T
    eta = [np.ones(g.Nc), (-1.0) ** x1, (-1.0) ** (x1 + x2)]
    eta4 = (-1.0) ** (x1 + x2 + x3)            # evalue sur l'enfant (ligne de S)
    eps = (-1.0) ** (x1 + x2 + x3 + g.depth)
    return eta, eta4, eps

def H_ks_kin(g, th, deriv=False, eta4_mode='child'):
    T, S = g.ops(th, deriv); eta, eta4, eps = ks_phases(g)
    H = sp.csr_matrix((g.Nc, g.Nc), dtype=complex)
    for j in range(3):
        A = sp.diags(eta[j]) @ T[j]; H = H + (A - A.getH()) / 2j
    e4 = eta4 if eta4_mode == 'child' else np.ones(g.Nc)
    V = sp.diags(e4) @ S; H = H + (V - V.getH()) / 2j
    return H.tocsr()

def mass_profile(g, M, wall):
    return np.where(g.depth >= g.N - 1, -M, M) if wall else np.full(g.Nc, float(M))

def gamma5_singlet(g, th, a=(1, 0, 1, 0)):
    """f(x) * moyenne sur les 24 ordres de C1 C2 C3 C4, partie hermitienne ; x4 = n"""
    T, S = g.ops(th)
    C = [(T[j] + T[j].getH()) / 2 for j in range(3)] + [(S + S.getH()) / 2]
    P = None
    for perm in itertools.permutations(range(4)):
        Q = C[perm[0]] @ C[perm[1]] @ C[perm[2]] @ C[perm[3]]
        P = Q if P is None else P + Q
    x = np.column_stack([g.coords, g.depth])
    f = (-1.0) ** (x @ np.array(a))
    G = sp.diags(f) @ (P / 24)
    return ((G + G.getH()) / 2).tocsr()

def layer_weights(g, V, ncomp):
    """V : (dim, k) -> poids par couche (N+1, k)"""
    w = np.abs(V) ** 2
    w = w.reshape(g.Nc, ncomp, -1).sum(1)
    out = np.zeros((g.N + 1, w.shape[1]))
    for n in range(g.N + 1):
        out[n] = w[g.depth == n].sum(0)
    return out
