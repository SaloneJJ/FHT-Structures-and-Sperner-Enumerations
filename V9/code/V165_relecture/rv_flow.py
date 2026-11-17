"""Flot spectral par couche sous flux (Phi = 1, torsion theta le long de z), trois estimateurs independants :
 M1 : fenetre de Hellmann-Feynman (comme V158, re-ecrit) : F(n) = int dth sum_{|E|<eps} E' w(n) / (2 eps) ; aussi pour eps/2 ;
 M2 : pompage adiabatique (Kubo, toutes les valeurs propres) : F(n) = - int dth sum_i w_i'(n) G(E_i),
      G(E) = clamp(E / 2eps, -1/2, 1/2) ; egal a M1 dans la limite continue (integration par parties), discretisation differente ;
 M3 : suivi des niveaux par recouvrement et comptage des traversees de E = 0 (signe x poids par couche au passage).
Usage : python3 rv_flow.py MODELE GEOM N [param] [NTH] [methodes] [jauge]
  MODELE : wil (V151, m = param) | ksO / ksW (KS + masse eps, M = param) | sgO / sgW (KS + masse singulet, M = param)
  methodes : sous-ensemble de '123' (defaut '13')"""
import sys, os, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_core import *
t0 = time.time()
MOD, GEOM, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
PAR = float(sys.argv[4]) if len(sys.argv) > 4 else (-1.0 if MOD == 'wil' else 1.0)
NTH = int(sys.argv[5]) if len(sys.argv) > 5 else 64
METH = sys.argv[6] if len(sys.argv) > 6 else '13'
GAUGE = sys.argv[7] if len(sys.argv) > 7 else 'G2'
g = Geo(N, GEOM == 'flat', GAUGE)
ncomp = 4 if MOD == 'wil' else 1
eta, eta4, eps_st = ks_phases(g)

def HD(th):
    if MOD == 'wil':
        return H_wilson(g, th, PAR), H_wilson(g, th, PAR, deriv=True)
    K, Kd = H_ks_kin(g, th), H_ks_kin(g, th, deriv=True)
    Mn = mass_profile(g, PAR, MOD.endswith('W'))
    if MOD.startswith('ks'):
        return (K + sp.diags(Mn * eps_st)).tocsr(), Kd
    h = 1e-5
    G5 = gamma5_singlet(g, th); dG5 = (gamma5_singlet(g, th + h) - gamma5_singlet(g, th - h)) / (2 * h)
    return (K + sp.diags(Mn) @ G5).tocsr(), (Kd + sp.diags(Mn) @ dG5).tocsr()

# reduction de Bloch en z (plat, modele wil : invariance par translation d'une cellule en z)
blocks = None
if GEOM == 'flat' and MOD == 'wil':
    L = g.L; z = g.coords[:, 2]; base = np.where(z == 0)[0]
    # cellule de base (n, x, y, 0) -> cellules (n, x, y, z)
    lookup = {(g.depth[i], g.coords[i, 0], g.coords[i, 1]): b for b, i in enumerate(base)}
    col = np.array([lookup[(g.depth[i], g.coords[i, 0], g.coords[i, 1])] for i in range(g.Nc)])
    blocks = []
    for mz in range(L):
        ph = np.exp(2j * np.pi * mz * z / L) / np.sqrt(L)
        F = sp.csr_matrix((np.repeat(ph, ncomp), (np.arange(g.Nc * ncomp), (col[:, None] * ncomp + np.arange(ncomp)).ravel())),
                          shape=(g.Nc * ncomp, len(base) * ncomp))
        blocks.append(F)
    depth_red = g.depth[base]
    def layer_w_red(V):
        w = (np.abs(V) ** 2).reshape(len(base), ncomp, -1).sum(1)
        return np.array([w[depth_red == n].sum(0) for n in range(N + 1)])

EPS = 0.25 * 2 * np.pi / g.L
def Gs(E, e): return np.clip(E / (2 * e), -0.5, 0.5)

def spectrum(th, full):
    """renvoie E, V (vecteurs dans l'espace complet), D en base propre partielle (diag), et, si full, matrice D_ij et P_n"""
    H, D = HD(th)
    if blocks is not None:
        Es, Vs, Ds = [], [], []
        for F in blocks:
            Hk = (F.getH() @ H @ F).toarray(); Dk = (F.getH() @ D @ F).toarray()
            E, W = np.linalg.eigh(Hk)
            if not full: s_ = np.abs(E) < 3 * EPS; E, W = E[s_], W[:, s_]
            Es.append(E); Vs.append(W); Ds.append(Dk)
        return Es, Vs, Ds
    if full or H.shape[0] <= 1200:
        E, V = np.linalg.eigh(H.toarray())
        if not full: s_ = np.abs(E) < 3 * EPS; E, V = E[s_], V[:, s_]
        return [E], [V], [D]
    k = 40
    E, V = sla.eigsh(H, k=k, sigma=1e-5, which='LM')
    if np.sum(np.abs(E) < 3 * EPS) >= k - 2: print('  ATTENTION fenetre saturee', flush=True)
    s_ = np.abs(E) < 3 * EPS
    return [E[s_]], [V[:, s_]], [D]

ths = np.linspace(0, 2 * np.pi, NTH + 1)          # point final 2 pi inclus (M3 a besoin de fermer la boucle)
dth = ths[1] - ths[0]
F1 = np.zeros(N + 1); F1h = np.zeros(N + 1); G1b = np.zeros(N + 1); F2 = np.zeros(N + 1); F3 = np.zeros(N + 1); n3 = [0, 0]
prev = None; emin_edge = np.inf
full = '2' in METH
for it, th in enumerate(ths):
    Es, Vs, Ds = spectrum(th, full)
    cur = []
    for E, V, D in zip(Es, Vs, Ds):
        if len(E) == 0: cur.append((E, V, np.zeros((N + 1, 0)))); continue
        w = layer_weights(g, V, ncomp) if blocks is None else layer_w_red(V)
        cur.append((E, V, w))
        if it < NTH:
            DV = D @ V
            sl = np.real(np.einsum('ij,ij->j', V.conj(), DV)) if not full else None
            win = np.abs(E) < EPS
            if win.any():
                slw = np.real(np.einsum('ij,ij->j', V[:, win].conj(), DV[:, win]))
                F1 += dth / (2 * EPS) * (w[:, win] * slw).sum(1); G1b += dth / (2 * EPS) * (w[:, win] * np.abs(slw)).sum(1)
            win2 = np.abs(E) < EPS / 2
            if win2.any():
                slw = np.real(np.einsum('ij,ij->j', V[:, win2].conj(), DV[:, win2]))
                F1h += dth / EPS * (w[:, win2] * slw).sum(1)
            if full:
                # M2 : - sum_i G(E_i) w_i'(n),  w_i'(n) = 2 Re sum_{j != i} <i|P_n|j><j|D|i>/(E_i - E_j)
                Dm = V.conj().T @ DV
                dE = E[:, None] - E[None, :]; GG = Gs(E, EPS); Gp = np.where(np.abs(E) < EPS, 1 / (2 * EPS), 0.0)
                deg = np.abs(dE) < 1e-10
                fac = np.where(deg, 0.5 * (Gp[:, None] + Gp[None, :]), (GG[:, None] - GG[None, :]) / np.where(deg, 1.0, dE))
                np.fill_diagonal(fac, 0.0)   # sum_{i,j} G_i (...)/(Ei-Ej) symetrise
                for n in range(N + 1):
                    Vn = V[(np.repeat(g.depth if blocks is None else depth_red, ncomp) == n), :]
                    Pn = Vn.conj().T @ Vn
                    # sum_i G_i 2Re sum_j Pn_ij D_ji/(E_i-E_j) = sum_{i<j} 2Re(Pn_ij D_ji)(G_i-G_j)/(E_i-E_j)
                    F2[n] -= dth * np.sum(np.real(Pn * Dm.T) * fac)
    # M3 : appariement par recouvrement avec le pas precedent, comptage des changements de signe
    if prev is not None:
        for (Ep, Vp, wp), (Ec, Vc, wc) in zip(prev, cur):
            sp_ = np.abs(Ep) < 2 * EPS; sc_ = np.abs(Ec) < 2 * EPS
            if not sp_.any() or not sc_.any(): continue
            O = np.abs(Vp[:, sp_].conj().T @ Vc[:, sc_]) ** 2
            r_, c_ = linear_sum_assignment(-O)
            ip, ic = np.where(sp_)[0][r_], np.where(sc_)[0][c_]
            for a, b, o in zip(ip, ic, O[r_, c_]):
                if np.sign(Ep[a]) != np.sign(Ec[b]) and o > 0.3:
                    sgn = 1.0 if Ec[b] > Ep[a] else -1.0
                    F3 += sgn * 0.5 * (wp[:, a] + wc[:, b]); n3[0 if sgn > 0 else 1] += 1
    prev = cur
    for (E, V, w) in cur:
        if len(E):
            edge = (w[N] + (w[N - 1] if N >= 1 else 0)) > 0.5
            if edge.any(): emin_edge = min(emin_edge, np.min(np.abs(E[edge])))
print(f'=== relecture : modele {MOD}, {GEOM}, N = {N}, param {PAR}, jauge {GAUGE}, {NTH} pas, eps = {EPS:.3f} ===')
fmt = lambda F: '[' + ' '.join(f'{x:+.3f}' for x in F) + f'] total {F.sum():+.3f}'
print('  M1 (fenetre eps)   :', fmt(F1), ' brut', '[' + ' '.join(f'{x:.3f}' for x in G1b) + ']')
print('  M1 (fenetre eps/2) :', fmt(F1h))
if full: print('  M2 (pompage Kubo)  :', fmt(F2))
print('  M3 (traversees)    :', fmt(F3), f' montees {n3[0]}, descentes {n3[1]}')
print(f'  plus petit |E| porte par n >= N-1 (sur les theta, fenetre 3 eps) : {emin_edge:.4f}')
print(f'  temps {time.time() - t0:.0f} s', flush=True)
