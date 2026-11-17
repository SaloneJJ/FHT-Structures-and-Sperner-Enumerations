"""V202 (derive de V158) : piste B native, paroi par dimerisation freres/cousins (MODELE ksD : PAR = rho_fin, 6e arg = rho_gros ; ksU : rho uniforme = PAR).
V158 : R2 (spin natif) -- flot spectral sous flux magnetique, (a) validation sur le modele emprunte de V151,
(b) modele natif de Kawamoto-Smit (une composante par cellule) sur l'arbre augmente et sur le controle plat.
Criteres et predictions : V158 (ecrits avant calcul).
Usage : python3 V158_FHT_R2_staggered.py MODELE GEOM N [param] [NTH]
   MODELE = bor (Gamma empruntes, V151, m = param, r = 1) | ksO (KS, masse uniforme M = param) | ksW (KS, paroi n >= N-1)
   GEOM = tree | flat
Flux : Phi = 1 quantum a travers (x, y) ; torsion theta le long de z, 0 -> 2 pi ; secteur gamma_F = +1.
Flot : F(n) = int dtheta sum_i delta_eps(E_i) dE_i/dtheta w_i(n), delta_eps = boite de demi-largeur eps (Hellmann-Feynman).
JJ SALONE, FHT 2026"""
import sys, time, itertools
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
t0 = time.time()
MOD, GEOM, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
PAR = float(sys.argv[4]) if len(sys.argv) > 4 else (-1.0 if MOD == 'bor' else 1.0)
NTH = int(sys.argv[5]) if len(sys.argv) > 5 else 96
RHOG = float(sys.argv[6]) if len(sys.argv) > 6 else 2.0
flat = GEOM == 'flat'; L = 2 ** N; B = 2 * np.pi / L ** 2; r = 1.0
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
GAM = [np.kron(sx, sx), np.kron(sx, sy), np.kron(sx, sz), np.kron(sy, I2), np.kron(sz, I2)]   # comme V151

# ---------------------------------------------------------------- geometrie
sides = [L if flat else 2 ** n for n in range(N + 1)]
offs = np.concatenate([[0], np.cumsum([s ** 3 for s in sides])]).astype(int); Nc = int(offs[-1])
def idx(n, x): s = sides[n]; return offs[n] + (x[0] % s * s + x[1] % s) * s + x[2] % s
coords = np.zeros((Nc, 3), int); depth = np.zeros(Nc, int); cen = np.zeros((Nc, 3))
for n, s in enumerate(sides):
    d = L / s
    for x in itertools.product(range(s), repeat=3):
        i = idx(n, x); coords[i] = x; depth[i] = n; cen[i] = (np.array(x) + 0.5) * d

def links():
    """liens (i <- j, type, donnees) : lateraux (dir) et verticaux (enfant <- parent)"""
    lat = []
    for n, s in enumerate(sides):
        d = L / s
        for x in itertools.product(range(s), repeat=3):
            i = idx(n, x)
            for j in range(3):
                y = list(x); y[j] += 1; jj = idx(n, y)
                wrapx = (j == 0 and x[0] == s - 1)
                lat.append((i, jj, j, d, wrapx))
    ver = []
    for n in range(1, N + 1):
        for x in itertools.product(range(sides[n]), repeat=3):
            c = idx(n, x); p = idx(n - 1, x if flat else tuple(v // 2 for v in x))
            ver.append((c, p))
    return lat, ver
LAT, VER = links()

def phase_lat(i, jj, j, d, wrapx, th):
    """U_{i <- jj} = exp(i int_{jj}^{i} A) ; A = (0, B X, th/L) ; jj = i + d e_j (eventuellement enroule)"""
    X, Y = cen[i][0], cen[i][1]
    if j == 0: ph = 0.0 + (B * L * Y if wrapx else 0.0)       # translation magnetique psi(X+L) = e^{i B L Y} psi(X)
    elif j == 1: ph = -B * X * d
    else: ph = -th / L * d
    return np.exp(1j * ph)
def phase_ver(c, p, th):
    """U_{c <- p} = exp(i int_p^c A), segment droit entre centres"""
    Xp, Yp, Zp = cen[p]; Xc, Yc, Zc = cen[c]
    return np.exp(1j * (B * 0.5 * (Xp + Xc) * (Yc - Yp) + th / L * (Zc - Zp)))

# tableaux de liens vectorises : phase = base * exp(i theta cz)
LINK = []
for jdir in range(3):
    sel = [l for l in LAT if l[2] == jdir]
    R = np.array([l[0] for l in sel]); C = np.array([l[1] for l in sel])
    base = np.array([phase_lat(*l, 0.0) for l in sel]); cz = np.array([(-l[3] / L if jdir == 2 else 0.0) for l in sel])
    xj = np.array([coords[l[0]][jdir] for l in sel]); dl = np.array([depth[l[0]] for l in sel])
    LINK.append((R, C, base, cz, xj, dl))
_wv = 1.0 if flat else 1 / np.sqrt(8)
VR = np.array([c for c, p in VER]); VC = np.array([p for c, p in VER])
VB = np.array([_wv * phase_ver(c, p, 0.0) for c, p in VER]); VZ = np.array([(cen[c][2] - cen[p][2]) / L for c, p in VER])

def tdim(xj, dl):
    if MOD not in ('ksD', 'ksU'): return np.ones(len(xj))
    rho = np.where(dl >= N - 1, PAR, RHOG) if MOD == 'ksD' else np.full(len(xj), PAR)
    return np.where(xj % 2 == 0, rho, 1.0)
def ops(th, deriv=False):
    """T_j (T psi)(x) = U psi(x + e_j) ; S (S psi)(enfant) = w U psi(parent) ; deriv : d/dtheta"""
    T = []
    for (R, C, base, cz, xj, dl) in LINK:
        v = base * np.exp(1j * th * cz) * tdim(xj, dl); v = 1j * cz * v if deriv else v
        T.append(sp.csr_matrix((v, (R, C)), shape=(Nc, Nc)))
    v = VB * np.exp(1j * th * VZ); v = 1j * VZ * v if deriv else v
    return T, sp.csr_matrix((v, (VR, VC)), shape=(Nc, Nc))

def plaquettes():
    """controle : holonomie de chaque plaquette laterale (x, y) (sens x puis y) = exp(i B d^2), a tous les niveaux de cote >= 3,
    plaquettes enroulees comprises ; plaquettes (y, z) a theta = 0 : holonomie 1"""
    T, S = ops(0.0); err = 0.0; err2 = 0.0
    for n, s in enumerate(sides):
        if s < 3: continue
        d = L / s
        for x in itertools.product(range(s), repeat=3):
            a = idx(n, x); b = idx(n, (x[0] + 1, x[1], x[2])); c = idx(n, (x[0] + 1, x[1] + 1, x[2])); e = idx(n, (x[0], x[1] + 1, x[2]))
            holo = np.conj(T[0][a, b]) * np.conj(T[1][b, c]) * T[0][e, c] * T[1][a, e]
            err = max(err, abs(np.angle(holo * np.exp(-1j * B * d * d))))
            f = idx(n, (x[0], x[1], x[2] + 1)); g = idx(n, (x[0], x[1] + 1, x[2] + 1))
            holo2 = np.conj(T[1][a, e]) * np.conj(T[2][e, g]) * T[1][f, g] * T[2][a, f]
            err2 = max(err2, abs(np.angle(holo2)))
    return max(err, err2)

def kron(A, g): return sp.kron(A, sp.csr_matrix(g), format='csr')
def H_bor(th, m, deriv=False):
    T, S = ops(th, deriv); H = sp.csr_matrix((4 * Nc, 4 * Nc), dtype=complex)
    for j in range(3):
        H = H + kron((T[j] - T[j].getH()) / 2j, GAM[j]) - (r / 2) * kron(T[j] + T[j].getH(), GAM[4])
    H = H + kron((S - S.getH()) / 2j, GAM[3]) - (r / 2) * kron(S + S.getH(), GAM[4]) + (0.0 if deriv else (m + 4 * r)) * kron(sp.identity(Nc), GAM[4])
    return H.tocsr()
x1, x2, x3 = coords[:, 0], coords[:, 1], coords[:, 2]
eta = [np.ones(Nc), (-1.0) ** x1, (-1.0) ** (x1 + x2)]; eta4 = (-1.0) ** (x1 + x2 + x3)
eps = (-1.0) ** (x1 + x2 + x3 + depth)
def H_ks(th, M, wall, deriv=False):
    T, S = ops(th, deriv); H = sp.csr_matrix((Nc, Nc), dtype=complex)
    for j in range(3):
        A = sp.diags(eta[j]) @ T[j]; H = H + (A - A.getH()) / 2j
    V = sp.diags(eta4) @ S; H = H + (V - V.getH()) / 2j               # eta4 evalue sur l'enfant (ligne de S)
    Mn = np.where(depth >= N - 1, -M, M) if wall else np.full(Nc, M)
    return (H + (0 if deriv else sp.diags(Mn * eps))).tocsr()

def H_of(th, deriv=False):
    if MOD == 'bor': return H_bor(th, PAR, deriv)
    if MOD in ('ksD', 'ksU'): return H_ks(th, 0.0, False, deriv)
    return H_ks(th, PAR, MOD == 'ksW', deriv)
ncomp = 4 if MOD == 'bor' else 1
dHdth = lambda th: H_of(th, True)

# ---------------------------------------------------------------- controles
print(f'=== V202 : rho_gros {RHOG} ; modele {MOD}, geometrie {GEOM}, N = {N} (L = {L}), parametre {PAR}, Phi = 1, {NTH} pas de theta ===', flush=True)
print(f'  cellules {Nc}, dimension {ncomp * Nc}')
print(f'  flux des plaquettes (x, y) = B d^2 : ecart max {plaquettes():.1e}')
H0 = H_of(0.3); print(f'  hermiticite : {abs(H0 - H0.getH()).max():.1e}')
Dn = (H_of(0.3 + 1e-6) - H_of(0.3 - 1e-6)) / 2e-6; print(f'  dH/dtheta analytique vs differences finies : {abs(Dn - dHdth(0.3)).max():.1e}')
H2p = H_of(2 * np.pi)
# equivalence de jauge H(2 pi) ~ H(0) : spectres egaux
def lowest(H, k=40):
    if H.shape[0] <= 6000:
        return np.linalg.eigvalsh(H.toarray())
    return np.sort(sla.eigsh(H, k=k, sigma=1e-4, which='LM', return_eigenvectors=False))
e0 = lowest(H_of(0.0)); e2 = lowest(H2p)
k0 = min(len(e0), len(e2)); print(f'  spectre(theta = 2 pi) = spectre(0) : ecart {np.abs(np.sort(e0)[:k0] - np.sort(e2)[:k0]).max():.1e}')

# ---------------------------------------------------------------- flot spectral
K = 20
EPS = 0.25 * 2 * np.pi / L          # quart d'espacement des niveaux chiraux (2 pi v / L, v ~ 1) : flot net exact a travers la fenetre
ths = np.linspace(0, 2 * np.pi, NTH, endpoint=False); dth = ths[1] - ths[0]
flow = np.zeros(N + 1); gross = np.zeros(N + 1); emin = []
for th in ths:
    H = H_of(th); D = dHdth(th)
    if K is None:
        E, V = np.linalg.eigh(H.toarray())
        sel = np.abs(E) < EPS; E, V = E[sel], V[:, sel]
    else:
        E, V = sla.eigsh(H, k=K, sigma=1e-4, which='LM'); sel = np.abs(E) < EPS; E, V = E[sel], V[:, sel]
        if np.sum(np.abs(E) < EPS) == K: print('  ATTENTION : fenetre saturee', flush=True)
    emin.append(np.min(np.abs(E)) if len(E) else np.nan)
    for a in range(len(E)):
        v = V[:, a]; slope = np.real(np.vdot(v, D @ v))
        w = np.bincount(depth, weights=np.sum(np.abs(v.reshape(Nc, ncomp)) ** 2, axis=1), minlength=N + 1)
        flow += dth / (2 * EPS) * slope * w; gross += dth / (2 * EPS) * abs(slope) * w
print(f'  fenetre eps = {EPS:.3f}')
print(f'\n  flot net par couche (n = 0 .. {N}) : [{" ".join(f"{f:+.3f}" for f in flow)}] ; total {flow.sum():+.3f}')
print(f'  flot brut par couche             : [{" ".join(f"{g:.3f}" for g in gross)}] ; total {gross.sum():.3f}')
print(f'  couche fine (n = N) : net {flow[N]:+.3f}, brut {gross[N]:.3f} ; couches n <= N-1 : net {flow[:N].sum():+.3f}, brut {gross[:N].sum():.3f}')
print(f'  |E| minimal sur les theta : {np.nanmin(emin):.2e} ; theta ou aucun etat sous eps : {int(np.sum(np.isnan(emin)))}')
print(f'  temps {time.time() - t0:.0f} s')
