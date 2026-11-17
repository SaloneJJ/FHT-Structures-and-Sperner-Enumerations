"""V146 : chantier 2 (V144 §5) -- triplet produit (bord Sigma = T3) x (partie finie F, KO 6), impair x pair.
H = H_F (x) L2(Sigma, S) ; D = D_F (x) 1 + gamma_F (x) D_Sigma ; J = (J_F gamma_F) (x) J_Sigma ; pas de gamma global.
F : modele de meme structure que le triplet fini de K26 (articles V-VII) : H_F = I + II (48 + 48),
    I = A (x) B (x) generations, A = (s_a, d1, d2, s_b), gamma_A = (+,-,-,+), B = (lepton, 3 couleurs), 3 generations ;
    A_F = C + H + M3(C) agit sur l'indice A dans I (lambda, q, conj(lambda)) et sur l'indice B dans II (lambda, m) ;
    J_F echange I et II avec conjugaison complexe ; D_F = Yukawa (s_a<->d1, s_b<->d2, matrices 3x3 aleatoires) + Majorana nu_R.
Sigma : tore T3 = [0, 2 pi]^3, structure spin eps in {0,1}^3 (k = n + eps/2), D_Sigma(k) = sigma.k, J_Sigma = i sigma_y o conj (k -> -k).
Verifications (V144 §5.2) : (a) KO-dimension 1 ; (b) ordre zero et ordre un (exacts mode par mode) ; (c) spectre
pythagoricien, cone commun, noyau ; (d) chiralite -> helicite ; (e) antisymetrie de <J psi, D phi> (Majorana) ;
(f) controles (sans gamma_F ; J = J_F (x) J_Sigma) ; loi de Weyl ; structures spin. JJ SALONE, FHT 2026"""
import numpy as np, itertools, sys
rng = np.random.default_rng(146)
np.set_printoptions(precision=4, suppress=True)
# ------------------------------------------------------------------ partie finie F
nA, nB, nG = 4, 4, 3
dI = nA * nB * nG; dF = 2 * dI                          # 48 + 48 = 96
def iI(a, b, g): return (a * nB + b) * nG + g          # indice dans I (et dans II, decale de dI)
gA = np.array([1, -1, -1, 1])
gamma_F = np.zeros(dF)
for a, b, g in itertools.product(range(nA), range(nB), range(nG)):
    gamma_F[iI(a, b, g)] = gA[a]; gamma_F[dI + iI(a, b, g)] = -gA[a]
GF = np.diag(gamma_F)
UF = np.zeros((dF, dF)); UF[:dI, dI:] = np.eye(dI); UF[dI:, :dI] = np.eye(dI)   # J_F = UF o conj
def rep(lam, q, m):
    """pi(lambda, q, m) sur H_F"""
    P = np.zeros((dF, dF), complex)
    MA = np.zeros((nA, nA), complex); MA[0, 0] = lam; MA[3, 3] = np.conj(lam); MA[1:3, 1:3] = q
    P[:dI, :dI] = np.kron(np.kron(MA, np.eye(nB)), np.eye(nG))
    MB = np.zeros((nB, nB), complex); MB[0, 0] = lam; MB[1:, 1:] = m
    P[dI:, dI:] = np.kron(np.kron(np.eye(nA), MB), np.eye(nG))
    return P
def alea_alg():
    lam = rng.normal() + 1j * rng.normal()
    al, be = rng.normal(size=2) + 1j * rng.normal(size=2)
    q = np.array([[al, be], [-np.conj(be), np.conj(al)]])
    m = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    return rep(lam, q, m)
def herm(n): X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); return X
Yu, Yd, Yn, Ye = herm(3), herm(3), herm(3), herm(3)
YR = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); YR = YR + YR.T          # Majorana, symetrique
S = np.zeros((dI, dI), complex)
for b in range(nB):
    Yup, Ydn = (Yn, Ye) if b == 0 else (Yu, Yd)
    for g, h in itertools.product(range(nG), repeat=2):
        S[iI(1, b, g), iI(0, b, h)] += Yup[g, h]          # d1 <- s_a
        S[iI(2, b, g), iI(3, b, h)] += Ydn[g, h]          # d2 <- s_b
S = S + S.conj().T
T = np.zeros((dI, dI), complex)                           # I -> II : nu_R -> conj(nu_R)
for g, h in itertools.product(range(nG), repeat=2): T[iI(0, 0, g), iI(0, 0, h)] = YR[g, h]
DF = np.zeros((dF, dF), complex)
DF[:dI, :dI] = S; DF[dI:, dI:] = S.conj(); DF[dI:, :dI] = T; DF[:dI, dI:] = T.conj().T
def Japp(U, X): return U @ X.conj() @ U.T                  # J X J^{-1} pour J = U o conj (U reel orthogonal)
print('=== Partie finie F (modele de structure K26 / CCM) : dim', dF)
print(f'  J_F^2 = 1 : {np.allclose(UF @ UF.conj(), np.eye(dF))} ; J_F gamma_F = -gamma_F J_F : {np.allclose(UF @ GF, -GF @ UF)}')
print(f'  D_F hermitien : {np.allclose(DF, DF.conj().T)} ; J_F D_F = D_F J_F : {np.allclose(Japp(UF, DF), DF)} ; D_F gamma_F = -gamma_F D_F : {np.allclose(DF @ GF, -GF @ DF)}')
e0 = e1 = 0.0
for _ in range(20):
    a, b = alea_alg(), alea_alg(); bo = Japp(UF, b)
    e0 = max(e0, np.abs(a @ bo - bo @ a).max()); C1 = DF @ a - a @ DF
    e1 = max(e1, np.abs(C1 @ bo - bo @ C1).max())
    assert np.allclose(a @ GF, GF @ a)
print(f'  ordre zero F : {e0:.1e} ; ordre un F : {e1:.1e} ; gamma_F commute a A_F : True')
# ------------------------------------------------------------------ espace Sigma = T3
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
SIG = [sx, sy, sz]; iSy = 1j * sy                          # reel : [[0,1],[-1,0]]
def DS(k): return k[0] * sx + k[1] * sy + k[2] * sz
def modes(N, eps):
    axes = [list(range(-N, N + 1)) if e == 0 else [x + 0.5 for x in range(-N - 1, N + 1)] for e in eps]
    return [tuple(float(c) for c in m) for m in itertools.product(*axes)]   # boite symetrique par k -> -k
I2 = np.eye(2)
def Dtot(k, DFm=DF, avec_gamma=True):
    return np.kron(DFm, I2) + (np.kron(GF, DS(k)) if avec_gamma else np.kron(np.eye(dF), DS(k)))
U = np.kron(UF @ GF, iSy).real                             # J = U o conj o (k -> -k)
print('\n=== Produit F x T3 (impair x pair) ===')
print(f'  U reel orthogonal : {np.allclose(U @ U.T, np.eye(2*dF))} ; J^2 = U U* = {np.round(np.diag(U @ U.conj())[:3].real, 6)} (attendu +1)')
ks = modes(2, (0, 0, 0)) + modes(2, (1, 0, 1))
errJD = max(np.abs(U @ Dtot(k).conj() @ U.T + Dtot(tuple(-x for x in k))).max() for k in ks)
errJDp = max(np.abs(U @ Dtot(k).conj() @ U.T - Dtot(tuple(-x for x in k))).max() for k in ks)
print(f'  (a) JD = -DJ : ecart {errJD:.1e} ; (JD = +DJ : ecart {errJDp:.2f}) => KO-dimension 1 (J^2 = +1, eps\' = -1)')
# (b) ordre zero et ordre un, exacts mode par mode : a = e^{i q.x} (x) a_F, b = e^{i q'.x} (x) b_F
e0 = e1 = 0.0
for _ in range(12):
    q = rng.integers(-1, 2, 3); qp = rng.integers(-1, 2, 3)
    aF, bF = alea_alg(), alea_alg()
    A2 = np.kron(aF, I2)
    Bo = U @ np.kron(bF, I2).conj() @ U.T                  # J b J^{-1} : multiplication par e^{-i q'.x} (k -> k - q')
    for k in modes(1, (1, 0, 0)):
        k = np.array(k)
        Ra = lambda kap: Dtot(kap + q) @ A2 - A2 @ Dtot(kap)   # [D, a] de kap vers kap + q
        lhs = Ra(k - qp) @ Bo - Bo @ Ra(k)
        e1 = max(e1, np.abs(lhs).max()); e0 = max(e0, np.abs(A2 @ Bo - Bo @ A2).max())
print(f'  (b) ordre zero : {e0:.1e} ; ordre un : {e1:.1e} (12 paires (a, b) aleatoires, decalages q, q\' dans {{-1,0,1}}^3)')
# (c) spectre pythagoricien
mu2 = np.sort(np.linalg.eigvalsh(DF @ DF)); err = 0.0
for k in ks[:40]:
    ev = np.sort(np.linalg.eigvalsh(Dtot(k)))
    kk = np.dot(k, k)
    err = max(err, np.abs(np.sort(np.abs(ev)) - np.sort(np.repeat(np.sqrt(mu2 + kk), 2))).max())
print(f'  (c) |valeurs propres| = sqrt(mu^2 + |k|^2) pour toutes les especes : ecart max {err:.1e} (meme cone, vitesse 1 pour toutes)')
# (d) chiralite -> helicite (D_F = 0, puis massif a grand |k|)
print('  (d) helicite des etats d\'energie positive (D_F = 0), par secteur :')
k = np.array([0.3, -0.7, 1.1]); kh = k / np.linalg.norm(k)
H = Dtot(k, DFm=np.zeros_like(DF)); w, Vv = np.linalg.eigh(H)
hel = np.kron(np.eye(dF), DS(kh))
noms = ['s_a (singulet)', 'd1 (doublet)', 'd2 (doublet)', 's_b (singulet)']
for sect, off in (('I', 0), ('II', dI)):
    for a in range(nA):
        proj = np.zeros(dF);
        for b, g in itertools.product(range(nB), range(nG)): proj[off + iI(a, b, g)] = 1
        P = np.kron(np.diag(proj), I2)
        pos = Vv[:, w > 1e-9]
        Pp = pos.conj().T @ P @ pos
        # helicite moyenne sur le sous-espace d'energie positive porte par ce secteur
        ww, vv = np.linalg.eigh(Pp); sub = pos @ vv[:, ww > 0.5]
        h = np.real(np.diag(sub.conj().T @ hel @ sub))
        print(f'      {sect} {noms[a]:16s} gamma_F = {int(gamma_F[off + iI(a,0,0)]):+d} : {sub.shape[1]} etats, helicite {np.round(np.unique(np.round(h, 6)), 4)}')
for kn in (0.5, 5, 50, 500):
    kk = kn * kh; Hm = Dtot(kk); w, Vv = np.linalg.eigh(Hm); pos = Vv[:, w > 0]
    proj = np.zeros(dF)
    for b, g in itertools.product(range(nB), range(nG)): proj[iI(1, b, g)] = 1
    P = np.kron(np.diag(proj), I2); hel = np.kron(np.eye(dF), DS(kh))
    num = np.real(np.trace(pos.conj().T @ P @ hel @ P @ pos)); den = np.real(np.trace(pos.conj().T @ P @ pos))
    print(f'      massif, |k| = {kn:6.1f} : helicite moyenne du poids d1 (doublet, I) dans les etats d\'energie positive = {num/den:+.4f}')
# (e) antisymetrie de A(psi, phi) = <J psi, D phi> sur le tore tronque (boite symetrique)
eps = (0, 1, 0); Ms = modes(1, eps); mi = {m: i for i, m in enumerate(Ms)}
d = 2 * dF
psi = rng.normal(size=(len(Ms), d)) + 1j * rng.normal(size=(len(Ms), d))
phi = rng.normal(size=(len(Ms), d)) + 1j * rng.normal(size=(len(Ms), d))
def Jv(x): return np.array([U @ x[mi[tuple(-c for c in m)]].conj() for m in Ms])
def Dv(x): return np.array([Dtot(np.array(m)) @ x[i] for i, m in enumerate(Ms)])
def form(x, y): return np.vdot(Jv(x).ravel(), Dv(y).ravel())
a1, a2 = form(psi, phi), form(phi, psi)
print(f'  (e) <J psi, D phi> = {a1:.4f} ; <J phi, D psi> = {a2:.4f} ; somme {abs(a1 + a2):.1e} => antisymetrique (eps eps\' = -1)')
print(f'      J^2 = +1 : la condition de Majorana J psi = psi definit un sous-espace reel de dimension reelle {d} par point (= dim_C H),'
      f' soit {d//2} composantes complexes = 48 champs de Weyl x 2 ; pas de quadruplement ni de projection de chiralite')
# (f) controles
errA = max(np.abs(U @ Dtot(k, avec_gamma=False).conj() @ U.T + Dtot(tuple(-x for x in k), avec_gamma=False)).max() for k in ks)
errB = max(np.abs(U @ Dtot(k, avec_gamma=False).conj() @ U.T - Dtot(tuple(-x for x in k), avec_gamma=False)).max() for k in ks)
k = np.array([0.4, 0.2, -0.9]); ev = np.sort(np.abs(np.linalg.eigvalsh(Dtot(k, avec_gamma=False))))
pyth = np.sort(np.repeat(np.sqrt(mu2 + k @ k), 2))
print(f'  (f) sans gamma_F : JD+DJ ecart {errA:.2f}, JD-DJ ecart {errB:.2f} (aucun signe) ; ecart au spectre pythagoricien {np.abs(ev - pyth).max():.2f}')
U2 = np.kron(UF, iSy).real
print(f'      J = J_F (x) J_Sigma : J^2 = {np.round(np.diag(U2 @ U2.conj())[0].real, 6)} ; '
      f'JD+DJ ecart {max(np.abs(U2 @ Dtot(k).conj() @ U2.T + Dtot(tuple(-x for x in k))).max() for k in ks):.2f}, '
      f'JD-DJ ecart {max(np.abs(U2 @ Dtot(k).conj() @ U2.T - Dtot(tuple(-x for x in k))).max() for k in ks):.2f} (signes meles)')
# loi de Weyl (D_F = 0) et structures spin
print('  Loi de Weyl de D_Sigma (une espece, 2 composantes) : N(Lambda) / [2 (4 pi / 3) Lambda^3] ->')
for Lam in (4, 8, 16):
    N = 0
    rng_n = range(-Lam - 1, Lam + 2)
    nn = np.array(list(itertools.product(rng_n, repeat=3)), float)
    N = 2 * int((np.linalg.norm(nn, axis=1) <= Lam).sum())
    print(f'      Lambda = {Lam:3d} : {N / (2 * 4 * np.pi / 3 * Lam**3):.4f}')
print('  Structures spin de T3 (eps in {0,1}^3) : dimension du noyau de D_Sigma (une espece)')
for eps in itertools.product((0, 1), repeat=3):
    ker = sum(2 for m in modes(1, eps) if np.allclose(m, 0))
    print(f'      eps = {eps} : noyau {ker}')
