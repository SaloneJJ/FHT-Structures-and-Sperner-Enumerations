"""V148 (e) : le raccord (V148 §3).  Ea metrique induite ; Eb compression de Haar ; Ec doubleurs, Wilson, SLAC, Schur ; Ed comptage.
Criteres et predictions ecrits AVANT calcul dans V148 §3.   Usage : python3 V148_FHT_e_raccord.py [Ea Eb Ec Ed]
JJ SALONE, FHT 2026"""
import itertools, sys, time
import numpy as np
from scipy.special import erf
t00 = time.time()
np.set_printoptions(precision=4, suppress=True, linewidth=150)
parts = set(sys.argv[1:]) or {'Ea', 'Eb', 'Ec', 'Ed'}

# ================================================================== Ea : metrique induite
def Ea():
    print('=== Ea : fermeture par chaines de la quasi-distance dyadique sur la grille de profondeur n ===', flush=True)
    for n in (1, 2, 3, 4):
        side = 2 ** n + 1
        P = np.array(list(itertools.product(range(side), repeat=3)), dtype=np.int64)
        N = len(P)
        # m(a,b) = plus grande profondeur m' dans 0..n telle que a et b sont dans une meme cellule fermee de cote s = 2^(n-m')
        m = np.zeros((N, N), dtype=np.int8)
        for mp in range(n + 1):
            s = 2 ** (n - mp)
            lo = np.minimum(P[:, None, :], P[None, :, :]); hi = np.maximum(P[:, None, :], P[None, :, :])
            ok = np.all(-(-hi // s) - lo // s <= 1, axis=2)      # ceil(hi/s) - floor(lo/s) <= 1
            m = np.where(ok, mp, m)
        delta = (2.0 ** (-m.astype(np.float64))).astype(np.float64)
        np.fill_diagonal(delta, 0.0)
        # Dijkstra dense depuis des sources
        rng = np.random.default_rng(7)
        srcs = range(N) if n <= 3 else rng.choice(N, 12, replace=False)
        worst = 0.0; ratios = []
        for s0 in srcs:
            dist = np.full(N, np.inf); dist[s0] = 0.0; done = np.zeros(N, bool)
            for _ in range(N):
                u = np.argmin(np.where(done, np.inf, dist)); done[u] = True
                dist = np.minimum(dist, dist[u] + delta[u])
            inf_norm = np.abs(P - P[s0]).max(axis=1) / 2 ** n
            mask = inf_norm > 0
            r = dist[mask] / inf_norm[mask]; ratios.append((r.min(), r.max()))
            worst = max(worst, np.abs(dist - inf_norm).max())
        lo = min(a for a, b in ratios); hi = max(b for a, b in ratios)
        print(f'  n = {n} : {N} sommets, {"toutes les sources" if n<=3 else "12 sources"} ; rapport (delta_hat / ||x-y||_inf) dans [{lo:.12f}, {hi:.12f}] ; ecart absolu max {worst:.2e}', flush=True)
    print('  => delta_hat = ||x - y||_inf exactement (cube unite)  [E, n <= 4]  ; ||.||_inf <= ||.||_2 <= sqrt(3) ||.||_inf', flush=True)

# ================================================================== Eb : compression de Haar
def Eb():
    print('\n=== Eb : <u, d(eta_eps * v)> pour u = 1_[i h, (i+1) h), v = 1_[j h, (j+1) h), eta_eps gaussien symetrique ===', flush=True)
    h = 1.0
    def U(i, y, eps):   # int_{i h}^{(i+1) h} eta_eps(x - y) dx
        a, b = i * h, (i + 1) * h
        return 0.5 * (erf((b - y) / (eps * np.sqrt(2))) - erf((a - y) / (eps * np.sqrt(2))))
    for eps in (0.3, 0.03, 0.003):
        M = np.zeros((5, 5)); idx = range(-2, 3)
        for ii, i in enumerate(idx):
            for jj, j in enumerate(idx):
                M[ii, jj] = U(i, j * h, eps) - U(i, (j + 1) * h, eps)
        print(f'  eps = {eps} : ligne i = 0, colonnes j = -2..2 : {M[2]}  (attendu ... +1/2 en j = 1 et -1/2 en j = -1)')
    print('  => limite eps -> 0 : (D)_{ij} = (1/2)(delta_{i,j-1} - delta_{i,j+1}) : difference centrale (S - S^-1)/2 [E, convention de valeur moyenne = limite symetrique]', flush=True)

# ================================================================== modele F (copie de V146, meme graine)
rng = np.random.default_rng(146)
nA, nB, nG = 4, 4, 3
dI = nA * nB * nG; dF = 2 * dI
def iI(a, b, g): return (a * nB + b) * nG + g
gA = np.array([1, -1, -1, 1])
gamma_F = np.zeros(dF)
for a, b, g in itertools.product(range(nA), range(nB), range(nG)):
    gamma_F[iI(a, b, g)] = gA[a]; gamma_F[dI + iI(a, b, g)] = -gA[a]
GF = np.diag(gamma_F)
UF = np.zeros((dF, dF)); UF[:dI, dI:] = np.eye(dI); UF[dI:, :dI] = np.eye(dI)
def rep(lam, q, m):
    Pm = np.zeros((dF, dF), complex)
    MA = np.zeros((nA, nA), complex); MA[0, 0] = lam; MA[3, 3] = np.conj(lam); MA[1:3, 1:3] = q
    Pm[:dI, :dI] = np.kron(np.kron(MA, np.eye(nB)), np.eye(nG))
    MB = np.zeros((nB, nB), complex); MB[0, 0] = lam; MB[1:, 1:] = m
    Pm[dI:, dI:] = np.kron(np.kron(np.eye(nA), MB), np.eye(nG))
    return Pm
def alea_alg():
    lam = rng.normal() + 1j * rng.normal()
    al, be = rng.normal(size=2) + 1j * rng.normal(size=2)
    q = np.array([[al, be], [-np.conj(be), np.conj(al)]])
    m = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    return rep(lam, q, m)
def herm(n): return rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
Yu, Yd, Yn, Ye = herm(3), herm(3), herm(3), herm(3)
YR = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); YR = YR + YR.T
def build_DF(Yu, Yd, Yn, Ye, YR):
    S = np.zeros((dI, dI), complex)
    for b in range(nB):
        Yup, Ydn = (Yn, Ye) if b == 0 else (Yu, Yd)
        for g, h in itertools.product(range(nG), repeat=2):
            S[iI(1, b, g), iI(0, b, h)] += Yup[g, h]
            S[iI(2, b, g), iI(3, b, h)] += Ydn[g, h]
    S = S + S.conj().T
    T = np.zeros((dI, dI), complex)
    for g, h in itertools.product(range(nG), repeat=2): T[iI(0, 0, g), iI(0, 0, h)] = YR[g, h]
    D = np.zeros((dF, dF), complex)
    D[:dI, :dI] = S; D[dI:, dI:] = S.conj(); D[dI:, :dI] = T; D[:dI, dI:] = T.conj().T
    return D
DF = build_DF(Yu, Yd, Yn, Ye, YR)
Z3 = np.zeros((3, 3), complex)
DF0 = build_DF(Yu, Yd, Z3, Ye, Z3)          # variante : neutrinos sans masse (Yn = 0, pas de Majorana)
def Japp(U, X): return U @ X.conj() @ U.T
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
SIG = [sx, sy, sz]; iSy = 1j * sy
I2 = np.eye(2)
U = np.kron(UF @ GF, iSy).real                 # J = U o conj o (k -> -k)

# --- symboles de derivee et de Wilson sur le tore a N cellules par cote (unites de reseau, h = 1)
def symbols(N, kind, r=0.0):
    ms = np.arange(-N // 2, N // 2)
    p = 2 * np.pi * ms / N
    P1, P2, P3 = np.meshgrid(p, p, p, indexing='ij')
    PP = [P1.ravel(), P2.ravel(), P3.ravel()]
    if kind == 'naif':   s = [np.sin(x) for x in PP]
    elif kind == 'slac': s = [np.where(np.isclose(np.abs(x), np.pi), 0.0, x) for x in PP]   # branche principale, mode de Nyquist -> 0
    w = r * sum(2 * (1 - np.cos(x)) for x in PP)
    return PP, np.array(s).T, w                 # (3 x N^3 impulsions), s (N^3, 3), w (N^3)

def Dk(s, w, DFm, mode):
    """bloc de moment (192 x 192) : D_F (x) 1 [(x) (1 + w) si Wilson-Yukawa] + gamma_F (x) (sigma.s [+ w])"""
    DS = s[0] * sx + s[1] * sy + s[2] * sz
    if mode == 'W':   return np.kron(DFm, I2) + np.kron(GF, DS + w * I2)
    if mode == 'Y':   return np.kron(DFm * (1 + w), I2) + np.kron(GF, DS)
    return np.kron(DFm, I2) + np.kron(GF, DS)

# ================================================================== Ec
def Ec():
    N = 8
    print(f'\n=== Ec : D_n sur H_F (x) C^2 (x) l^2(T^3, {N}^3 cellules), blocs de moment (192 x 192), unites de reseau ===', flush=True)
    PP, s_nai, _ = symbols(N, 'naif'); _, s_sl, _ = symbols(N, 'slac')
    nk = N ** 3
    # indice de -k
    ms = np.arange(-N // 2, N // 2)
    def minusk(i):
        a = np.unravel_index(i, (N, N, N)); b = tuple((-(x - N // 2)) % N for x in a)   # m -> -m  (m = x - N/2)
        b = tuple(((-(x - N // 2)) + N // 2) % N for x in a)
        return np.ravel_multi_index(b, (N, N, N))
    # ---------- Ec2 : relations KO-1 (sur 40 impulsions)
    def ko_err(s, w, mode, DFm=DF, nsamp=40):
        errs = []; ids = np.random.default_rng(1).choice(nk, nsamp, replace=False)
        for i in ids:
            j = minusk(i)
            Di = Dk(s[i], w[i], DFm, mode); Dj = Dk(s[j], w[j], DFm, mode)
            errs.append(np.abs(U @ Di.conj() @ U.T + Dj).max())
            assert np.allclose(Di, Di.conj().T)
        return max(errs)
    r = 0.7
    print('  Ec2 : ecart max de JD + DJ (sur 40 impulsions) et hermiticite, J^2 = +1 :', f'J^2 = {np.allclose(U @ U, np.eye(2*dF))}')
    for nom, sk, mode in (('naif (N)', 'naif', 'N'), ('SLAC (S)', 'slac', 'N'), ('Wilson scalaire (W), r=0.7', 'naif', 'W'), ('Wilson-Yukawa (Y), r=0.7', 'naif', 'Y'), ('SLAC + Wilson-Yukawa', 'slac', 'Y')):
        _, s, w = symbols(N, sk, r)
        print(f'      {nom:30s} : {ko_err(s, w, mode):.1e}')
    # ---------- ordre un pour A_F : [[D, a (x) 1], b0 (x) 1] avec a, b dans A_F (agit sur F seul)
    A2 = np.kron(alea_alg(), I2); B2 = Japp(np.kron(UF, np.eye(2)), np.kron(alea_alg(), I2))
    worst = 0.0
    for nom, sk, mode in (('naif', 'naif', 'N'), ('SLAC', 'slac', 'N'), ('W', 'naif', 'W'), ('Y', 'naif', 'Y')):
        _, s, w = symbols(N, sk, r)
        for i in np.random.default_rng(2).choice(nk, 10, replace=False):
            Di = Dk(s[i], w[i], DF, mode); C = Di @ A2 - A2 @ Di
            worst = max(worst, np.abs(C @ B2 - B2 @ C).max())
    print(f'  ordre un pour A_F (a (x) 1 avec F seul), toutes variantes : ecart max {worst:.1e}')
    # ---------- Ec1, Ec3, Ec5 : comptage des modes bas dans le secteur D_F = 0 et avec masses
    Ec_cut = 0.5
    ball = (4 * np.pi / 3) * Ec_cut ** 3 * (N / (2 * np.pi)) ** 3        # nombre d'impulsions du continu avec |p| < Ec_cut
    cont = 2 * dF / 2 * 0 + dF * ball                                    # chaque etat de F donne w +- |s| : 1 etat |E|<cut par (etat F, signe) ; ici dF etats
    print(f'\n  Comptage des valeurs propres |E| < {Ec_cut} (secteur D_F = 0, 96 etats de F x 2 composantes de spin) ; continu attendu : {96*2*ball/ (1):.0f}')
    # pour D_F = 0 : valeurs propres gamma_F * (w +- |s|), chacune avec multiplicite 96 (etats de F) :
    def n_low_free(s, w, cut):
        a = np.abs(w + np.linalg.norm(s, axis=1)); b = np.abs(w - np.linalg.norm(s, axis=1))
        return 96 * (np.sum(a < cut) + np.sum(b < cut))
    ref = None
    for nom, sk, rr in (('naif (N)', 'naif', 0.0), ('SLAC (S)', 'slac', 0.0), ('Wilson scalaire r=0.05', 'naif', 0.05), ('Wilson scalaire r=0.2', 'naif', 0.2), ('Wilson scalaire r=0.5', 'naif', 0.5), ('Wilson scalaire r=1', 'naif', 1.0)):
        _, s, w = symbols(N, sk, rr)
        c = n_low_free(s, w, Ec_cut)
        ref = ref or c
        print(f'      {nom:26s} : {c:7d} modes ; rapport au continu {c/(96*2*ball):6.2f} ; min |E| hors k=0 : {min(np.sort(np.abs(w+np.linalg.norm(s,axis=1)))[1:3].min(), np.sort(np.abs(w-np.linalg.norm(s,axis=1)))[1:3].min()):.3f}')
    print('  Wilson scalaire : nappe w = |s| (modes d\'energie nulle) sur un reseau plus fin, N = 32 :')
    for rr in (0.05, 0.2, 0.5, 1.0):
        _, s, w = symbols(32, 'naif', rr)
        g = np.abs(w - np.linalg.norm(s, axis=1))
        print(f'      r = {rr:4.2f} : min |w - |s|| sur k != 0 : {np.sort(g)[1]:.2e} ; nombre de k avec |w - |s|| < 0.05 : {np.sum(g<0.05)} (k=0 compris)')
    # Ec5 : Wilson-Yukawa : doubleurs des especes de masse nulle (variante DF0 : nu sans masse)
    print('\n  Wilson-Yukawa (Y) : valeurs propres |E|^2 = mu^2 (1+w)^2 + |s|^2 ; especes de masse nulle (variante D_F : Yn = 0, pas de Majorana)')
    mu2 = np.sort(np.linalg.eigvalsh(DF0 @ DF0))
    n0 = int(np.sum(mu2 < 1e-10)); print(f'      etats de masse nulle dans la variante : {n0} sur {dF}')
    for rr in (0.0, 0.5):
        _, s, w = symbols(N, 'naif', rr)
        cnt = 0
        for i in range(nk):
            # D^2 = D_F^2 (1+w)^2 + |s|^2 : seules les especes de masse nulle ont E = |s|
            E = np.sqrt(mu2[mu2 < 1e-10] * 0 + np.dot(s[i], s[i]))
            cnt += np.sum(E < Ec_cut)
        print(f'      r = {rr} : modes |E| < {Ec_cut} des especes sans masse (2 spins) : {2*cnt:d}  (continu : {2*n0*ball:.0f}) ; rapport {2*cnt/(2*n0*ball):.2f}')
    print('      => masse des doubleurs massifs = mu (1 + w_doubler), avec w_doubler =', (symbols(N, "naif", 0.5)[2].max()).round(2), 'pour r = 0.5 : un doubleur d\'espece massive est lourd si mu (1+w) >> 1/N seulement')
    # ---------- Ec4 : Schur
    print('\n  Ec4 (Schur) : dimension reelle de {B : B^+ = B, B gamma_F = -gamma_F B, J_F B J_F^-1 = B, [B, pi(a)] = 0, [B, J pi(b) J^-1] = 0}')
    X = []
    for _ in range(3):
        a = alea_alg(); X.append(a + a.conj().T); X.append(Japp(UF, a) + Japp(UF, a).conj().T)
    # deux elements hermitiens generiques de l'algebre engendree par pi(A_F) et J pi(A_F) J^-1
    c1, c2 = rng.normal(size=6), rng.normal(size=6)
    X1 = sum(c * x for c, x in zip(c1, X)); X2 = sum(c * x for c, x in zip(c2, X))
    # commutant de X1 : base par les sous-espaces propres
    w1, V1 = np.linalg.eigh(X1)
    groups = []; cur = [0]
    for i in range(1, dF):
        if abs(w1[i] - w1[cur[0]]) < 1e-8: cur.append(i)
        else: groups.append(cur); cur = [i]
    groups.append(cur)
    unknowns = [(g, a_, b_) for g in groups for a_ in g for b_ in g]
    nu = len(unknowns)
    print(f'      sous-espaces propres de X1 : {len(groups)}, inconnues complexes : {nu}')
    # base de B : E_{ab} dans la base propre ; contraintes lineaires reelles sur (Re, Im)
    basis = []
    for (g, a_, b_) in unknowns:
        E = np.zeros((dF, dF), complex); E[a_, b_] = 1; basis.append(V1 @ E @ V1.conj().T)
    cols = []
    for Bm in basis:
        for coef in (1.0, 1.0j):
            Bc = coef * Bm
            cons = np.concatenate([
                (Bc - Bc.conj().T).ravel(),                      # hermiticite
                (Bc @ GF + GF @ Bc).ravel(),                     # impair
                (Japp(UF, Bc) - Bc).ravel(),                     # J_F B J_F^-1 = B
                (Bc @ X2 - X2 @ Bc).ravel(),
            ])
            cols.append(np.concatenate([cons.real, cons.imag]))
    Mc = np.array(cols).T
    sv = np.linalg.svd(Mc, compute_uv=False)
    dimnull = int(np.sum(sv < 1e-9 * sv.max())) + (Mc.shape[1] - len(sv))
    print(f'      inconnues reelles {Mc.shape[1]}, rang {int(np.sum(sv > 1e-9 * sv.max()))}, DIMENSION DU NOYAU = {dimnull}')
    # temoin : sans la condition \"impair\" (doit etre > 0) pour verifier que le calcul voit bien des solutions
    cols2 = []
    for Bm in basis:
        for coef in (1.0, 1.0j):
            Bc = coef * Bm
            cons = np.concatenate([(Bc - Bc.conj().T).ravel(), (Japp(UF, Bc) - Bc).ravel(), (Bc @ X2 - X2 @ Bc).ravel()])
            cols2.append(np.concatenate([cons.real, cons.imag]))
    M2 = np.array(cols2).T; sv2 = np.linalg.svd(M2, compute_uv=False)
    print(f'      temoin (sans la condition impair) : dimension du noyau = {int(np.sum(sv2 < 1e-9 * sv2.max())) + (M2.shape[1] - len(sv2))}  (operateurs pairs, hermitiens, J-invariants, dans le commutant)')
    # ---------- Ec6 : SLAC, localite et commutateurs en 1D
    print('\n  Ec6 : commutateur [D_n, f] en 1D (D_n = (i)-derivee ; f multiplication), unites physiques x in [0,1], D_phys = N * D_reseau')
    def comm_norm(N, f, kind):
        x = (np.arange(N) + 0.5) / N
        m = np.arange(-N // 2, N // 2)
        p = 2 * np.pi * m / N
        sd = np.sin(p) if kind == 'naif' else np.where(np.isclose(np.abs(p), np.pi), 0.0, p)
        F = np.fft.fft(np.eye(N), axis=0) / np.sqrt(N)
        # Fourier : ordre standard de numpy
        freqs = np.fft.fftfreq(N, d=1.0 / N)                 # entiers -N/2..N/2-1 en ordre numpy
        pp = 2 * np.pi * freqs / N
        sd = np.sin(pp) if kind == 'naif' else np.where(np.isclose(np.abs(pp), np.pi), 0.0, pp)
        Dm = np.fft.ifft(np.diag(sd)[None, :, :].squeeze() @ np.fft.fft(np.eye(N), axis=0), axis=0)
        Dm = (np.fft.ifft(np.eye(N), axis=0) @ np.diag(sd) @ np.fft.fft(np.eye(N), axis=0))
        Dm = N * Dm
        Fm = np.diag(f(x))
        return np.linalg.norm(Dm @ Fm - Fm @ Dm, 2)
    fcos = lambda x: np.cos(2 * np.pi * x); ftri = lambda x: np.abs(2 * x - 1) * 1.0; fsaw = lambda x: np.sin(2 * np.pi * x) ** 1
    print('      N      : ' + '  '.join(f'{N:7d}' for N in (16, 32, 64, 128, 256)))
    for nom, f, lipf in (('cos 2 pi x (Lip = 2 pi)', fcos, 2 * np.pi), ('onde triangulaire |2x-1| (Lip = 2)', ftri, 2.0)):
        for kind in ('naif', 'slac'):
            vals = [comm_norm(N, f, kind) for N in (16, 32, 64, 128, 256)]
            print(f'      {nom:36s} {kind:5s}: ' + '  '.join(f'{v:7.3f}' for v in vals) + f'   (Lip = {lipf:.3f})')

def Ed():
    print('\n=== Ed : comptage ===')
    print(f'  H_F (x) S : {dF} x 2 = {2*dF} composantes complexes par cellule ; Majorana J psi = psi : {2*dF} reels = {dF} complexes equivalents ; tuile de K26 : {dF} etats ; 48 Weyl x 2 = {48*2}')
    print('  coincidence de comptage : ', dF == 48 * 2 == 2 * dF // 2)

for p, fn in (('Ea', Ea), ('Eb', Eb), ('Ec', Ec), ('Ed', Ed)):
    if p in parts: fn()
print(f'\ntemps {time.time()-t00:.0f} s')
