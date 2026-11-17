"""Point 5 (V163) : transport de jauge par blocs (b) re-ecrit independamment, contre echantillonnage (a).
Convention : u_n(x, j) = transport x -> x + e_j au niveau n (se transforme en g(x+e) u g(x)^*) ; T_j[t, t+e_j] = conj(u(t, j)) ;
S[enfant, parent] = u(parent -> enfant) / sqrt(8).
Fond : flux Phi = 1 (jauge G2 : A = (-B Y, 0, theta/L)) + perturbation A_y = a cos(k X - phi0), k = 2 pi 2^j / L.
Usage : python3 rv_jauge.py N mesure [options]
  mesure = base (covariance, flux, racine, facteur de forme, test de robustesse de f_n = 0)
         | p4 NTH (couplage des etats de basse energie, plusieurs mesures)
         | flot NTH TR (flot spectral, fenetres eps et eps/2)"""
import sys, os, itertools, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_core import Geo, GAM, layer_weights
t0 = time.time()
N = int(sys.argv[1]); MES = sys.argv[2]
L = 2 ** N; PHI0 = np.pi / 4; m_, r_ = -1.0, 1.0
geo = Geo(N, False, 'G2')
B = geo.B

PDIR = 'yX'
def fine_links(th, a, k, extra=None, q=1.0):
    """u_N(x, j) a partir du continu : A = (-B Y, a cos(k X - phi0), th/L) [+ extra(X, Y) ajoute a A_x]
    PDIR change la perturbation : 'yX' (defaut, V162), 'xY' : A_x += a cos(kY - phi0), 'zX' : A_z += a cos(kX - phi0)"""
    s = L; d = 1.0; I, J, K = np.meshgrid(np.arange(s), np.arange(s), np.arange(s), indexing='ij')
    X, Y = I + 0.5, J + 0.5
    py = a * np.cos(k * X - PHI0) if PDIR == 'yX' else 0.0
    px = a * np.cos(k * Y - PHI0) if PDIR == 'xY' else 0.0
    pz = a * np.cos(k * X - PHI0) if PDIR == 'zX' else 0.0
    # x : int_{X}^{X+1} -B Y dX = -B Y ; enroulement : image x+e_x = (0, y, z) + L e_x, section G2 triviale en x
    ux = np.exp(1j * (-B * Y * d + px * d + (extra(X, Y) * d if extra is not None else 0.0)))
    # y : int A_y dY = a cos(kX - phi0) d ; enroulement en y : psi(r + L e_y) = exp(-i B L X) psi(r)
    uy = np.exp(1j * (py * d)) * np.where(J == s - 1, np.exp(1j * B * L * X), 1.0)
    uz = np.exp(1j * (th / L + pz) * d) * np.ones_like(X)
    return [ux ** q, uy ** q, uz ** q]
# NB : avec la section G2, u(x -> x + e_y) a l'enroulement vaut exp(i int A) * conj(g) ou psi(image) = g psi(stocke) :
# T[t, s] = U_{t<-s} = exp(i int_{s~}^{t} A) g(s) ; u(t -> s) = conj(T[t, s]) = exp(i int_t^{s~} A) conj(g). Verifie ci-dessous contre Geo.

def block_level(Uf, s):
    """niveau grossier (cote s) a partir du niveau fin (cote 2s) : renvoie liens projetes, moyennes brutes, transports internes P[c]"""
    S2 = 2 * s
    def uf(j, p): return Uf[j][p[0] % S2, p[1] % S2, p[2] % S2]
    P = {}; U = [np.zeros((s, s, s), complex) for _ in range(3)]; RAW = [np.zeros((s, s, s), complex) for _ in range(3)]
    for x in itertools.product(range(s), repeat=3):
        for c in itertools.product((0, 1), repeat=3):
            pos = [2 * x[0], 2 * x[1], 2 * x[2]]; val = 1.0 + 0j
            for j in range(3):                     # chemin x, puis y, puis z
                if c[j]: val *= uf(j, pos); pos[j] += 1
            P[(x, c)] = val
    for x in itertools.product(range(s), repeat=3):
        for j in range(3):
            xe = list(x); xe[j] = (xe[j] + 1) % s; xe = tuple(xe); acc = 0j
            for c in itertools.product((0, 1), repeat=3):
                p = [2 * x[0] + c[0], 2 * x[1] + c[1], 2 * x[2] + c[2]]
                a1 = uf(j, p); p2 = list(p); p2[j] += 1; a2 = uf(j, p2)
                acc += P[(x, c)] * a1 * a2 * np.conj(P[(xe, c)])
            RAW[j][x] = acc / 8
            U[j][x] = acc / abs(acc) if abs(acc) / 8 > 1e-6 else 0.0
    return U, RAW, P

def sampled_level(n, th, a, k, q=1.0):
    """transport (a) : integrale du continu entre centres au niveau n"""
    s = 2 ** n; d = L / s; I, J, K = np.meshgrid(np.arange(s), np.arange(s), np.arange(s), indexing='ij')
    X, Y = (I + 0.5) * d, (J + 0.5) * d
    ux = np.exp(1j * (-B * Y * d))
    uy = np.exp(1j * a * np.cos(k * X - PHI0) * d) * np.where(J == s - 1, np.exp(1j * B * L * X), 1.0)
    uz = np.exp(1j * th / L * d) * np.ones_like(X)
    return [ux ** q, uy ** q, uz ** q]

def build(th, a, k, TR, gfine=None, extra=None, q=1.0):
    U = [None] * (N + 1); Pin = [None] * N; RAW = [None] * (N + 1)
    U[N] = fine_links(th, a, k, extra, q)
    if gfine is not None:
        U[N] = [np.roll(gfine, -1, axis=j) * U[N][j] * np.conj(gfine) for j in range(3)]
    for n in range(N - 1, -1, -1):
        Ub, raw, P = block_level(U[n + 1], 2 ** n)
        U[n] = Ub if TR == 'b' else sampled_level(n, th, a, k, q); RAW[n] = raw; Pin[n] = P
    # matrices
    rows = [[], [], []]; cols = [[], [], []]; vals = [[], [], []]; vr, vc, vv = [], [], []
    for i, (n, x) in enumerate(geo.cells):
        s = 2 ** n
        for j in range(3):
            y = list(x); y[j] = (y[j] + 1) % s
            rows[j].append(i); cols[j].append(geo.index[(n, tuple(y))]); vals[j].append(np.conj(U[n][j][x]))
        if n > 0:
            px = tuple(v // 2 for v in x); c = tuple(v % 2 for v in x); p = geo.index[(n - 1, px)]
            if TR == 'b': u = Pin[n - 1][(px, c)]
            else:
                # (a) : segment droit parent -> enfant, A lineaire en Y pour A_x, perturbation integree exactement
                Xp, Yp, Zp = geo.cen[p]; Xc, Yc, Zc = geo.cen[i]; dX, dY = Xc - Xp, Yc - Yp
                ph = -B * 0.5 * (Yp + Yc) * dX + th / L * (Zc - Zp) + a * dY * (np.sin(k * Xc - PHI0) - np.sin(k * Xp - PHI0)) / (k * dX)
                u = np.exp(1j * q * ph)
            vr.append(i); vc.append(p); vv.append(u / np.sqrt(8))
    Nc = geo.Nc
    T = [sp.csr_matrix((vals[j], (rows[j], cols[j])), shape=(Nc, Nc)) for j in range(3)]
    S = sp.csr_matrix((vv, (vr, vc)), shape=(Nc, Nc))
    kr = lambda A, G: sp.kron(A, sp.csr_matrix(G), format='csr')
    H = sp.csr_matrix((4 * Nc, 4 * Nc), dtype=complex)
    for j in range(3):
        H = H + kr((T[j] - T[j].getH()) / 2j, GAM[j]) - (r_ / 2) * kr(T[j] + T[j].getH(), GAM[4])
    H = H + kr((S - S.getH()) / 2j, GAM[3]) - (r_ / 2) * kr(S + S.getH(), GAM[4]) + (m_ + 4 * r_) * kr(sp.identity(Nc), GAM[4])
    return H.tocsr(), U, RAW, T, S

def kj(j): return 2 * np.pi * 2 ** j / L

if MES == 'base':
    # 0) (a) a theta, a = 0 doit coincider avec Geo (jauge G2) : verifie les conventions de ce script
    Hb, U, RAW, T, S = build(0.4, 0.0, kj(1), 'a')
    Tg, Sg = geo.ops(0.4)
    print(f'(0) transport (a) contre Geo (jauge G2) : |T - T_geo| = {max(abs(T[j] - Tg[j]).max() for j in range(3)):.1e}, |S - S_geo| = {abs(S - Sg).max():.1e}')
    # 1) covariance de (b)
    rng = np.random.default_rng(3); gf = np.exp(2j * np.pi * rng.random((L, L, L)))
    for q in (1.0, 2.0):
        E1 = np.linalg.eigvalsh(build(0.3, 0.2, kj(1), 'b', q=q)[0].toarray())
        E2 = np.linalg.eigvalsh(build(0.3, 0.2, kj(1), 'b', gfine=gf, q=q)[0].toarray())
        print(f'(1) covariance (b), q = {q:.0f}, a = 0.2 : ecart spectral {np.abs(E1 - E2).max():.1e}')
    # 2) flux par niveau et moyenne brute a la racine / niveau 1
    H, U, RAW, T, S = build(0.3, 0.0, kj(1), 'b')
    for n in range(N, -1, -1):
        s = 2 ** n
        if s >= 2:
            pl = U[n][0] * np.roll(U[n][1], -1, axis=0) * np.conj(np.roll(U[n][0], -1, axis=1)) * np.conj(U[n][1])
            fl = np.angle(pl).sum() / (2 * np.pi) / s
        else: fl = np.nan
        rr = '' if n == N else ' ; |moyenne brute| x,y,z : min ' + ' '.join(f'{np.abs(RAW[n][j]).min():.4f}' for j in range(3))
        print(f'(2) niveau {n} : flux par tranche z = {fl:.4f}{rr}')
    # 3) facteur de forme : RMS par cellule (V162) et norme d'operateur de P_n V P_n (definition de V162 §2)
    print('(3) facteur de forme f_n / f_N, n = 0..N  [RMS par cellule | norme d operateur]')
    dep4 = np.repeat(geo.depth, 4)
    for TR in ('a', 'b'):
        for j in range(1, N + 1):
            h = 1e-6
            V = (build(0.3, h, kj(j), TR)[0] - build(0.3, -h, kj(j), TR)[0]) / (2 * h)
            rms, opn = [], []
            for n in range(N + 1):
                Pn = sp.diags((dep4 == n).astype(float)); Mn = (Pn @ V @ Pn).tocsr()
                rms.append(np.sqrt(abs(Mn.multiply(Mn.conj())).sum() / (4 * 8 ** n)))
                opn.append(np.max(np.abs(np.linalg.eigvalsh(Mn[dep4 == n][:, dep4 == n].toarray()))) if Mn.nnz else 0.0)
            rms, opn = np.array(rms), np.array(opn)
            print(f'    ({TR}) j = {j} (l = L/{2 ** j}) : RMS [' + ' '.join(f'{x / rms[N]:.3g}' for x in rms) + '] | op [' + ' '.join(f'{x / opn[N]:.3g}' for x in opn) + ']')
    # 4) robustesse de l'annulation exacte : fond non uniforme (bruit de flux sur le niveau fin) ; perturbation sur A_x
    print('(4) robustesse de f_n = 0 (transport b, j = 2) : max |P_n V P_n| par niveau')
    rng = np.random.default_rng(5)
    noise = 0.3 * rng.standard_normal((L, L))
    def extra_noise(X, Y):
        return noise[(X - 0.5).astype(int) % L, (Y - 0.5).astype(int) % L]
    for lab, ex in (('fond uniforme', None), ('fond + bruit A_x(X,Y) amplitude 0.3', extra_noise)):
        h = 1e-6
        V = (build(0.3, h, kj(2), 'b', extra=ex)[0] - build(0.3, -h, kj(2), 'b', extra=ex)[0]) / (2 * h)
        mx = []
        for n in range(N + 1):
            Pn = sp.diags((dep4 == n).astype(float)); Mn = (Pn @ V @ Pn).tocsr(); mx.append(abs(Mn).max() if Mn.nnz else 0.0)
        print(f'    {lab} : [' + ' '.join(f'{x:.2e}' for x in mx) + ']')
    print(f'  temps {time.time() - t0:.0f} s')

def eig_low(H, k=24):
    if H.shape[0] <= 2500 and k >= H.shape[0]:
        return np.linalg.eigh(H.toarray())
    E, V = sla.eigsh(H, k=k, sigma=1e-5, which='LM'); o = np.argsort(np.abs(E)); return E[o], V[:, o]

if MES == 'p4':
    NTH = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    EPS = 0.25 * 2 * np.pi / L; h = 1e-6
    ths = np.linspace(0, 2 * np.pi, NTH, endpoint=False)
    for j in ([1, 2] if N == 3 else [2, 3]):
        for TR in ('a', 'b'):
            recs = []
            for th in ths:
                H = build(th, 0.0, kj(j), TR)[0]
                V = (build(th, h, kj(j), TR)[0] - build(th, -h, kj(j), TR)[0]) / (2 * h)
                V0 = (build(th, h, 0.0, TR)[0] - build(th, -h, 0.0, TR)[0]) / (2 * h) if False else None
                E, W = eig_low(H, 24); win = np.abs(E) < EPS; E, W = E[win], W[:, win]
                if len(E) == 0: continue
                w = layer_weights(geo, W, 4); nb = np.arange(N + 1) @ w
                VW = V @ W; low = W.conj().T @ VW
                for i in range(len(E)):
                    recs.append((th, E[i], nb[i], np.linalg.norm(VW[:, i]), np.sum(np.abs(low[:, i]) ** 2), w[:, i]))
            nbs = np.array([x[2] for x in recs]); c1 = np.array([x[3] for x in recs]); c2 = np.array([x[4] for x in recs])
            def stat(mask):
                return (mask.sum(), np.median(c1[mask]) if mask.any() else np.nan, np.median(c2[mask]) if mask.any() else np.nan)
            e1 = stat(nbs >= N - 0.5); m1 = stat(nbs <= N - 1.5); m2 = stat(nbs < N - 0.5); e2 = stat(nbs >= N - 1.0); m3 = stat(nbs < N - 1.0)
            print(f'N={N} j={j} (l=L/{2**j}) ({TR}) {NTH} theta, {len(recs)} etats dans la fenetre ; n moyens : ' + ' '.join(f'{x:.2f}' for x in sorted(nbs)))
            print(f'   classes V162 : bord n>=N-0.5 : {e1[0]} etats, ||V psi|| med {e1[1]:.4f}, sum_fenetre |<psi\'|V|psi>|^2 med {e1[2]:.2e} ;'
                  f' miroir n<=N-1.5 : {m1[0]} etats, ||V psi|| med {m1[1]:.4f}, fenetre {m1[2]:.2e} ; rapports {m1[1]/e1[1]:.3f} | {m1[2]/e1[2] if e1[2] > 0 else np.nan:.3g}')
            print(f'   classes completes (n<N-0.5 = miroir) : {m2[0]} etats, ||V psi|| med {m2[1]:.4f} ; rapport {m2[1]/e1[1]:.3f} ;'
                  f' coupure N-1 : bord {e2[0]} ({e2[1]:.4f}), miroir {m3[0]} ({m3[1]:.4f}), rapport {m3[1]/e2[1]:.3f}')
            # poids de chaque classe par couche (moyen)
            for lab, mask in (('bord', nbs >= N - 0.5), ('miroir<=N-1.5', nbs <= N - 1.5), ('intermediaire', (nbs > N - 1.5) & (nbs < N - 0.5))):
                if mask.any():
                    wm = np.mean([recs[i][5] for i in np.where(mask)[0]], axis=0)
                    print(f'     profil moyen {lab} : [' + ' '.join(f'{x:.2f}' for x in wm) + ']  theta : ' + ' '.join(f'{recs[i][0]:.2f}' for i in np.where(mask)[0][:12]))
            print(f'   temps {time.time() - t0:.0f} s', flush=True)

if MES == 'flot':
    NTH = int(sys.argv[3]); TR = sys.argv[4]; q = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    EPS = 0.25 * 2 * np.pi / L; h = 1e-6
    ths = np.linspace(0, 2 * np.pi, NTH, endpoint=False); dth = ths[1] - ths[0]
    F = np.zeros(N + 1); Fh = np.zeros(N + 1)
    for th in ths:
        H = build(th, 0.0, kj(1), TR, q=q)[0]; D = (build(th + h, 0.0, kj(1), TR, q=q)[0] - build(th - h, 0.0, kj(1), TR, q=q)[0]) / (2 * h)
        E, W = eig_low(H, 30)
        if np.sum(np.abs(E) < EPS) >= 28: print('ATTENTION fenetre saturee')
        w = layer_weights(geo, W, 4); sl = np.real(np.einsum('ij,ij->j', W.conj(), D @ W))
        a1 = np.abs(E) < EPS; a2 = np.abs(E) < EPS / 2
        F += dth / (2 * EPS) * (w[:, a1] * sl[a1]).sum(1); Fh += dth / EPS * (w[:, a2] * sl[a2]).sum(1)
    print(f'flot ({TR}), q = {q:.0f}, N = {N}, {NTH} pas : fenetre eps [' + ' '.join(f'{x:+.3f}' for x in F) + f'] total {F.sum():+.3f} ;'
          f' fenetre eps/2 [' + ' '.join(f'{x:+.3f}' for x in Fh) + f'] total {Fh.sum():+.3f} ; temps {time.time() - t0:.0f} s', flush=True)

if MES == 'ffgen':
    # annulation exacte de f_n pour d'autres perturbations et un fond non uniforme (transport b seulement)
    dep4 = np.repeat(geo.depth, 4); rng = np.random.default_rng(5); noise = 0.3 * rng.standard_normal((L, L))
    def extra_noise(X, Y): return noise[(X - 0.5).astype(int) % L, (Y - 0.5).astype(int) % L]
    for pdir in ('yX', 'xY', 'zX'):
        PDIR = pdir
        for lab, ex in (('uniforme', None), ('bruit', extra_noise)):
            for j in range(1, N + 1):
                h = 1e-6
                V = (build(0.3, h, kj(j), 'b', extra=ex)[0] - build(0.3, -h, kj(j), 'b', extra=ex)[0]) / (2 * h)
                mx = []
                for n in range(N + 1):
                    Pn = sp.diags((dep4 == n).astype(float)); Mn = (Pn @ V @ Pn).tocsr(); mx.append(abs(Mn).max() if Mn.nnz else 0.0)
                mx = np.array(mx)
                print(f'  {pdir} fond {lab:8s} j = {j} : max|P_n V P_n| / niveau N = [' + ' '.join(f'{x / mx[N]:.1e}' for x in mx) + ']', flush=True)
    print(f'  temps {time.time() - t0:.0f} s')
