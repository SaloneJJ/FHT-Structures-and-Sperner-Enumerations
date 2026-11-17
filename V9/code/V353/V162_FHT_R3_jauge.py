"""V162 : R3 -- jauge U(1) de fond sur l'arbre augmente (piste B, Gamma empruntes de V151) ; transport vers les niveaux grossiers
par (a) echantillonnage du continu (V158) ou (b) blocs (moyenne sur les 8 enfants, projection U(1), jauge du parent = enfant 0).
Criteres et predictions : V162 (ecrits avant calcul).
Usage : python3 V162_FHT_R3_jauge.py N m transport(a|b) q mesure(flot|ff|p4) [j] [NTH]
JJ SALONE, FHT 2026"""
import sys, time, itertools
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
t0 = time.time()
N, m, TR, q, MES = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], float(sys.argv[4]), sys.argv[5]
JW = int(sys.argv[6]) if len(sys.argv) > 6 else 1
NTH = int(sys.argv[7]) if len(sys.argv) > 7 else 64
ROOT = sys.argv[8] if len(sys.argv) > 8 else 'r0'          # r0 : projection U(1), lien 0 si moyenne nulle ; r1 : moyenne non projetee
PART = sys.argv[9] if len(sys.argv) > 9 else '0/1'          # repartition des pas de theta : k/K
L = 2 ** N; B = 2 * np.pi / L ** 2; r = 1.0; PHI0 = np.pi / 4
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
GAM = [np.kron(sx, sx), np.kron(sx, sy), np.kron(sx, sz), np.kron(sy, I2), np.kron(sz, I2)]
sides = [2 ** n for n in range(N + 1)]
offs = np.concatenate([[0], np.cumsum([s ** 3 for s in sides])]).astype(int); Nc = int(offs[-1])
depth = np.concatenate([np.full(s ** 3, n) for n, s in enumerate(sides)])
def grid(s): return np.meshgrid(np.arange(s), np.arange(s), np.arange(s), indexing='ij')
def lin(n, I, J, K): s = sides[n]; return offs[n] + ((I % s) * s + J % s) * s + K % s

def fine_and_sampled(n, th, a, k):
    """u(x -> x + e_j) au niveau n par integrale du continu le long du segment entre centres (transport a)"""
    s = sides[n]; d = L / s; I, J, K = grid(s); X, Y = (I + 0.5) * d, (J + 0.5) * d
    u0 = np.where(I == s - 1, np.exp(-1j * B * L * Y), 1.0 + 0j)           # translation magnetique a l'enroulement en x
    u1 = np.exp(1j * (B * X + a * np.cos(k * X - PHI0)) * d)
    u2 = np.exp(1j * th / L * d) * np.ones_like(X)
    return [u0 ** q, u1 ** q, u2 ** q]

def T_in(Uf, n):
    """transport, au niveau n + 1, de l'enfant representant 2x vers l'enfant 2x + c (chemin x puis y puis z) ; tableau [c][...] sur x du niveau n"""
    s = sides[n]; I, J, K = grid(s); out = {}
    for c in itertools.product((0, 1), repeat=3):
        P = np.ones((s, s, s), complex); pos = [2 * I, 2 * J, 2 * K]
        for jd in range(3):
            if c[jd]:
                P = P * Uf[jd][pos[0] % (2 * s), pos[1] % (2 * s), pos[2] % (2 * s)]; pos[jd] = pos[jd] + 1
        out[c] = P
    return out

def link_fields(th, a):
    k = 2 * np.pi * 2 ** JW / L
    U = [None] * (N + 1); TIN = [None] * (N + 1)
    U[N] = fine_and_sampled(N, th, a, k)
    for n in range(N - 1, -1, -1):
        if TR == 'a':
            U[n] = fine_and_sampled(n, th, a, k)
        Uf = U[n + 1]; TI = T_in(Uf, n); TIN[n] = TI
        if TR == 'b':
            s = sides[n]; I, J, K = grid(s); S2 = 2 * s; new = []
            for jd in range(3):
                acc = np.zeros((s, s, s), complex)
                for c, P in TI.items():
                    pos = [(2 * I + c[0]), (2 * J + c[1]), (2 * K + c[2])]
                    a1 = Uf[jd][pos[0] % S2, pos[1] % S2, pos[2] % S2]
                    pos[jd] = pos[jd] + 1
                    a2 = Uf[jd][pos[0] % S2, pos[1] % S2, pos[2] % S2]
                    sh = [0, 0, 0]; sh[jd] = -1
                    Pn = np.roll(P, sh, axis=(0, 1, 2))                       # T_in(x + e_j, c)
                    acc += P * a1 * a2 * np.conj(Pn)
                if ROOT == 'r1': new.append(acc / 8); continue
                mag = np.abs(acc); new.append(np.where(mag / 8 > 1e-6, acc / np.where(mag > 0, mag, 1), 0.0))   # moyenne nulle (racine sous flux) : lien 0
            U[n] = new
    return U, TIN

def vertical_sampled(th, a):
    """u(p -> c) par le segment droit (transport a)"""
    k = 2 * np.pi * 2 ** JW / L; out = []
    for n in range(1, N + 1):
        s = sides[n]; d = L / s; I, J, K = grid(s)
        Xc, Yc, Zc = (I + .5) * d, (J + .5) * d, (K + .5) * d
        Xp, Yp, Zp = (I // 2 + .5) * 2 * d, (J // 2 + .5) * 2 * d, (K // 2 + .5) * 2 * d
        dX, dY = Xc - Xp, Yc - Yp
        pert = a * dY * (np.sin(k * Xc - PHI0) - np.sin(k * Xp - PHI0)) / (k * dX)
        out.append(np.exp(1j * q * (B * 0.5 * (Xp + Xc) * dY + pert + th / L * (Zc - Zp))))
    return out

def hamiltonian(th, a):
    U, TIN = link_fields(th, a)
    R, C, Vv = [[], [], []], [[], [], []], [[], [], []]
    for n in range(N + 1):
        s = sides[n]; I, J, K = grid(s); i = lin(n, I, J, K).ravel()
        for jd, sh in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
            nb = lin(n, I + sh[0], J + sh[1], K + sh[2]).ravel()
            R[jd].append(i); C[jd].append(nb); Vv[jd].append(np.conj(U[n][jd]).ravel())   # T[i, nb] = u(nb -> i)
    T = [sp.csr_matrix((np.concatenate(Vv[jd]), (np.concatenate(R[jd]), np.concatenate(C[jd]))), shape=(Nc, Nc)) for jd in range(3)]
    rs, cs, vs = [], [], []
    VS = vertical_sampled(th, a) if TR == 'a' else None
    for n in range(1, N + 1):
        s = sides[n]; I, J, K = grid(s); ch = lin(n, I, J, K).ravel(); pa = lin(n - 1, I // 2, J // 2, K // 2).ravel()
        if TR == 'a': u = VS[n - 1]
        else:
            u = np.zeros((s, s, s), complex)
            for c, P in TIN[n - 1].items():
                u[c[0]::2, c[1]::2, c[2]::2] = P
        rs.append(ch); cs.append(pa); vs.append(u.ravel() / np.sqrt(8))
    S = sp.csr_matrix((np.concatenate(vs), (np.concatenate(rs), np.concatenate(cs))), shape=(Nc, Nc))
    kr = lambda A, g: sp.kron(A, sp.csr_matrix(g), format='csr')
    H = sp.csr_matrix((4 * Nc, 4 * Nc), dtype=complex)
    for jd in range(3):
        H = H + kr((T[jd] - T[jd].getH()) / 2j, GAM[jd]) - (r / 2) * kr(T[jd] + T[jd].getH(), GAM[4])
    H = H + kr((S - S.getH()) / 2j, GAM[3]) - (r / 2) * kr(S + S.getH(), GAM[4]) + (m + 4 * r) * kr(sp.identity(Nc), GAM[4])
    return H.tocsr(), U

def level_flux(U):
    out = []
    for n in range(N + 1):
        s = sides[n]
        if s < 2: out.append(np.nan); continue
        u0, u1 = U[n][0], U[n][1]
        pl = u0 * np.roll(u1, -1, axis=0) * np.conj(np.roll(u0, -1, axis=1)) * np.conj(u1)
        out.append(np.angle(pl).sum() / (2 * np.pi) / s)                    # par tranche z
    return out

print(f'=== V162 : N = {N} (L = {L}), m = {m}, transport ({TR}{"/" + ROOT if TR == "b" else ""}), q = {q}, mesure {MES}, j = {JW} (l = L/{2**JW}), pas {NTH}, part {PART} ===', flush=True)
H0, U0 = hamiltonian(0.3, 0.0)
print(f'  hermiticite {abs(H0 - H0.getH()).max():.1e} ; flux par niveau (tranche z) : {["%.4f" % f for f in level_flux(U0)]}')
# B-P0 : covariance sous jauge aleatoire du niveau fin (transport b) : u(x -> x+e) -> g(x+e) u g(x)^*, puis reconstruction des niveaux grossiers
if TR == 'b' and MES == 'flot' and 4 * Nc < 6000:
    rng = np.random.default_rng(1); s = sides[N]; g = np.exp(2j * np.pi * rng.random((s, s, s)))
    _orig = fine_and_sampled
    def _gauged(n, th, a, k):
        Uf = _orig(n, th, a, k)
        if n != N: return Uf
        return [np.roll(g, -1, axis=jd) * Uf[jd] * np.conj(g) for jd in range(3)]
    E1 = np.sort(np.linalg.eigvalsh(hamiltonian(0.3, 0.0)[0].toarray()))
    fine_and_sampled = _gauged
    E2 = np.sort(np.linalg.eigvalsh(hamiltonian(0.3, 0.0)[0].toarray()))
    fine_and_sampled = _orig
    print(f'  B-P0 covariance (jauge aleatoire du niveau fin) : ecart spectral {np.abs(E1 - E2).max():.1e}', flush=True)

def eig_window(H, eps, K=24):
    if H.shape[0] <= 6000:
        E, V = np.linalg.eigh(H.toarray())
    else:
        E, V = sla.eigsh(H, k=K, sigma=1e-4, which='LM')
    sel = np.abs(E) < eps; return E[sel], V[:, sel]
EPS = 0.25 * 2 * np.pi / L
def prof(v): return np.bincount(depth, weights=np.sum(np.abs(v.reshape(Nc, 4)) ** 2, axis=1), minlength=N + 1)

if MES == 'flot':
    ths = np.linspace(0, 2 * np.pi, NTH, endpoint=False); dth = ths[1] - ths[0]; flow = np.zeros(N + 1); gross = np.zeros(N + 1); h = 1e-6
    kp, Kp = map(int, PART.split('/')); ths = ths[kp::Kp]
    for th in ths:
        H, _ = hamiltonian(th, 0.0); D = (hamiltonian(th + h, 0.0)[0] - hamiltonian(th - h, 0.0)[0]) / (2 * h)
        E, V = eig_window(H, EPS)
        for i in range(len(E)):
            v = V[:, i]; sl = np.real(np.vdot(v, D @ v)); w = prof(v)
            flow += dth / (2 * EPS) * sl * w; gross += dth / (2 * EPS) * abs(sl) * w
    print(f'  fenetre {EPS:.3f} ; flot net par couche : [{" ".join(f"{x:+.3f}" for x in flow)}] ; total {flow.sum():+.3f}')
    print(f'  flot brut par couche : [{" ".join(f"{x:.3f}" for x in gross)}] ; couche fine : net {flow[N]:+.3f} ; couches <= N-1 : net {flow[:N].sum():+.3f}')
    print('  BRUT_NET', ' '.join(f'{x:.6f}' for x in flow))
elif MES == 'ff':
    h = 1e-6; V = (hamiltonian(0.3, h)[0] - hamiltonian(0.3, -h)[0]) / (2 * h)
    rows = np.repeat(np.arange(Nc), 4)
    ff = []
    for n in range(N + 1):
        P = sp.diags((depth[np.arange(4 * Nc) // 4] == n).astype(float))
        M = P @ V @ P; ff.append(np.sqrt(abs(M.multiply(M.conj())).sum() / (4 * sides[n] ** 3)))
    ff = np.array(ff); print(f'  facteur de forme RMS par cellule f_n / f_N (n = 0 .. {N}) : [{" ".join(f"{x / ff[N]:.3f}" for x in ff)}]')
elif MES == 'p4b':
    h = 1e-6; edge, mir, ALL = [], [], []
    for th in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2):
        H, _ = hamiltonian(th, 0.0); V = (hamiltonian(th, h)[0] - hamiltonian(th, -h)[0]) / (2 * h)
        if H.shape[0] <= 6000:
            E, W = np.linalg.eigh(H.toarray()); o = np.argsort(np.abs(E))[:24]; E, W = E[o], W[:, o]
        else:
            E, W = sla.eigsh(H, k=24, sigma=1e-4, which='LM')
        for i in range(len(E)):
            v = W[:, i]; nb = np.dot(np.arange(N + 1), prof(v)); c = np.linalg.norm(V @ v)
            (edge if nb >= N - 0.5 else mir if nb <= N - 1.5 else []).append((c, E[i], nb)); ALL.append((c, E[i], nb))
    ce = [x[0] for x in edge]; cm = [x[0] for x in mir]
    e, mm = (np.median(ce) if ce else np.nan), (np.median(cm) if cm else np.nan)
    print(f'  24 plus bas etats x 4 theta : bord {len(ce)} etats (|E| max {max([abs(x[1]) for x in edge], default=np.nan):.3f}), mediane {e:.4f} ;'
          f' miroir {len(cm)} etats (|E| max {max([abs(x[1]) for x in mir], default=np.nan):.3f}, n moyen {np.mean([x[2] for x in mir]) if mir else np.nan:.2f}), mediane {mm:.4f} ; rapport {mm / e:.4f}')
    allst = edge + mir + []
    print('  (descriptif, hors critere) les 6 etats de plus petite profondeur moyenne parmi les 96 : ' + ' '.join(f'n{x[2]:.2f}/E{x[1]:+.3f}/c{x[0]:.3f}' for x in sorted(ALL, key=lambda t: t[2])[:6]))
elif MES == 'p4':
    h = 1e-6; edge, mir = [], []
    for th in np.linspace(0, 2 * np.pi, min(NTH, 16), endpoint=False):
        H, _ = hamiltonian(th, 0.0); V = (hamiltonian(th, h)[0] - hamiltonian(th, -h)[0]) / (2 * h)
        E, W = eig_window(H, EPS)
        for i in range(len(E)):
            v = W[:, i]; nb = np.dot(np.arange(N + 1), prof(v)); c = np.linalg.norm(V @ v)
            (edge if nb >= N - 0.5 else mir if nb <= N - 1.5 else []).append(c)
    e, mm = np.median(edge) if edge else np.nan, np.median(mir) if mir else np.nan
    print(f'  couplage ||V psi|| : bord (n >= N-0.5) {len(edge)} etats, mediane {e:.4f} ; compensation (n <= N-1.5) {len(mir)} etats, mediane {mm:.4f} ; rapport {mm / e:.4f}')
print(f'  temps {time.time() - t0:.0f} s', flush=True)
