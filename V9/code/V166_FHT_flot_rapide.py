"""V166 : flot spectral rapide (composantes de Fourier exactes de H(theta), derivee analytique). Criteres : V166 (ecrits avant calcul).
Modele et transport : ceux de V162_FHT_R3_jauge.py (Gamma empruntes de V151, transport (a) ou (b) avec regle r0/r1).
Usage : python3 V166_FHT_flot_rapide.py N m transport q regle NTH part(k/K)    JJ SALONE, FHT 2026"""
import sys, time
N_, m_, TR_, q_, ROOT_, NTH_, PART_ = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]), sys.argv[7]
sys.argv = ['x', N_, m_, TR_, q_, 'none', '1', str(NTH_), ROOT_, PART_]
exec(open('V162_FHT_R3_jauge.py').read().split("print(f'=== V162")[0])
t0 = time.time()
n4 = 4 * Nc
kmax = int(round(2 * L * abs(q))); M = 2 * kmax + 1
thm = 2 * np.pi * 2 * L * np.arange(M) / M
# memoire (V166) : on accumule la transformee de Fourier au fil des constructions, sans garder les H(theta_m)
h0 = hamiltonian(0.1234, 0.0)[0].tocoo(); keys = np.unique(h0.row.astype(np.int64) * n4 + h0.col)
ks = np.arange(-kmax, kmax + 1); F = np.zeros((len(ks), len(keys)), complex)
for i, t in enumerate(thm):
    h = hamiltonian(t, 0.0)[0].tocoo(); kk = h.row.astype(np.int64) * n4 + h.col; pos = np.searchsorted(keys, kk)
    assert np.all(keys[np.minimum(pos, len(keys) - 1)] == kk), 'motif creux variable'
    d = np.zeros(len(keys), complex); np.add.at(d, pos, h.data)
    F += np.outer(np.exp(-2j * np.pi * ks * i / M), d) / M
    del h, d
keep = np.abs(F).max(axis=1) > 1e-13; ks, F = ks[keep], F[keep]
rows, cols = (keys // n4).astype(int), (keys % n4).astype(int)
def Hfast(th, deriv=False):
    ph = np.exp(1j * th * ks / (2 * L)); w = (1j * ks / (2 * L)) * ph if deriv else ph
    return sp.csr_matrix((w @ F, (rows, cols)), shape=(n4, n4))
print(f'=== V166 : N = {N} (L = {L}), m = {m}, transport {TR}{"/" + ROOT if TR == "b" else ""}, q = {q}, {NTH} pas, part {PART} ===', flush=True)
print(f'  composantes de Fourier : {len(ks)} non nulles sur {M} (construction {time.time() - t0:.0f} s)')
rng = np.random.default_rng(3); err = 0; errd = 0
for t in rng.random(3) * 2 * np.pi:
    err = max(err, abs(Hfast(t) - hamiltonian(t, 0.0)[0]).max())
    Dn = (hamiltonian(t + 1e-6, 0.0)[0] - hamiltonian(t - 1e-6, 0.0)[0]) / 2e-6; errd = max(errd, abs(Hfast(t, True) - Dn).max())
print(f'  controle : |H reconstruit - H(V162)| max {err:.1e} ; derivee analytique vs differences finies {errd:.1e}', flush=True)
EPS = 0.25 * 2 * np.pi / L
def prof(v): return np.bincount(depth, weights=np.sum(np.abs(v.reshape(Nc, 4)) ** 2, axis=1), minlength=N + 1)
ths = np.linspace(0, 2 * np.pi, NTH, endpoint=False); dth = ths[1] - ths[0]
kp, Kp = map(int, PART.split('/')); ths = ths[kp::Kp]
flow = np.zeros(N + 1); gross = np.zeros(N + 1); sat = 0
for th in ths:
    H = Hfast(th); Dm = Hfast(th, True)
    E, V = sla.eigsh(H, k=16, sigma=1e-4, which='LM')                    # creux pour toutes les tailles (V166)
    if np.sum(np.abs(E) < EPS) >= 15: sat += 1
    sl_ = np.abs(E) < EPS
    for i in np.where(sl_)[0]:
        v = V[:, i]; s_ = np.real(np.vdot(v, Dm @ v)); w = prof(v)
        flow += dth / (2 * EPS) * s_ * w; gross += dth / (2 * EPS) * abs(s_) * w
print(f'  fenetre {EPS:.3f} ; flot net par couche : [{" ".join(f"{x:+.3f}" for x in flow)}] ; total {flow.sum():+.3f} ; fenetres saturees {sat}')
print(f'  flot brut par couche : [{" ".join(f"{x:.3f}" for x in gross)}]')
print('  BRUT_NET', ' '.join(f'{x:.6f}' for x in flow)); print('  BRUT_GROSS', ' '.join(f'{x:.6f}' for x in gross))
print(f'  temps {time.time() - t0:.0f} s', flush=True)
