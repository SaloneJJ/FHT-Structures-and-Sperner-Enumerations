"""V339 relecture, partie 5b : à N = 4, valeurs propres de |H| les plus petites par ARPACK (Lanczos) sur H² sans inversion (méthode indépendante de LOBPCG),
comparées à LOBPCG (procédure de V338, départ aléatoire puis départ chaud) sur les mêmes θ. Usage : python3 rl5b_arpack_N4.py ITH1 ITH2 ..."""
import sys, os, time, numpy as np, scipy.sparse.linalg as sla
here = os.path.dirname(os.path.abspath(__file__))
idx = [int(x) for x in sys.argv[1:]] or [0, 1]
ns = {'__name__': 'rl2', '__file__': os.path.join(here, 'rl2_flot.py')}
sys.argv = ['rl2', '4', '64', 'alt', '1', '0', '1', '1']
exec(open(os.path.join(here, 'rl2_flot.py')).read(), ns)
Hs, EPS, prof, N = ns['Hs'], ns['EPS'], ns['prof'], ns['N']
ths = np.linspace(0, 2 * np.pi, 64, endpoint=False); h = 1e-6; Xprev = None; KL = 64
def contrib(E, V, D):
    out = []
    for i in np.argsort(np.abs(E)):
        if abs(E[i]) < EPS:
            v = V[:, i]; out.append((E[i], np.real(np.vdot(v, D @ v)), prof(v)[N]))
    return out
for j in idx:
    th = ths[j]; H = Hs(th); D = (Hs(th + h) - Hs(th - h)) / (2 * h)
    H2 = sla.LinearOperator(H.shape, matvec=lambda x: H @ (H @ x), matmat=lambda X_: H @ (H @ X_), dtype=complex)
    t = time.time()
    X0 = Xprev if Xprev is not None else np.random.default_rng(0).normal(size=(H.shape[0], KL)) + 0j
    w_, X_ = sla.lobpcg(H2, X0, largest=False, tol=1e-7, maxiter=400)
    Q_, _ = np.linalg.qr(X_); Hsub = Q_.conj().T @ (H @ Q_); e_, c_ = np.linalg.eigh((Hsub + Hsub.conj().T) / 2)
    El, Vl = e_, Q_ @ c_; Xprev = X_; tl = time.time() - t
    t = time.time()
    w2, Z = sla.eigsh(H2, k=40, which='SA', tol=1e-10, ncv=160, maxiter=200000)
    Hsub = Z.conj().T @ (H @ Z); e2, c2 = np.linalg.eigh((Hsub + Hsub.conj().T) / 2); Ea, Va = e2, Z @ c2; ta = time.time() - t
    res = np.linalg.norm(H @ Va - Va * Ea, axis=0)
    cl, ca = contrib(El, Vl, D), contrib(Ea, Va, D)
    print(f'θ = {th:.4f} (pas {j}) : LOBPCG {len(cl)} états |E| < ε ({tl:.0f} s) ; ARPACK-SA {len(ca)} ({ta:.0f} s, résidu max {res.max():.1e}) ; '
          f'plus petits |E| : LOBPCG {np.round(np.sort(np.abs(El))[:6], 5)} ; ARPACK {np.round(np.sort(np.abs(Ea))[:6], 5)}', flush=True)
    for (e1, s1, w1), (e2_, s2, w2) in zip(cl, ca):
        print(f'    E {e1:+.6f} / {e2_:+.6f}   pente {s1:+.5f} / {s2:+.5f}   poids couche fine {w1:.4f} / {w2:.4f}')
    print(f'    Σ pente × poids(couche fine) : LOBPCG {sum(s*w for _, s, w in cl):+.5f} ; ARPACK {sum(s*w for _, s, w in ca):+.5f}', flush=True)
