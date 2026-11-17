"""Point 3 (suite) : pour chaque etat propre (E, v) de H(theta), l etat de H(-theta) d energie -E a-t-il la meme pente
dE/dtheta et le meme profil par couche ? Quelle reflexion spatiale relie les densites ? (diagonalisation complete)
Usage : python3 rv3_sym2.py N q"""
import sys, itertools, numpy as np
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
N, q = int(sys.argv[1]), float(sys.argv[2]); g = C.setup(N, q); Nc = g['Nc']; L = g['L']
ks, Cf, keys, n4, *_ = C.fourier_fit(kmax_factor=1.0); fast = C.Fast(ks, Cf, keys, n4)
sides, lin, grid = g['sides'], g['lin'], g['grid']
def perm(refl):
    p = np.zeros(Nc, int)
    for n in range(N + 1):
        s = sides[n]; I, J, K = grid(s); src = lin(n, I, J, K).ravel()
        I2, J2, K2 = [(s - 1 - A) if f else A for A, f in zip((I, J, K), refl)]
        p[src] = lin(n, I2, J2, K2).ravel()
    return p
for th in (0.37, 2.2):
    Ep, Vp = np.linalg.eigh(fast.H(th).toarray()); Em, Vm = np.linalg.eigh(fast.H(-th).toarray())
    Dp, Dm = fast.H(th, 1), fast.H(-th, 1)
    sp_ = np.real(np.einsum('ij,ij->j', Vp.conj(), Dp @ Vp)); sm_ = np.real(np.einsum('ij,ij->j', Vm.conj(), Dm @ Vm))
    # appariement E_i(th) <-> -E_j(-th) : Em croissant <-> -Ep decroissant
    j = np.arange(len(Ep))[::-1]
    gap = np.min(np.diff(Ep))
    print(f'theta = {th} : |E(-th) + E(th)| max {np.abs(Em[j] + Ep).max():.1e} (plus petit ecart entre niveaux {gap:.1e})')
    print(f'   |pente(-th, -E) - pente(th, E)| max {np.abs(sm_[j] - sp_).max():.1e} ; |pente(-th,-E) + pente(th,E)| max {np.abs(sm_[j] + sp_).max():.1e}')
    Pp, Pm = C.profs(Vp), C.profs(Vm)
    print(f'   |profil par couche(-th, -E) - profil(th, E)| max {np.abs(Pm[:, j] - Pp).max():.1e}')
    rho_p = np.sum(np.abs(Vp.reshape(Nc, 4, -1)) ** 2, axis=1); rho_m = np.sum(np.abs(Vm.reshape(Nc, 4, -1)) ** 2, axis=1)
    sel = np.argsort(np.abs(Ep))[:40]
    for refl in itertools.product((0, 1), repeat=3):
        p = perm(refl); rr = np.zeros_like(rho_p); rr[p] = rho_p
        print(f'   densite par cellule, reflexion {refl} : ecart max (40 etats les plus bas) {np.abs(rr[:, sel] - rho_m[:, j[sel]]).max():.1e}')
