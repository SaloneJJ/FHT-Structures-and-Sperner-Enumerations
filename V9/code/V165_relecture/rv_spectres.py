"""Points 3 et 4 (hors flot) :
(A) unicite de l'operateur singulet sur le tore plat 4D (L = 4 et L = 6), code independant ;
(B) modele KS + masse eps (V159) : plus petit |E| porte par les couches n >= N-1, avec et sans flux, arbre N = 3, 4 et plat N = 3 ;
(C) modele KS + masse singulet (V161) : spectre bas et profondeur moyenne, avec et sans flux (arbre N = 3, plat N = 3)."""
import sys, os, itertools, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_core import *
t0 = time.time()
PART = sys.argv[1] if len(sys.argv) > 1 else 'ABC'

if 'A' in PART:
    for L in (4, 6):
        sites = list(itertools.product(range(L), repeat=4)); ix = {s_: i for i, s_ in enumerate(sites)}; X = np.array(sites); n = len(sites)
        def shift(mu):
            R = np.arange(n); C = np.array([ix[tuple((np.array(s_) + np.eye(4, dtype=int)[mu]) % L)] for s_ in sites])
            return sp.csr_matrix((np.ones(n), (R, C)), shape=(n, n))
        T = [shift(mu) for mu in range(4)]
        eta = [(-1.0) ** X[:, :mu].sum(1) for mu in range(4)]
        K = sum((sp.diags(eta[mu]) @ T[mu] - (sp.diags(eta[mu]) @ T[mu]).getH()) / 2j for mu in range(4))
        zeta = [(-1.0) ** X[:, mu + 1:].sum(1) for mu in range(4)]
        Ssym = [sp.diags(zeta[mu]) @ T[mu] for mu in range(4)]
        C = [(T[mu] + T[mu].T) / 2 for mu in range(4)]; P = C[0] @ C[1] @ C[2] @ C[3]
        nrm = lambda A: abs(A).max() if A.nnz else 0.0
        print(f'(A) L = {L} : [S_mu, K] max = {max(nrm(Ssym[m] @ K - K @ Ssym[m]) for m in range(4)):.1e}')
        ok = []
        for a in itertools.product((0, 1), repeat=4):
            f = (-1.0) ** (X @ np.array(a))
            for ph in (1, 1j):
                G = ph * sp.diags(f) @ P
                h = nrm(G - G.getH()); ac = nrm(G @ K + K @ G); cs = max(nrm(G @ Ssym[m] - Ssym[m] @ G) for m in range(4))
                if h < 1e-12 and ac < 1e-12 and cs < 1e-12: ok.append((a, ph))
        print(f'      vecteurs a (phase) satisfaisant (i) hermitien (ii) anticommute (iii) commute aux 4 decalages : {ok}')
        # spectre de K + M Gamma5 sur le tore (gap de volume)
        a = ok[0][0]; f = (-1.0) ** (X @ np.array(a)); G5 = sp.diags(f) @ P
        if n <= 1296:
            E = np.linalg.eigvalsh((K + G5).toarray()); print(f'      tore L = {L}, M = 1 : min |E| = {np.min(np.abs(E)):.4f} ; {{G5, K}} = {nrm(G5 @ K + K @ G5):.1e} ; G5^2 propre : {np.unique(np.round(np.linalg.eigvalsh((G5 @ G5).toarray()), 6))[:5]}')

def low(H, k=24):
    if H.shape[0] <= 2500:
        E, V = np.linalg.eigh(H.toarray())
    else:
        E, V = sla.eigsh(H.tocsc(), k=k, sigma=1e-6, which='LM')
    o = np.argsort(np.abs(E)); return E[o], V[:, o]

def report(g, H, label, k=24, show=8):
    E, V = low(H, k); w = layer_weights(g, V, 1); nb = np.arange(g.N + 1) @ w
    edge = (w[g.N] + w[g.N - 1]) > 0.5
    me = np.min(np.abs(E[edge])) if edge.any() else np.nan
    print(f'  {label} : min|E| global {np.min(np.abs(E)):.4f} ; min|E| (poids n >= N-1 > 1/2) {me:.4f} ; bas du spectre : ' +
          ' '.join(f'{E[i]:+.3f}(n{nb[i]:.1f})' for i in range(min(show, len(E)))))
    return E, nb

if 'B' in PART:
    print('(B) KS + masse eps (V159), theta = 0')
    for flat, N in ((False, 3), (False, 4), (True, 3)):
        for Phi in (1.0, 0.0):
            g = Geo(N, flat, 'G2', Phi=Phi); _, _, eps = ks_phases(g)
            for M, wall in ((1.0, False), (1.0, True), (0.3, False)):
                if flat and M == 0.3: continue
                H = H_ks_kin(g, 0.0) + sp.diags(mass_profile(g, M, wall) * eps)
                report(g, H.tocsr(), f'{"plat " if flat else "arbre"} N={N} Phi={Phi:.0f} M={M} {"W" if wall else "O"}')
    print(f'  temps {time.time() - t0:.0f} s', flush=True)

if 'C' in PART:
    print('(C) KS + masse singulet (V161), theta = 0')
    for flat, N in ((False, 3), (True, 3), (False, 4)):
        for Phi in (1.0, 0.0):
            g = Geo(N, flat, 'G2', Phi=Phi)
            for wall in (False, True):
                G5 = gamma5_singlet(g, 0.0)
                H = (H_ks_kin(g, 0.0) + sp.diags(mass_profile(g, 1.0, wall)) @ G5).tocsr()
                E, nb = report(g, H, f'{"plat " if flat else "arbre"} N={N} Phi={Phi:.0f} {"W" if wall else "O"}', k=30, show=10)
                print(f'      nombre d etats |E| < 0.2 : {np.sum(np.abs(E) < 0.2)} (dont n moyen <= N-1.5 : {np.sum((np.abs(E) < 0.2) & (nb <= N - 1.5))})'
                      f' ; |G5| max element {abs(G5).max():.4f} ; norme d operateur ~ {np.max(np.abs(sla.eigsh(G5, k=2, which="LM", return_eigenvectors=False))):.3f}')
        print(f'  temps {time.time() - t0:.0f} s', flush=True)
