"""Controle de la question B2 (vitesse, exposant z) sur un ANALOGUE 1+1 D, calculable a grand N :
arbre binaire ferme en cercle (profondeur n : 2^n cellules, parent = floor(c/2), poids 1/sqrt(2) -> S^+S = 1),
Wilson-Dirac 2D : Gamma_lat = sigma_x, Gamma_vert = sigma_y, masse sigma_z, terme sur place (m + 2r).
Bord : fermion chiral E = +/- v p. On mesure v_N = E_1 / sin(2 pi / 2^N) et z_eff = log2(E_1(N-1)/E_1(N)),
compares au reseau plat (meme pochoir, couches de 2^N cellules).
Remarque (apres coup) : ce n'est pas un simple analogue. Pour la couche n = (1,0,0), c'est la REDUCTION EXACTE du modele 3D des auteurs
au sous-espace invariant des fonctions independantes de (y, z) (voir rel_B2_reduction.py, verifie contre le spectre 3D complet a N = 3) ;
les spineurs a 4 composantes s'y scindent en deux copies du modele a 2 composantes ci-dessous."""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex)
r = 1.0
def H_tree(N, m, flat=False):
    sides = [2 ** N if flat else 2 ** n for n in range(N + 1)]; offs = np.concatenate([[0], np.cumsum(sides)]); Nc = offs[-1]
    rT, cT, rS, cS, vS = [], [], [], [], []
    for n, s in enumerate(sides):
        c = np.arange(s); rT.append(offs[n] + c); cT.append(offs[n] + (c + 1) % s)
        if n >= 1:
            rS.append(offs[n] + c); cS.append(offs[n - 1] + (c if flat else c // 2)); vS.append(np.full(s, 1.0 if flat else 1 / np.sqrt(2)))
    T = sp.csr_matrix((np.ones(Nc), (np.concatenate(rT), np.concatenate(cT))), shape=(Nc, Nc))
    S = sp.csr_matrix((np.concatenate(vS), (np.concatenate(rS), np.concatenate(cS))), shape=(Nc, Nc))
    kr = lambda A, g: sp.kron(A, sp.csr_matrix(g))
    H = kr((T - T.T) / 2j, sx) - r / 2 * kr(T + T.T, sz) + kr((S - S.T) / 2j, sy) - r / 2 * kr(S + S.T, sz) + (m + 2 * r) * kr(sp.identity(Nc), sz)
    return H.tocsc(), offs
def odd_sector(offs, N, flat):
    # secteur impair sous la demi-translation : contient l'impulsion comobile n = 1, exclut k = 0 (paire bord-racine hybridee)
    rows, cols, vals = [], [], []; col = 0
    for n in range(1, N + 1):
        s = offs[n + 1] - offs[n]; h = s // 2
        for c in range(h):
            rows += [offs[n] + c, offs[n] + c + h]; cols += [col, col]; vals += [1 / np.sqrt(2), -1 / np.sqrt(2)]; col += 1
    P = sp.csr_matrix((vals, (rows, cols)), shape=(offs[-1], col))
    return sp.kron(P, sp.identity(2)).tocsr()
def E1(N, m, flat=False):
    H, offs = H_tree(N, m, flat)
    P = odd_sector(offs, N, flat); H = (P.getH() @ H @ P).tocsc()
    if H.shape[0] <= 3000:
        E = np.linalg.eigvalsh(H.toarray())
    else:
        guess = 1.3 * np.sin(2 * np.pi / 2 ** N)
        E = sla.eigsh(H, k=12, sigma=guess * 0.7, which='LM', return_eigenvectors=False)
    E = np.sort(np.abs(E)); return E[E > 1e-7][0]
Ns = list(range(3, int(sys.argv[1]) + 1)) if len(sys.argv) > 1 else list(range(3, 15))
for m in (-1.0, -0.7, -0.5, -1.5, -0.25, -1.75):
    lam = abs(1 + m / r)
    pred = (1 - lam ** 2) / (1 - 2 * lam ** 2) if 2 * lam ** 2 < 1 else float('inf')
    vals = [(N, E1(N, m)) for N in Ns]
    v = [e / np.sin(2 * np.pi / 2 ** N) for N, e in vals]
    z = [np.log2(vals[i - 1][1] / vals[i][1]) for i in range(1, len(vals))]
    print(f'm = {m:+.2f} (|1+m/r| = {lam:.2f} ; 2 lambda^2 = {2*lam**2:.3f} ; v predit au 1er ordre {pred:.3f})')
    print('   N    : ' + ' '.join(f'{N:6d}' for N in Ns))
    print('   v_N  : ' + ' '.join(f'{x:6.3f}' for x in v))
    print('   z_eff:        ' + ' '.join(f'{x:6.3f}' for x in z), flush=True)
vf = [(N, E1(N, -1.0, flat=True)) for N in Ns[:6]]
print('plat m = -1 : v_N = ' + ' '.join(f'{e/np.sin(2*np.pi/2**N):.4f}' for N, e in vf) + ' ; z_eff = ' + ' '.join(f'{np.log2(vf[i-1][1]/vf[i][1]):.3f}' for i in range(1, len(vf))))
