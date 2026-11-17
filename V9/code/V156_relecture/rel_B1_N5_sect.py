"""N = 5, m = -1 : secteurs 100 et 000 recalcules avec le code de la relecture (arbre 3D complet, projection propre sur les secteurs)."""
import sys, time, numpy as np, scipy.sparse.linalg as sla
from rel_B2_reduction import tree_dD, sector, CHI
t0 = time.time(); N = 5; cells, H = tree_dD(N, -1.0, 3); depth = np.repeat([n for n, c in cells], 4)
win = 1.1 * np.sin(2 * np.pi * 2 / 2 ** N)
for eps in [tuple(int(x) for x in a) for a in sys.argv[1:]] or [(1, 0, 0), (0, 0, 0)]:
    P = sector(cells, 3, eps) if eps != (0, 0, 0) else None
    if P is None:   # secteur 000 : ajouter la racine
        import scipy.sparse as sp
        Pn = sector(cells, 3, eps); r0 = sp.csr_matrix(([1.0] * 4, (list(range(4)), list(range(4)))), shape=(H.shape[0], 4)); P = sp.hstack([r0, Pn]).tocsr()
    Hs = (P.getH() @ H @ P).tocsc()
    E, V = sla.eigsh(Hs, k=70, sigma=1e-4, which='LM'); Vf = P @ V
    o = np.argsort(np.abs(E)); cnt = {}
    for a in o:
        if abs(E[a]) >= win: continue
        p = Vf[:, a].reshape(-1, 4); chi = np.real(np.einsum('ia,ab,ib->', p.conj(), CHI, p)); nb = np.real(np.sum(np.abs(Vf[:, a]) ** 2 * depth))
        key = ('bord' if nb >= N - 1.5 else ('racine' if nb <= 1 else 'interieur'), '+' if chi > 0.3 else ('-' if chi < -0.3 else '0')); cnt[key] = cnt.get(key, 0) + 1
    neg = sorted([abs(E[a]) for a in o if np.real(np.einsum('ia,ab,ib->', Vf[:, a].reshape(-1, 4).conj(), CHI, Vf[:, a].reshape(-1, 4))) < -0.3 and np.real(np.sum(np.abs(Vf[:, a]) ** 2 * depth)) > 1])
    print(f'secteur {eps} : dim {Hs.shape[0]} ; max |E| calcule {np.abs(E).max():.4f} ; sous la fenetre {win:.4f} : {cnt} ; helicite - hors racine la plus basse : {neg[:3]}  ({time.time()-t0:.0f} s)', flush=True)
