"""Profil en profondeur du partenaire (premier etat chi < -0.3, hors racine) de n = (1,0,0), via la reduction exacte (arbre binaire, secteur impair, diagonalisation dense)."""
import numpy as np
from rel_B2_reduction import tree_dD, sector, CHI
for m in (-1.0, -0.5):
    for N in (5, 7, 9):
        cells, H = tree_dD(N, m, 1); P = sector(cells, 1, (1,)); Hs = (P.getH() @ H @ P).toarray()
        E, V = np.linalg.eigh(Hs); Vf = P @ V; depth = np.repeat([n for n, c in cells], 4)
        for a in np.argsort(np.abs(E)):
            p = Vf[:, a].reshape(-1, 4); chi = np.real(np.einsum('ia,ab,ib->', p.conj(), CHI, p))
            if chi < -0.3:
                w = np.bincount(depth, weights=np.abs(Vf[:, a]) ** 2, minlength=N + 1)
                print(f'm = {m}, N = {N} : partenaire E = {abs(E[a]):.5f}, chi = {chi:+.2f}, w(n) = {np.round(w, 2)}'); break
        nlow = int(np.sum(np.abs(E) < 0.7)); nneg = 0
        for a in range(len(E)):
            if abs(E[a]) < 0.7:
                p = Vf[:, a].reshape(-1, 4)
                if np.real(np.einsum('ia,ab,ib->', p.conj(), CHI, p)) < -0.3: nneg += 1
        print(f'        etats |E| < 0.7 dans ce sous-espace : {nlow}, dont chi < -0.3 : {nneg}')
