"""Paire k = 0 hybridee bord-racine : ecart Delta(N) compare a E(1,0,0) (unites physiques), via la reduction exacte au sous-espace
independant de (y, z) puis au secteur pair (contient k = 0 et la racine)."""
import numpy as np, scipy.sparse as sp, itertools
from rel_B2_reduction import tree_dD, sector, low
def sector_even_with_root(cells):
    rows, cols, vals = [0], [0], [1.0]; col = 1; idx = {cc: i for i, cc in enumerate(cells)}
    for (n, c) in cells:
        if n == 0: continue
        h = 2 ** (n - 1)
        if c[0] >= h: continue
        rows += [idx[(n, c)], idx[(n, (c[0] + h,))]]; cols += [col, col]; vals += [2 ** -0.5, 2 ** -0.5]; col += 1
    return sp.kron(sp.csr_matrix((vals, (rows, cols)), shape=(len(cells), col)), sp.identity(4)).tocsr()
for m in (-0.7, -0.5, -1.5, -0.35):
    out = []
    for N in range(4, 10):
        cells, H = tree_dD(N, m, 1); P = sector_even_with_root(cells); E = np.sort(np.abs(np.linalg.eigvalsh((P.getH() @ H @ P).toarray())))
        d = E[0]; E1 = [x for x in low(N, m, 1, (1,)) if x[1] > 0.3][0][0]
        out.append(f'N={N}: Delta={d:.2e}, Delta/E1={d/E1:.3f}')
    print(f'm = {m} (|1+m| = {abs(1+m):.2f}, (2|1+m|) = {2*abs(1+m):.2f}) : ' + ' ; '.join(out), flush=True)
