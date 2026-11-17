# Controle E3 (V145 §3) : modele calcule (contacts par coins seulement entre noeuds "spatiaux", + adjacence fraternelle)
# contre le modele decrit dans la note (les coins de CHAQUE copie engendrent un petit cube : (ii)-(iii) appliques partout).
import io, contextlib, itertools, numpy as np
from collections import defaultdict, Counter
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
with contextlib.redirect_stdout(io.StringIO()):
    import V145_FHT_contact_bord as B
C, Fp, tau, voisin, V, nV, idx, adj = B.C, B.Fp, B.tau, B.voisin, B.V, B.nV, B.idx, B.adj
Cs = set(C)
def graphe(n, partout):
    W = list(itertools.product(V, repeat=n)); pos = {w: i for i, w in enumerate(W)}
    uf = B.UF(); glued = []
    for m in range(1, n + 1):
        Wm = list(itertools.product(V, repeat=m - 1))
        new = []
        for w in Wm:                                   # (ii) dans la copie posee en w
            if partout or all(x in Cs for x in w):
                for (c, k), c2 in voisin.items(): new.append((w + (c,), w + (c2,), k))
        for (u, v, k) in glued:                        # (iii)
            for c in Fp[k]: new.append((u + (c,), v + (tau[k][c],), k))
        glued = new
    for (u, v, k) in glued:
        for x in Fp[k]: uf.u((u, x), (v, tau[k][x]))
    rows, cols = [], []
    byc = defaultdict(list)
    for (u, v, k) in glued:
        for w in (u, v):
            for x in C: byc[uf.f((w, x))].append(pos[w])
    for lst in byc.values():
        lst = sorted(set(lst))
        for a in lst:
            for b in lst:
                if a != b: rows.append(a); cols.append(b)
    for w in W:                                        # adjacence fraternelle (hote)
        for b in V:
            if b != w[-1] and adj[idx[w[-1]], idx[b]]: rows.append(pos[w]); cols.append(pos[w[:-1] + (b,)])
    G = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(W),) * 2).tocsr()
    nc, lab = connected_components(G, directed=False)
    return nc, sorted(Counter(Counter(lab).values()).items(), reverse=True)
for n in (2, 3):
    for partout in (False, True):
        nc, tailles = graphe(n, partout)
        print(f"n={n} ({'(ii)-(iii) dans toutes les copies' if partout else 'modele calcule en V145 (noeuds spatiaux seulement)'}) : {nc} composantes ; (taille, nombre) = {tailles}")
