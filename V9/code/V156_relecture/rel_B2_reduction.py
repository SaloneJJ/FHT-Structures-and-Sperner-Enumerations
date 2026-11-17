"""Reduction EXACTE du modele 3D des auteurs au sous-espace des fonctions independantes de (y, z) [resp. de z].
Ce sous-espace est invariant par H (sauts lateraux et S, S^+ preservent les constantes en y, z) ; sur lui, sin k_y = sin k_z = 0,
le Wilson en y, z s'annule, et le poids parent-enfant devient 1/sqrt(2) (resp. 1/2) : arbre binaire (resp. quaternaire),
Wilson-Dirac (1 ou 2 directions laterales + verticale). La couche |n|^2 = 1 (n = (1,0,0)) vit dans le secteur impair en x
de l'arbre binaire ; |n|^2 = 2 (n = (1,1,0)) dans le secteur (impair, impair) de l'arbre quaternaire.
Verification : comparaison directe avec l'arbre 3D (rel_tree.build) a N = 3."""
import sys, itertools, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
s0 = np.eye(2); s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1., -1.]).astype(complex)
G = [np.kron(s1, s3), np.kron(s2, s3), np.kron(s3, s3), np.kron(s0, s1), np.kron(s0, s2)]; CHI = -1j * G[3] @ G[4]
r = 1.0
def tree_dD(N, m, D):
    """arbre de dimension laterale D (D = 1, 2, 3), fermions a 4 composantes, Wilson total m + 4r (les directions absentes ont k = 0)."""
    cells = []; index = {}
    for n in range(N + 1):
        for c in itertools.product(range(2 ** n), repeat=D): index[(n, c)] = len(cells); cells.append((n, c))
    Nc = len(cells); R, C, V = [], [], []
    def put(i, j, M):
        for a in range(4):
            for b in range(4):
                if M[a, b] != 0: R.append(4 * i + a); C.append(4 * j + b); V.append(M[a, b])
    w = 1 / np.sqrt(2 ** D)
    for i, (n, c) in enumerate(cells):
        s = 2 ** n
        put(i, i, (m + 4 * r - (3 - D) * r) * G[4])
        for d in range(D):
            c2 = list(c); c2[d] = (c2[d] + 1) % s; j = index[(n, tuple(c2))]
            M = G[d] / 2j - r / 2 * G[4]; put(i, j, M); put(j, i, M.conj().T)
        if n >= 1:
            p = index[(n - 1, tuple(x // 2 for x in c))]; M = w * (G[3] / 2j - r / 2 * G[4]); put(i, p, M); put(p, i, M.conj().T)
    H = sp.csr_matrix((V, (R, C)), shape=(4 * Nc, 4 * Nc))
    return cells, H
def sector(cells, D, eps):
    rows, cols, vals = [], [], []; col = 0; idx = {cc: i for i, cc in enumerate(cells)}
    for (n, c) in cells:
        if n == 0: continue
        h = 2 ** (n - 1)
        if any(x >= h for x in c): continue
        for dl in itertools.product((0, 1), repeat=D):
            rows.append(idx[(n, tuple(x + d_ * h for x, d_ in zip(c, dl)))]); cols.append(col); vals.append((-1) ** sum(e * d_ for e, d_ in zip(eps, dl)) / np.sqrt(2 ** D))
        col += 1
    P = sp.csr_matrix((vals, (rows, cols)), shape=(len(cells), col)); return sp.kron(P, sp.identity(4)).tocsr()
def low(N, m, D, eps, k=16):
    cells, H = tree_dD(N, m, D); P = sector(cells, D, eps); Hs = (P.getH() @ H @ P).tocsc()
    if Hs.shape[0] < 2500: E, V = np.linalg.eigh(Hs.toarray())
    else:
        E, V = sla.eigsh(Hs, k=k, sigma=0.05 * np.sin(2 * np.pi / 2 ** N), which='LM')
        E2_, V2_ = sla.eigsh(Hs, k=12, sigma=0.55, which='LM')       # region des partenaires
        E = np.concatenate([E, E2_]); V = np.concatenate([V, V2_], axis=1)
    Vf = P @ V; depth = np.repeat([n for n, c in cells], 4)
    out = []
    for a in np.argsort(np.abs(E)):
        if E[a] <= 1e-9: continue
        v = Vf[:, a]; p = v.reshape(-1, 4); chi = np.real(np.einsum('ia,ab,ib->', p.conj(), CHI, p)); nb = np.real(np.sum(np.abs(v) ** 2 * depth))
        out.append((E[a], chi, nb))
    return out
if __name__ == '__main__':
    # (1) verification de la reduction contre l'arbre 3D complet (N = 3), secteur 100 de la relecture
    from rel_tree import build
    cells3, idx3, K, W, depth3 = build(3, -0.5); E3 = np.sort(np.linalg.eigvalsh((K + W).toarray()))
    o1 = low(3, -0.5, 1, (1,)); o2 = low(3, -0.5, 2, (1, 1))
    for lab, o in (('n=(1,0,0) via arbre binaire', o1), ('n=(1,1,0) via arbre quaternaire', o2)):
        e = o[0][0]; print(f'N = 3, m = -0.5 : {lab} : E = {e:.6f} ; present dans le spectre 3D complet : {np.min(np.abs(E3 - e)):.1e}')
    Nmax1 = int(sys.argv[1]) if len(sys.argv) > 1 else 12; Nmax2 = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    for m in (-1.0, -0.5, -1.5):
        print(f'\nm = {m}')
        res1 = {}
        for N in range(3, Nmax1 + 1):
            o = low(N, m, 1, (1,)); E1 = [x for x in o if x[1] > 0.3][0][0]; part = [x for x in o if x[1] < -0.3]
            res1[N] = E1
            pt = f'partenaire (chi -) : E = {part[0][0]:.5f}, n_moy = {part[0][2]:.2f}, rapport {part[0][0]/E1:.2f}' if part else 'partenaire non trouve parmi les etats calcules'
            print(f'  N = {N:2d} : E(1,0,0) = {E1:.5f} ; v = {E1/np.sin(2*np.pi/2**N):.4f} ; {pt}', flush=True)
        for N in range(3, Nmax2 + 1):
            o = low(N, m, 2, (1, 1)); E2 = [x for x in o if x[1] > 0.3][0][0]
            print(f'  N = {N:2d} : E(1,1,0) / E(1,0,0) = {E2/res1[N]:.4f}  (continu : 1.4142 ; reseau plat : {np.sqrt(2)*np.sin(2*np.pi/2**N)/np.sin(2*np.pi/2**N):.4f} -> sin(k)sqrt2/sin k ; exact plat = {np.sqrt(2*np.sin(2*np.pi/2**N)**2)/np.sin(2*np.pi/2**N):.4f})', flush=True)
