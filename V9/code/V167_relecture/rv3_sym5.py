"""Point 3 : recherche elargie de l operateur reliant H(th) et +-H(-th) : operations spatiales du groupe du carre en (x, y)
(dont l echange x <-> y) x reflexion z, conjugaison complexe ou non, matrice 4x4 = produit de Gamma x rotation (1 + G1 G2)/sqrt2 ^k,
jauge diagonale par cellule. Usage : python3 rv3_sym5.py N q"""
import sys, itertools, numpy as np, scipy.sparse as sp
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
N, q = int(sys.argv[1]), float(sys.argv[2]); g = C.setup(N, q); Nc = g['Nc']; GAM = g['GAM']; sides, lin, grid = g['sides'], g['lin'], g['grid']
def perm(op):
    swap, rx, ry, rz = op; p = np.zeros(Nc, int)
    for n in range(N + 1):
        s = sides[n]; I, J, K = grid(s); src = lin(n, I, J, K).ravel()
        A, B_ = (J, I) if swap else (I, J)
        A = s - 1 - A if rx else A; B_ = s - 1 - B_ if ry else B_; K2 = s - 1 - K if rz else K
        p[src] = lin(n, A, B_, K2).ravel()
    return p
def solve(X, T):
    Tc = T.tocoo(); lam = np.full(Nc, np.nan + 0j); lam[0] = 1; stack = [0]; adj = {}
    for r, c, v in zip(Tc.row, Tc.col, Tc.data):
        if abs(v) > 1e-9 and r // 4 != c // 4: adj.setdefault(r // 4, []).append((c // 4, r, c, v))
    while stack:
        a = stack.pop()
        for b, r, c, v in adj.get(a, []):
            if np.isnan(lam[b]):
                x = X[r, c]
                if abs(x) < 1e-12: return np.inf
                lam[b] = np.conj(x / (lam[a] * v)); stack.append(b)
    if np.isnan(lam).any(): return np.inf
    Lm = sp.kron(sp.diags(lam), sp.identity(4)); return abs(Lm @ T @ Lm.getH() - X).max()
R12 = (np.eye(4) + GAM[0] @ GAM[1]) / np.sqrt(2)
mats = []
for sub in itertools.product((0, 1), repeat=5):
    M = np.eye(4, dtype=complex)
    for j, b in enumerate(sub):
        if b: M = M @ GAM[j]
    for k in range(2): mats.append((sub, k, M @ np.linalg.matrix_power(R12, k)))
th = 0.37; Hp = C.Href(th).tocsr(); Hm = C.Href(-th).tocsr(); found = 0
for conj in (0, 1):
    H1 = Hp.conj() if conj else Hp
    for op in itertools.product((0, 1), repeat=4):
        P = sp.csr_matrix((np.ones(Nc), (perm(op), np.arange(Nc))), shape=(Nc, Nc))
        for sub, k, M in mats:
            U = sp.kron(P, sp.csr_matrix(M)).tocsr(); X = (U @ H1 @ U.getH()).tocsr()
            for sign in (1, -1):
                T = sign * Hm
                if abs(abs(X).sum() - abs(T).sum()) > 1e-8: continue
                e = solve(X, T)
                if e < 1e-9:
                    found += 1; print(f'  conj {conj}, (echange, Rx, Ry, Rz) = {op}, Gamma{[i+1 for i, b in enumerate(sub) if b]} x R12^{k}, signe {sign:+d} : ecart {e:.1e}', flush=True)
print(f'{found} operateurs trouves (N = {N}, q = {q})')
