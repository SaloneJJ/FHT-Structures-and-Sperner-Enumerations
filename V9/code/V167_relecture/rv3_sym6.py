"""Point 3 : operateur explicite de la symetrie theta -> -theta. U = Lambda (R'_z x Gamma_3), R'_z : reflexion z -> L/2 - z
(K -> s/2 - 1 - K mod s au niveau de cote s ; racine fixe), qui laisse fixes les cellules du niveau 1 (enfant representant inchange).
Test : U H(th) U^+ = -H(-th), Lambda jauge diagonale par cellule. Usage : python3 rv3_sym6.py N q"""
import sys, numpy as np, scipy.sparse as sp
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
N, q = int(sys.argv[1]), float(sys.argv[2]); g = C.setup(N, q); Nc = g['Nc']; GAM = g['GAM']; sides, lin, grid = g['sides'], g['lin'], g['grid']
src = open('/home/claude/fht/relecture3/rv3_sym5.py').read(); exec(src.split('def solve')[1].join(['def solve', '']).split('R12 =')[0])
p = np.zeros(Nc, int)
for n in range(N + 1):
    s = sides[n]; I, J, K = grid(s); p[lin(n, I, J, K).ravel()] = lin(n, I, J, (s // 2 - 1 - K) % s if s > 1 else K).ravel()
P = sp.csr_matrix((np.ones(Nc), (p, np.arange(Nc))), shape=(Nc, Nc)); U = sp.kron(P, sp.csr_matrix(GAM[2])).tocsr()
for th in (0.0, 0.37, 2.2, -1.1):
    X = (U @ C.Href(th) @ U.getH()).tocsr(); T = (-C.Href(-th)).tocsr()
    print(f'N = {N}, q = {q}, th = {th} : |Lambda U H(th) U^+ Lambda^+ + H(-th)| = {solve(X, T):.1e}', flush=True)
