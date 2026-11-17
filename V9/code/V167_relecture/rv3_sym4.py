"""Point 3 : verification directe, N = 3 et 4 : (i) U = Lambda (R_z x Gamma_3) : U H(th) U^+ = -H(-th) ;
(ii) periodicite : Lambda' H(th) Lambda'^+ = H(th + 2 pi / q) (jauge diagonale). Usage : python3 rv3_sym4.py N q"""
import sys, numpy as np
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
src = open('/home/claude/fht/relecture3/rv3_sym3.py').read().split('gp = {}')[0]
N, q = int(sys.argv[1]), float(sys.argv[2]); sys.argv = ['x', str(N), str(q)]; exec(src)
import scipy.sparse as sp
def test_shift(th, dth):
    X = C.Href(th + dth).tocsr(); T = C.Href(th).tocsr(); Tc = T.tocoo()
    lam = np.full(Nc, np.nan + 0j); lam[0] = 1; stack = [0]; adj = {}
    for r, c, v in zip(Tc.row, Tc.col, Tc.data):
        if abs(v) > 1e-9 and r // 4 != c // 4: adj.setdefault(r // 4, []).append((c // 4, r, c, v))
    while stack:
        a = stack.pop()
        for b, r, c, v in adj.get(a, []):
            if np.isnan(lam[b]): lam[b] = np.conj(X[r, c] / (lam[a] * v)); stack.append(b)
    Lm = sp.kron(sp.diags(lam), sp.identity(4)); return abs(Lm @ T @ Lm.getH() - X).max()
G3 = GAM[2]
for th in (0.0, 0.37, 2.2):
    print(f'N = {N}, q = {q}, th = {th} : (i) R_z x Gamma3, -H(-th) : ecart {test(th, (0, 0, 1), G3, 0, -1):.1e} ; (ii) H(th + 2pi/q) ~ H(th) : ecart {test_shift(th, 2 * np.pi / q):.1e}', flush=True)
