"""Point 3 : operateur explicite. Teste U = Lambda (R_z (x) Gamma_3), Lambda jauge diagonale par cellule a determiner :
U H(theta) U^+ = -H(-theta) ? Usage : python3 rv3_sym3.py N q"""
import sys, itertools, numpy as np, scipy.sparse as sp
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
N, q = int(sys.argv[1]), float(sys.argv[2]); g = C.setup(N, q); Nc = g['Nc']; GAM = g['GAM']
sides, lin, grid = g['sides'], g['lin'], g['grid']
def perm(refl):
    p = np.zeros(Nc, int)
    for n in range(N + 1):
        s = sides[n]; I, J, K = grid(s); src = lin(n, I, J, K).ravel()
        I2, J2, K2 = [(s - 1 - A) if f else A for A, f in zip((I, J, K), refl)]
        p[src] = lin(n, I2, J2, K2).ravel()
    return p
def test(th, refl, gam, conj, sign):
    p = perm(refl); P = sp.csr_matrix((np.ones(Nc), (p, np.arange(Nc))), shape=(Nc, Nc))
    U = sp.kron(P, sp.csr_matrix(gam)).tocsr()
    H = C.Href(th); H = H.conj() if conj else H
    X = (U @ H @ U.getH()).tocsr(); T = (sign * C.Href(-th)).tocsr()
    # jauge par cellule : X[4a+i, 4b+j] = lam_a conj(lam_b) T[4a+i, 4b+j]
    lam = np.full(Nc, np.nan + 0j); lam[0] = 1; stack = [0]; Tc = T.tocoo()
    adj = {}
    for r, c, v in zip(Tc.row, Tc.col, Tc.data):
        if abs(v) > 1e-9 and r // 4 != c // 4: adj.setdefault(r // 4, []).append((c // 4, r, c, v))
    while stack:
        a = stack.pop()
        for b, r, c, v in adj.get(a, []):
            if np.isnan(lam[b]):
                x = X[r, c]
                if abs(x) < 1e-9: return np.inf
                lam[b] = np.conj(x / (lam[a] * v)); stack.append(b)
    if np.isnan(lam).any(): return np.inf
    Lm = sp.kron(sp.diags(lam), sp.identity(4))
    return abs(Lm @ T @ Lm.getH() - X).max()
gp = {}
for sub in itertools.product((0, 1), repeat=5):
    M = np.eye(4, dtype=complex)
    for j, b in enumerate(sub):
        if b: M = M @ GAM[j]
    gp[sub] = M
res = []
for conj in (0, 1):
    for refl in itertools.product((0, 1), repeat=3):
        for sub, gm in gp.items():
            for sign in (1, -1):
                e = test(0.37, refl, gm, conj, sign)
                if e < 1e-8: res.append((conj, refl, sub, sign, e))
for r_ in res:
    e2 = test(2.2, r_[1], gp[r_[2]], r_[0], r_[3])
    print(f'  U = Lambda (reflexion {r_[1]} x Gamma{[i+1 for i, b in enumerate(r_[2]) if b]}) {"K " if r_[0] else ""}: U H(th) U^-1 = {r_[3]:+d} H(-th), ecart {r_[4]:.1e} (th = 0.37), {e2:.1e} (th = 2.2)')
print(f'{len(res)} operateurs trouves')
