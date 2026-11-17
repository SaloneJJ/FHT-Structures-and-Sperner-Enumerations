"""Q4 : triplet a trois generations (H = H_32 (x) C^3), espace des D admissibles avec photon,
et fluctuations internes d'un D generique."""
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__)) + '/'

import contextlib, io, time
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(HERE + 'dspace.py').read().split("# V4 : commute with left rho(S4)")[0])
import numpy as np
t0 = time.time()
# 1-gen complex space : odd, order one, photon (complex-linear) -> basis B1c
X = left(1j * (e1 - es) + gprod(e2, u((1, 0, 2, 3))))
M1c = lin_constraints(lambda D: np.concatenate([order1(D), (D @ X - X @ D).ravel()]))
Bc, s = complex_basis(M1c)
print("1 generation : espace complexe (impair, ordre un, photon) dim_C =", len(Bc))
ng = 3
I3 = np.eye(ng)
Pg = np.kron(P, I3)                      # J = (J_32 (x) K) : X -> Pg conj(X) Pg
def Jop3(M): return Pg @ np.conj(M) @ Pg
cols = []
for B in Bc:
    for i in range(ng):
        for j in range(ng):
            E = np.zeros((ng, ng)); E[i, j] = 1
            cols.append(np.kron(B, E))
n = len(cols)
rows = []
for k in range(n):
    for part in (1, 1j):
        Bm = part * cols[k]
        rows.append(np.concatenate([(Bm - Bm.conj().T).ravel(), (Bm - Jop3(Bm)).ravel()]))
M = np.array(rows).T; Mr = np.vstack([M.real, M.imag])
u_, sv, vh = np.linalg.svd(Mr, full_matrices=False)
r = (sv > 1e-9 * sv[0]).sum(); ns = vh[r:]
print("3 generations : dim_R des D admissibles avec photon =", ns.shape[0], " (attendu 84)")
rng = np.random.default_rng(3)
coef = rng.normal(size=ns.shape[0]) @ ns
z = coef[0::2] + 1j * coef[1::2]
D0 = sum(z[k] * cols[k] for k in range(n))
L3 = [np.kron(left(c), I3) for c in AF_basis]
R3 = [Jop3(Lc) for Lc in L3]
assert np.abs(D0 - D0.conj().T).max() < 1e-9 and np.abs(D0 - Jop3(D0)).max() < 1e-9
ord1 = max(np.abs((D0 @ a - a @ D0) @ b - b @ (D0 @ a - a @ D0)).max() for a in L3[:6] for b in R3[:6])
print("D0 : autoadjoint, J'D0 = D0J', ordre un (echantillon) :", f"{ord1:.1e}")
forms = []
for a in L3:
    for b in L3:
        w = a @ (D0 @ b - b @ D0); forms.append(w + w.conj().T)
fl = [A + Jop3(A) for A in forms]
def vecR(Mm): return np.concatenate([Mm.real.ravel(), Mm.imag.ravel()])
F = np.array([vecR(A) for A in fl])
rf = np.linalg.matrix_rank(F, tol=1e-8)
print("fluctuations internes A + J'AJ'^-1 : dim_R =", rf, " (un doublet de Higgs : 4)")
ok = max(np.abs(A @ np.kron(gamma, I3) + np.kron(gamma, I3) @ A).max() for A in fl[:50])
print("fluctuations impaires :", f"{ok:.1e}", f"({time.time()-t0:.0f} s)")
