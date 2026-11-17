# Relecture V183 (2) : Cl(5) covariant dans End(C^8) de l'octet tordu ; argument « Γ6 scalaire » ; lien avec les α_d de V193.
import numpy as np, itertools as it
from scipy.optimize import least_squares
from rv_common import ROTS, lift_diag, nullspace
rng = np.random.default_rng(5)

# octet en jauge de Kawamoto–Smit (bits x ∈ {0,1}^3), positions ±1
X = [np.array(c) for c in it.product([0, 1], repeat=3)]; pos = [2. * x - 1 for x in X]
def ks(d, x): return [1, (-1) ** x[0], (-1) ** (x[0] + x[1])][d]
def hop(d, twist=True, phase=False):
    g = np.zeros((8, 8), complex)
    for a, x in enumerate(X):
        if x[d] == 0:
            y = x.copy(); y[d] = 1; b = X.index(list(y)) if False else [tuple(z) for z in X].index(tuple(y))
            s = ks(d, x) if twist else 1
            if phase: g[b, a] = 1j * s; g[a, b] = -1j * s     # « i γ_d Z_d »-type (saut de cousin dérivé)
            else: g[a, b] = g[b, a] = s
    return g
H = sum(hop(d) for d in range(3)); H0 = sum(hop(d, twist=False) for d in range(3))
I8 = np.eye(8)

def ad_data(Hc):
    Us = [lift_diag(Hc, pos, R)[0] for R in ROTS]
    A1 = nullspace(np.vstack([np.kron(U, U.conj()) - np.eye(64) for U in Us]))
    T1 = nullspace(np.vstack([np.kron(np.eye(3), np.kron(U, U.conj())) - np.kron(R.T, np.eye(64)) for U, R in zip(Us, ROTS)]))
    return Us, [a.reshape(8, 8) for a in A1], [t.reshape(3, 8, 8) for t in T1]

for lab, Hc in [("tordu", H), ("non tordu", H0)]:
    Us, A1, T1 = ad_data(Hc)
    comm = all(np.allclose(a @ b, b @ a) for a in A1 for b in A1)
    print("%-9s : dim A1 = %d, dim T1 = %d, A1 commutative : %s" % (lab, len(A1), len(T1), comm))

Us, A1, T1 = ad_data(H)
# --- base adaptée : A1 = 1_4 ⊗ M2 ; Γ5 = H/√3 ; Γ4 ∈ A1 hermitien anticommutant avec Γ5 ---
G5 = H / np.sqrt(3)
def herm(B):
    out = []
    for X_ in B:
        for Y in ((X_ + X_.conj().swapaxes(-1, -2)) / 2, (X_ - X_.conj().swapaxes(-1, -2)) / 2j):
            if np.linalg.norm(Y) > 1e-9: out.append(Y)
    M = np.array([y.ravel() for y in out]); M = np.hstack([M.real, M.imag])
    u, s, vt = np.linalg.svd(M, full_matrices=False); vt = vt[s > 1e-9]
    k = vt.shape[1] // 2; return [(v[:k] + 1j * v[k:]).reshape(B[0].shape) for v in vt]
A1h = herm(A1); T1h = herm(T1)
print("bases hermitiennes réelles : A1 %d, T1 %d ; H/√3 ∈ A1 : %s" % (len(A1h), len(T1h),
      np.linalg.matrix_rank(np.array([a.ravel() for a in A1h] + [G5.ravel()]), tol=1e-8) == len(A1h)))

# Recherche systématique : Γ4, Γ5 ∈ A1h, Γ1..3 ∈ T1h, contraintes de Clifford (moindres carrés, 200 départs)
def mats(p):
    a4, a5, b = p[:4], p[4:8], p[8:]
    G4 = sum(c * x for c, x in zip(a4, A1h)); G5_ = sum(c * x for c, x in zip(a5, A1h)); V = sum(c * x for c, x in zip(b, T1h))
    return [V[0], V[1], V[2], G4, G5_]
def res(p):
    G = mats(p); r = []
    for i in range(5):
        for j in range(i, 5):
            A = G[i] @ G[j] + G[j] @ G[i] - (2 * I8 if i == j else 0); r += [A.real.ravel(), A.imag.ravel()]
    return np.concatenate(r)
sols = []
for t in range(200):
    r = least_squares(res, rng.normal(size=8 + len(T1h)), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
    if np.linalg.norm(r.fun) < 1e-10: sols.append(mats(r.x))
g6 = []; inv123 = []
for G in sols:
    G6 = G[0] @ G[1] @ G[2] @ G[3] @ G[4]; g6.append(G6[0, 0] if np.allclose(G6, G6[0, 0] * I8, atol=1e-7) else None)
    P = G[0] @ G[1] @ G[2]; inv123.append(all(np.allclose(U @ P @ U.conj().T, P, atol=1e-7) for U in Us))
print("solutions Cl(5) covariantes : %d / 200 départs ; Γ6 = c·1 pour toutes : %s ; valeurs de c : %s ; Γ1Γ2Γ3 invariant (A1) : %s"
      % (len(sols), all(c is not None for c in g6), sorted({complex(np.round(c, 6)) for c in g6 if c is not None}, key=lambda z: z.real), all(inv123)))
# même recherche sans twist
Us0, A10, T10 = ad_data(H0); A1h0 = herm(A10); T1h0 = herm(T10)
def res0(p):
    a4, a5, b = p[:4], p[4:8], p[8:]
    G = [None] * 5; V = sum(c * x for c, x in zip(b, T1h0)); G[:3] = V
    G[3] = sum(c * x for c, x in zip(a4, A1h0)); G[4] = sum(c * x for c, x in zip(a5, A1h0))
    r = []
    for i in range(5):
        for j in range(i, 5):
            A = G[i] @ G[j] + G[j] @ G[i] - (2 * I8 if i == j else 0); r += [A.real.ravel(), A.imag.ravel()]
    return np.concatenate(r)
best = min(np.linalg.norm(least_squares(res0, rng.normal(size=8 + len(T1h0)), max_nfev=3000).fun) for t in range(60))
print("non tordu : meilleur résidu de Clifford sur 60 départs = %.3f (0 si solution)" % best)

# Commutant et action sur le commutant (doublet de goût)
G = sols[0]
com = nullspace(np.vstack([np.kron(g, I8) - np.kron(I8, g.T) for g in G]), tol=1e-7)
Cb = np.array([c for c in com]).T; Pp = np.linalg.pinv(Cb)
chars = {}
for U, R in zip(Us, ROTS):
    tr = np.real(np.trace(Pp @ np.kron(U, U.conj()) @ Cb)); ang = np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1)))
    lab = {0: "E", 120: "C3", 90: "C4", 180: "C2"}[int(round(ang))]
    if lab == "C2": lab = "C2(face)" if np.count_nonzero(np.diag(R) == 1) == 1 and np.count_nonzero(R - np.diag(np.diag(R))) == 0 else "C2'(arête)"
    chars.setdefault(lab, set()).add(round(tr, 6))
print("commutant : dim %d ; caractère de Ad : %s  (End(E) attendu : E 4, C3 1, C2 4, C4 0, C2' 0)" % (len(com), chars))

# Γ5 = H/√3 : famille des solutions ; Γ4 dans A1 anticommutant avec H
# D-P5 : γ⊥ = partie de γ_d anticommutant avec Γ5 ; représentation des γ_d
gam = [hop(d) for d in range(3)]
gp = [g - G5 / np.sqrt(3) for g in gam]
print("γ⊥ : Σ_d γ⊥_d = 0 : %s ; (γ⊥)² = %s·1 : %s" % (np.allclose(sum(gp), 0), round(np.real(np.trace(gp[0] @ gp[0])) / 8, 4), np.allclose(gp[0] @ gp[0], gp[0] @ gp[0][0, 0] if False else (2/3) * I8)))
# transformation des γ_d : signe sous un demi-tour x (e_y -> -e_y)
Ux = [U for U, R in zip(Us, ROTS) if np.array_equal(R, np.diag([1, -1, -1]))][0]
print("demi-tour x : U γ_y U† = %+.0f γ_y  (vecteur T1 : −1 ; permutation des axes A1+E : +1)" % np.real(np.trace(gam[1] @ Ux @ gam[1] @ Ux.conj().T) / 8))

# Lien V193 : α_d (dérivée des sauts de cousins en K) = i·(…) ; sont-ils dans T1, et complétables par un Γ4 ∈ A1 ?
alpha = [hop(d, phase=True) for d in range(3)]
alpha = [(a + a.conj().T) / 2 for a in alpha]
print("α_d : carré 1 : %s ; anticommutation mutuelle : %s ; {α_d, Γ5} = 0 : %s ; recouvrement avec γ_d : %s"
      % (all(np.allclose(a @ a, I8) for a in alpha), all(np.allclose(alpha[i] @ alpha[j] + alpha[j] @ alpha[i], 0) for i, j in [(0, 1), (0, 2), (1, 2)]),
         all(np.allclose(a @ G5 + G5 @ a, 0) for a in alpha), [round(abs(np.trace(alpha[d] @ gam[d])) / 8, 6) for d in range(3)]))
vT = np.concatenate([a.ravel() for a in alpha]); BT = np.array([t.ravel() for t in T1]).T
print("α dans T1 : résidu %.1e" % (np.linalg.norm(vT - BT @ np.linalg.lstsq(BT, vT, rcond=None)[0]) / np.linalg.norm(vT)))
# Γ4 = Σ c_k A1h_k avec {Γ4, Γ5} = 0 et {Γ4, α_d} = 0 (système linéaire réel)
rows = []
for Y in [G5] + alpha:
    cols = [(a @ Y + Y @ a).ravel() for a in A1h]; M = np.array(cols).T; rows.append(np.vstack([M.real, M.imag]))
ns = nullspace(np.vstack(rows))
print("Γ4 ∈ A1 anticommutant avec Γ5 et les α_d : dimension réelle des solutions = %d" % len(ns))
if len(ns):
    G4 = sum(c * a for c, a in zip(np.real(ns[0]), A1h)); G4 = G4 / np.sqrt(np.real(np.trace(G4 @ G4)) / 8)
    G6 = alpha[0] @ alpha[1] @ alpha[2] @ G4 @ G5
    print("  Γ4² = 1 : %s ; Γ6 = α1α2α3Γ4Γ5 scalaire : %s (valeur %s)" % (np.allclose(G4 @ G4, I8), np.allclose(G6, G6[0, 0] * I8), np.round(G6[0, 0], 6)))
