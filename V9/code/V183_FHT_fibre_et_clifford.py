# V183 — (1) fibre de spin sur deux niveaux (C0–C4) ; (2) Cl(5) covariant sur l'octet tordu (D-P1..D-P6). Critères : V182.
import numpy as np, itertools as it
from scipy.optimize import minimize
src = open("V180_FHT_spin_cube_monopole.py").read().split("# ---------- octet (cube) ----------")[0]
exec(src)                                   # ROT, Rx, Ry, build, outward_cycle, lift, comm, levels, flux_list
rng = np.random.default_rng(183)
# =============== (1) arbre à deux niveaux ===============
leafpos = [np.array(c, float) - 1.5 for c in it.product(range(4), repeat=3)]       # 4x4x4, centre en 0
L1pos = [np.array(s, float) for s in it.product([-1., 1.], repeat=3)]               # nœuds de niveau 1
par = [int(np.argmin([np.linalg.norm(p - q) for q in L1pos])) for p in leafpos]     # parent de chaque feuille
def leafH(twist):
    H = np.zeros((64, 64), complex)
    idx = {tuple(p): i for i, p in enumerate(leafpos)}
    for i, p in enumerate(leafpos):
        x = (p + 1.5).astype(int)
        for d in range(3):
            q = p.copy(); q[d] += 1
            if tuple(q) in idx:
                eta = [1, (-1)**x[0], (-1)**(x[0]+x[1])][d] if twist else 1
                t = 1.0 if par[i] == par[idx[tuple(q)]] else 0.8      # frères / cousins (poids génériques)
                H[i, idx[tuple(q)]] = H[idx[tuple(q)], i] = eta*t
    # potentiel radial générique (invariant par rotation) pour tuer les symétries accidentelles
    for i, p in enumerate(leafpos): H[i, i] = 0.13*np.dot(p, p)
    return H
def cubeH_scalar(pos, twist, t=1.0):
    n = len(pos); H = np.zeros((n, n), complex)
    for a, b in it.combinations(range(n), 2):
        dvec = pos[b] - pos[a]
        if np.isclose(np.sum(np.abs(dvec)), 2) and np.count_nonzero(np.abs(dvec) > 1e-9) == 1:
            H[a, b] = H[b, a] = t
    if twist:   # flux π sur les 6 faces : phases de Kawamoto–Smit sur le cube
        x = [((p + 1)/2).astype(int) for p in pos]
        for a, b in it.combinations(range(n), 2):
            if abs(H[a, b]) > 0:
                d = int(np.argmax(np.abs(pos[b] - pos[a]))); xa = x[a] if pos[b][d] > pos[a][d] else x[b]
                eta = [1, (-1)**xa[0], (-1)**(xa[0]+xa[1])][d]; H[a, b] *= eta; H[b, a] *= eta
    return H
def assemble(blocks, pos_of_block, Hfull):
    return dict(blocks=blocks, pos=pos_of_block, H=Hfull)
def solve_lift(sysd, R):
    blocks, pos, H = sysd["blocks"], sysd["pos"], sysd["H"]
    key = {tuple(np.round(p, 6)): b for b, p in enumerate(pos)}
    img = [key[tuple(np.round(R @ p, 6))] for p in pos]
    off = np.cumsum([0] + blocks)
    vars_ = []; 
    for b, d in enumerate(blocks):
        for i in range(d):
            for j in range(d): vars_.append((img[b], b, i, j))       # U[off[img]+i, off[b]+j]
    n = H.shape[0]; cols = []
    for (bi, b, i, j) in vars_:
        E = np.zeros((n, n), complex); E[off[bi]+i, off[b]+j] = 1
        cols.append((H @ E - E @ H).ravel())
    A = np.array(cols).T
    u, s, vt = np.linalg.svd(A, full_matrices=False)
    null = vt[s < 1e-8 * s[0]]
    if len(null) == 0: return None, 0
    v = null[0].conj() if False else null[0]
    U = np.zeros((n, n), complex)
    for c, (bi, b, i, j) in zip(v.conj(), vars_): U[off[bi]+i, off[b]+j] = c
    U = U / np.sqrt(np.real(np.trace(U.conj().T @ U))/n)
    return U, len(null)
def report(name, sysd, pred):
    Ux, kx = solve_lift(sysd, Rx); Uy, ky = solve_lift(sysd, Ry)
    nall = [solve_lift(sysd, R)[1] for R in ROT]
    if Ux is None or Uy is None or min(nall) == 0:
        print(" %-3s relevés : AUCUN pour %d/24 rotations   (prédit : %s)" % (name, sum(k == 0 for k in nall), pred)); return
    unit = np.allclose(Ux.conj().T @ Ux, np.eye(len(Ux)), atol=1e-8)
    C = Ux @ Uy @ np.linalg.inv(Ux) @ np.linalg.inv(Uy); lam = C[0, 0]
    sc = np.allclose(C, lam*np.eye(len(C)), atol=1e-7)
    w = np.linalg.eigvalsh(sysd["H"]); deg = []
    k = 0
    while k < len(w) and len(deg) < 6:
        l = k
        while l < len(w) and abs(w[l]-w[k]) < 1e-7: l += 1
        deg.append(l-k); k = l
    print(" %-3s relevés : 24/24 (dim noyau max %d), unitaire %s, commutateur %+.3f%+.3fi scalaire %s ; dégénérescences basses %s   (prédit : %s)"
          % (name, max(nall), unit, lam.real, lam.imag, sc, deg, pred))
def tree(conf):
    tw_leaf = conf in ("C1", "C2", "C3"); Hl = leafH(tw_leaf)
    blocks = [1]*64; pos = list(leafpos); Hparts = {}
    if conf in ("C0", "C1", "C4"):
        # niveau 1 scalaire, couplage uniforme w
        n = 64 + 8 + 1; H = np.zeros((n, n), complex); H[:64, :64] = Hl
        for i in range(64): H[i, 64+par[i]] = H[64+par[i], i] = 0.7
        H1 = cubeH_scalar(L1pos, twist=(conf == "C4"), t=0.9); H[64:72, 64:72] = H1
        for c in range(8): H[64+c, 64+c] = 0.31
        blocks += [1]*8; pos += L1pos
        if conf == "C4":     # racine = bloc sur le niveau propre −√3 (4 dim) du cube tordu de niveau 1
            w_, V = np.linalg.eigh(H1); Vs = V[:, :4]
            Hn = np.zeros((n+3, n+3), complex); Hn[:72, :72] = H[:72, :72]
            Hn[64:72, 72:76] = 0.6*Vs; Hn[72:76, 64:72] = 0.6*Vs.conj().T
            for k in range(4): Hn[72+k, 72+k] = 0.47
            H = Hn; blocks += [4]
        else:
            H[64:72, 72] = H[72, 64:72] = 0.6; H[72, 72] = 0.47; blocks += [1]
        pos += [np.zeros(3)]
        return dict(blocks=blocks, pos=pos, H=H)
    # C2, C3 : niveau 1 = blocs de 4 (isométrie sur le niveau propre −√3 du cube de frères tordu)
    W = []
    for c in range(8):
        kids = [i for i in range(64) if par[i] == c]
        Hc = Hl[np.ix_(kids, kids)] - np.diag(np.diag(Hl[np.ix_(kids, kids)]))
        w_, V = np.linalg.eigh(Hc); assert np.allclose(w_[:4], -np.sqrt(3)) and np.allclose(w_[4:], np.sqrt(3)), w_
        Wc = np.zeros((64, 4), complex); Wc[kids, :] = V[:, :4]; W.append(Wc)
    n = 64 + 32; H = np.zeros((n, n), complex); H[:64, :64] = Hl
    for c in range(8):
        H[:64, 64+4*c:68+4*c] = 0.7*W[c]; H[64+4*c:68+4*c, :64] = 0.7*W[c].conj().T
        H[64+4*c:68+4*c, 64+4*c:68+4*c] = 0.31*np.eye(4)
    # transport latéral induit (règle r1) : T_cc' = W_c† H_cousins W_c'
    for a, b in it.combinations(range(8), 2):
        if np.isclose(np.sum(np.abs(L1pos[a]-L1pos[b])), 2):
            T = W[a].conj().T @ Hl @ W[b]
            if conf == "C3":
                d = int(np.argmax(np.abs(L1pos[b]-L1pos[a]))); xa = ((np.minimum(L1pos[a], L1pos[b]) + 1)/2).astype(int)
                T = T * [1, (-1)**xa[0], (-1)**(xa[0]+xa[1])][d]          # twist ℤ/2 supplémentaire au niveau 1
            H[64+4*a:68+4*a, 64+4*b:68+4*b] = 1.5*T; H[64+4*b:68+4*b, 64+4*a:68+4*a] = 1.5*T.conj().T
    blocks += [4]*8; pos += L1pos
    # racine : bloc sur le niveau propre le plus bas (de dimension d) du niveau 1 seul
    H1 = H[64:96, 64:96]; w_, V = np.linalg.eigh(H1)
    d = 1
    while d < len(w_) and abs(w_[d]-w_[0]) < 1e-7: d += 1
    Hn = np.zeros((n+d, n+d), complex); Hn[:n, :n] = H
    Hn[64:96, n:n+d] = 0.6*V[:, :d]; Hn[n:n+d, 64:96] = 0.6*V[:, :d].conj().T; Hn[n:, n:] += 0.47*np.eye(d)
    print("   (%s : dimension du niveau le plus bas du niveau 1 = %d → bloc racine de %d)" % (conf, d, d))
    return dict(blocks=blocks + [d], pos=pos + [np.zeros(3)], H=Hn)
print("=== (1) fibre sur deux niveaux ===")
for conf, pred in [("C0", "+1"), ("C1", "aucun relevé"), ("C2", "−1"), ("C3", "aucun relevé"), ("C4", "−1")]:
    report(conf, tree(conf), pred)
# =============== (2) Cl(5) covariant sur l'octet ===============
print("=== (2) Cl(5) covariant dans End(C^8) ===")
cube = [np.array(c, float)*2-1 for c in it.product([0,1], repeat=3)]
ced = [(a, b) for a, b in it.combinations(range(8), 2) if np.sum(np.abs(cube[a]-cube[b])) == 2]
cfaces = [outward_cycle(cube, [i for i in range(8) if cube[i][ax] == s]) for ax in range(3) for s in (1, -1)]
def nullspace(rows, n, tol=1e-9):
    A = np.array(rows); u, s, vt = np.linalg.svd(A); r = np.sum(s > tol)
    return vt[r:].conj()
def herm_basis(vecs, shape):
    # base réelle d'un sous-espace complexe stable par †
    M = []
    for v in vecs:
        X = v.reshape(shape)
        for Y in ((X + X.conj().swapaxes(-1, -2))/2, (X - X.conj().swapaxes(-1, -2))/2j):
            M.append(Y.ravel())
    M = np.array(M); Mr = np.hstack([M.real, M.imag]); u, s, vt = np.linalg.svd(Mr, full_matrices=False)
    B = vt[s > 1e-9]; k = B.shape[1]//2
    return [ (b[:k] + 1j*b[k:]).reshape(shape) for b in B ]
for lab, m in [("tordu (π, m=3)", 3), ("non tordu", 0)]:
    H = build(cube, ced, cfaces, flux_list(6, 2*np.pi*m/6, m))
    Us = [lift(H, cube, R)[0] for R in ROT]
    I8 = np.eye(8)
    # A1 : U X U† = X
    rows = [ (np.kron(U, U.conj()) - np.eye(64)) for U in Us ]
    A1 = nullspace(np.vstack(rows), 64)
    # T1 : Ad(U) X_i = Σ_j R_ji X_j ; essai des deux conventions
    best = None
    for conv in ("R", "RT"):
        rows = []
        for U, R in zip(Us, ROT):
            Rm = R if conv == "R" else R.T
            AdU = np.kron(U, U.conj())                         # vec(U X U†) = (U ⊗ conj U) vec X (ordre ligne)
            rows.append(np.kron(np.eye(3), AdU) - np.kron(Rm.T, np.eye(64)))
        T1 = nullspace(np.vstack(rows), 192)
        if best is None or len(T1) > len(best[1]): best = (conv, T1)
    conv, T1 = best
    A1h = herm_basis(A1, (8, 8)); T1h = herm_basis(T1, (3, 8, 8))
    print(" %-15s dim A1 = %d, dim T1 = %d (convention %s) ; bases hermitiennes : %d, %d" % (lab, len(A1), len(T1), conv, len(A1h), len(T1h)))
    if len(A1h) < 2 or len(T1h) < 1:
        print("   → pas assez de scalaires ou de vecteurs : pas de Cl(5) covariant possible"); continue
    def mats(p):
        a4 = p[:len(A1h)]; a5 = p[len(A1h):2*len(A1h)]; b = p[2*len(A1h):]
        G4 = sum(c*X for c, X in zip(a4, A1h)); G5 = sum(c*X for c, X in zip(a5, A1h))
        Gv = sum(c*X for c, X in zip(b, T1h))
        return [Gv[0], Gv[1], Gv[2], G4, G5]
    def cost(p):
        G = mats(p); c = 0
        for i in range(5):
            for j in range(i, 5):
                A = G[i] @ G[j] + G[j] @ G[i] - (2*I8 if i == j else 0)
                c += np.sum(np.abs(A)**2)
        return c
    sols = []
    for trial in range(60):
        p0 = rng.normal(size=2*len(A1h)+len(T1h))
        r = minimize(cost, p0, method="BFGS", options=dict(gtol=1e-12, maxiter=4000))
        if r.fun < 1e-14: sols.append(mats(r.x))
    print("   solutions Clifford trouvées : %d / 60 essais" % len(sols), end="")
    if not sols:
        print(" — aucune"); continue
    print()
    types = []
    for G in sols:
        G6 = G[0] @ G[1] @ G[2] @ G[3] @ G[4]
        types.append(np.allclose(G6, G6[0, 0]*I8, atol=1e-6))
    print("   Γ6 scalaire pour toutes les solutions :", all(types), "(", sum(types), "/", len(types), ")")
    G = sols[0]
    # commutant de l'algèbre engendrée
    rows = [np.kron(X, I8) - np.kron(I8, X.T) for X in G]
    com = nullspace(np.vstack(rows), 64, tol=1e-5); print("   dimension du commutant (tolérance 1e-5, solutions précises à ~1e-8) :", len(com))
    # action des rotations sur le commutant (caractères par classe)
    Cb = np.array([c for c in com]).T                      # 64 x k
    P = np.linalg.pinv(Cb)
    cls = {}
    for U, R in zip(Us, ROT):
        tr = np.real(np.trace(P @ np.kron(U, U.conj()) @ Cb))
        key = (round(np.trace(R)), )
        cls.setdefault(key, []).append(round(tr, 6))
    print("   caractère de Ad sur le commutant, par trace de R (3:E, 0:C3, -1:C2/C2', 1:C4) :", {k: sorted(set(v)) for k, v in cls.items()})
    # D-P5 : candidat explicite
    G5c = H/np.sqrt(3)
    gam = []
    for ax in range(3):
        g = np.zeros_like(H)
        for a, b in ced:
            if abs(cube[a][ax]-cube[b][ax]) > 1: g[a, b] = H[a, b]; g[b, a] = H[b, a]
        gam.append(g)
    gp = [(g - G5c @ g @ G5c)/2 for g in gam]
    print("   D-P5 : γ⊥² =", [np.round(np.linalg.eigvalsh(x @ x), 4).tolist()[:2] for x in gp][:1],
          "; {γ⊥_i,γ⊥_j}=0 :", [np.allclose(gp[i]@gp[j]+gp[j]@gp[i], 0, atol=1e-9) for i, j in [(0,1),(0,2),(1,2)]],
          "; {γ⊥,Γ5}=0 :", all(np.allclose(x@G5c+G5c@x, 0) for x in gp))
