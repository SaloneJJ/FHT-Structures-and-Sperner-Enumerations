# V180 — spin sur le cube / germe K26 / tetraedre : flux uniforme sortant (monopole au centre), relevés projectifs des rotations.
import numpy as np, itertools as it
np.set_printoptions(precision=4, suppress=True)
# rotations du cube (O) : matrices de permutation signées de det +1
ROT = [np.array(M) for M in (np.diag(s) @ np.eye(3)[list(p)] for p in it.permutations(range(3)) for s in it.product([1,-1], repeat=3)) if round(np.linalg.det(M)) == 1]
assert len(ROT) == 24
Rx, Ry = np.diag([1,-1,-1]), np.diag([-1,1,-1])
def build(sites, edges, faces, flux):
    """sites: positions ; edges: paires (a,b) ; faces: cycles orientés (sortants) ; flux: flux par face (déjà corrigés pour sommer à 0)."""
    n, E = len(sites), len(edges)
    eid = {e: k for k, e in enumerate(edges)}
    Mf = np.zeros((len(faces), E))
    for f, cyc in enumerate(faces):
        for a, b in zip(cyc, cyc[1:] + cyc[:1]):
            if (a, b) in eid: Mf[f, eid[(a, b)]] += 1
            else: Mf[f, eid[(b, a)]] -= 1
    A, *_ = np.linalg.lstsq(Mf, flux, rcond=None)
    assert np.allclose(Mf @ A, flux, atol=1e-9), "flux incompatible"
    H = np.zeros((n, n), complex)
    for (a, b), t in zip(edges, A): H[a, b] = np.exp(1j*t); H[b, a] = np.exp(-1j*t)
    return H
def outward_cycle(pts, idx, center=np.zeros(3)):
    # ordonne les sommets d'une face autour de son centre, sens direct vu de l'extérieur
    P = np.array([pts[i] for i in idx]); c = P.mean(0); nrm = c - center; nrm /= np.linalg.norm(nrm)
    u = P[0] - c; u /= np.linalg.norm(u); v = np.cross(nrm, u)
    ang = [np.arctan2((p-c) @ v, (p-c) @ u) for p in P]
    return [idx[k] for k in np.argsort(ang)]
def lift(H, pts, R):
    key = {tuple(np.round(p, 6)): i for i, p in enumerate(pts)}
    perm = [key[tuple(np.round(R @ p, 6))] for p in pts]       # site i -> perm[i]
    P = np.zeros_like(H); 
    for i, j in enumerate(perm): P[j, i] = 1
    Hp = P @ H @ P.T
    n = len(pts); D = np.zeros(n, complex); D[0] = 1; todo = [0]
    while todo:
        a = todo.pop()
        for b in np.nonzero(np.abs(H[a]) > 1e-9)[0]:
            if D[b] == 0:
                D[b] = np.conj(H[a, b] / (D[a] * Hp[a, b])); todo.append(b)
    for b in range(n):            # composantes non connexes : phase libre
        if D[b] == 0:
            D[b] = 1; todo = [b]
            while todo:
                a = todo.pop()
                for c in np.nonzero(np.abs(H[a]) > 1e-9)[0]:
                    if D[c] == 0: D[c] = np.conj(H[a, c] / (D[a] * Hp[a, c])); todo.append(c)
    U = np.diag(D) @ P
    return U, np.allclose(U @ H @ U.conj().T, H, atol=1e-9)
def comm(H, pts, A=Rx, B=Ry):
    Ua, oka = lift(H, pts, A); Ub, okb = lift(H, pts, B)
    C = Ua @ Ub @ Ua.conj().T @ Ub.conj().T
    lam = C[0, 0]; return lam, oka and okb and np.allclose(C, lam*np.eye(len(H)), atol=1e-9)
def levels(H, pts, group):
    w, V = np.linalg.eigh(H); out = []; k = 0
    Us = [lift(H, pts, R)[0] for R in group]
    while k < len(w):
        l = k
        while l < len(w) and abs(w[l]-w[k]) < 1e-7: l += 1
        Vs = V[:, k:l]; s = sum(abs(np.trace(Vs.conj().T @ U @ Vs))**2 for U in Us)/len(group)
        out.append((round(w[k], 4), l-k, round(s, 3))); k = l
    return out
def flux_list(nf, phi, m_target):
    f = np.full(nf, phi); f[0] -= 2*np.pi*m_target   # somme exacte nulle (mod 2π, une face absorbe 2πm)
    return f
# ---------- octet (cube) ----------
cube = [np.array(c, float)*2-1 for c in it.product([0,1], repeat=3)]
ced = [(a, b) for a, b in it.combinations(range(8), 2) if np.sum(np.abs(cube[a]-cube[b])) == 2]
cfaces = [outward_cycle(cube, [i for i in range(8) if cube[i][ax] == s]) for ax in range(3) for s in (1, -1)]
print("=== Q1/Q2 : octet, flux uniforme sortant phi = 2 pi m / 6 ===")
for m in range(6):
    H = build(cube, ced, cfaces, flux_list(6, 2*np.pi*m/6, m))
    lam, ok = comm(H, cube)
    print(" m=%d  commutateur = %+.3f%+.3fi (scalaire: %s)   prédit %+d" % (m, lam.real, lam.imag, ok, (-1)**m))
    if m in (1, 3): print("   niveaux (E, dégénérescence, Σ|χ|²/24) :", levels(H, cube, ROT))
H3 = build(cube, ced, cfaces, flux_list(6, np.pi, 3))
print(" ℤ/2 (m=3) : sauts réels ?", np.allclose(H3.imag, 0) or "jauge complexe (classe ℤ/2 équivalente)")
# Q3 : triplet des sauts par direction
print("=== Q3 : opérateurs de saut γ_i (flux π) ===")
G = []
for ax in range(3):
    g = np.zeros_like(H3)
    for a, b in ced:
        if abs(cube[a][ax]-cube[b][ax]) > 1: g[a, b] = H3[a, b]; g[b, a] = H3[b, a]
    G.append(g)
print(" γ_i² = 1 :", [np.allclose(g@g, np.eye(8)) for g in G], " {γ_i,γ_j} = 0 :", [np.allclose(G[i]@G[j]+G[j]@G[i], 0) for i, j in [(0,1),(0,2),(1,2)]])
cov = True
for R in ROT:
    U, _ = lift(H3, cube, R)
    for i in range(3):
        j = int(np.argmax(np.abs(R[:, i]))); X = U @ G[i] @ U.conj().T
        s = np.vdot(G[j].ravel(), X.ravel())/8
        cov &= abs(abs(s)-1) < 1e-9 and np.allclose(X, s*G[j], atol=1e-9)
print(" covariance U_g γ_i U_g† = s γ_{g(i)} pour les 24 rotations :", cov)
# ---------- germe K26 (coquille) ----------
shell = [np.array(p, float) for p in it.product([-1,0,1], repeat=3) if any(p)]
sed = [(a, b) for a, b in it.combinations(range(26), 2) if np.isclose(np.linalg.norm(shell[a]-shell[b]), 1)]
sq = []
for ax in range(3):
    for s in (1, -1):
        for u0, v0 in it.product([-1, 0], repeat=2):
            o = [k for k in range(3) if k != ax]; idx = []
            for du, dv in [(0,0),(1,0),(1,1),(0,1)]:
                p = np.zeros(3); p[ax] = s; p[o[0]] = u0+du; p[o[1]] = v0+dv
                idx.append(next(i for i, q in enumerate(shell) if np.allclose(q, p)))
            sq.append(outward_cycle(shell, idx))
print("=== Q4 : germe K26 : %d sites, %d arêtes, %d carrés (Euler %d) ===" % (len(shell), len(sed), len(sq), len(shell)-len(sed)+len(sq)))
for lab, phi, m in [("ℤ/2, π par carré (m=12)", np.pi, 12), ("U(1), m=1", 2*np.pi/24, 1), ("U(1), m=3", 2*np.pi*3/24, 3), ("U(1), m=2", 2*np.pi*2/24, 2)]:
    H = build(shell, sed, sq, flux_list(24, phi, m)); lam, ok = comm(H, shell)
    print(" %-26s commutateur = %+.3f%+.3fi (scalaire: %s)" % (lab, lam.real, lam.imag, ok))
    if m == 1: print("   niveaux bas :", levels(H, shell, ROT)[:3])
# ---------- tétraèdre R (twist de A_F) ----------
R4 = [p for p in cube if np.prod(p) > 0]
ted = list(it.combinations(range(4), 2)); tfaces = [outward_cycle(R4, list(f)) for f in it.combinations(range(4), 3)]
TET = [R for R in ROT if all(any(np.allclose(R @ p, q) for q in R4) for p in R4)]
print("=== Q5 : tétraèdre R, groupe T (%d éléments) ===" % len(TET))
for lab, phi, m in [("π/2 par face (m=1)", np.pi/2, 1), ("0", 0, 0), ("π (m=2)", np.pi, 2)]:
    H = build(R4, ted, tfaces, flux_list(4, phi, m)); lam, ok = comm(H, R4)
    print(" %-18s commutateur = %+.3f (scalaire: %s) ; niveaux : %s" % (lab, lam.real, ok, levels(H, R4, TET)))
