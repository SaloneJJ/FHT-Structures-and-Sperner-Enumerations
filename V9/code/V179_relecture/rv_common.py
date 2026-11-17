# Outils indépendants du relecteur (aucune fonction des scripts relus n'est importée).
import numpy as np, itertools as it

def rotations():
    """Les 24 rotations propres du cube (matrices de permutation signées, det = +1)."""
    out = []
    for p in it.permutations(range(3)):
        P = np.eye(3)[list(p)]
        for s in it.product([1, -1], repeat=3):
            M = np.diag(s) @ P
            if np.linalg.det(M) > 0:
                out.append(M.astype(int))
    assert len(out) == 24
    return out

ROTS = rotations()
HX = np.diag([1, -1, -1]); HY = np.diag([-1, 1, -1])   # demi-tours x et y

def site_perm(pts, R, tol=1e-6):
    pts = np.asarray(pts, float)
    img = pts @ np.asarray(R, float).T
    perm = []
    for q in img:
        d = np.linalg.norm(pts - q, axis=1); j = int(np.argmin(d))
        if d[j] > tol: return None
        perm.append(j)
    return perm

def lift_diag(H, pts, R, tol=1e-9):
    """Relevé U = D P (D diagonale unitaire) avec U H U^† = H, obtenu comme noyau du système LINÉAIRE
    H D P - D P H = 0 (inconnues : les n coefficients de D). Retourne (U, dim du noyau)."""
    n = len(H); perm = site_perm(pts, R)
    if perm is None: return None, 0
    P = np.zeros((n, n)); P[perm, np.arange(n)] = 1          # P e_i = e_perm(i)
    M = np.zeros((n * n, n), complex)
    for k in range(n):
        E = np.zeros((n, n), complex); E[k, k] = 1
        M[:, k] = (H @ E @ P - E @ P @ H).ravel()
    u, s, vt = np.linalg.svd(M)
    null = vt[np.sum(s > tol * max(1, s[0])):].conj()
    if len(null) == 0: return None, 0
    d = null[0]; d = d / np.abs(d).max()
    if not np.allclose(np.abs(d), 1, atol=1e-7): return None, len(null)
    U = np.diag(d) @ P
    assert np.allclose(U @ H @ U.conj().T, H, atol=1e-8)
    return U, len(null)

def commutator_class(Ua, Ub):
    C = Ua @ Ub @ Ua.conj().T @ Ub.conj().T
    lam = C[0, 0]
    return lam, np.allclose(C, lam * np.eye(len(C)), atol=1e-8)

def level_structure(H, Us, tol=1e-7):
    """Niveaux : (E, dégénérescence, Σ_g |χ(g)|² / |G|) — ce dernier = Σ des multiplicités² (1 si irréductible)."""
    w, V = np.linalg.eigh(H); out = []; k = 0
    while k < len(w):
        l = k
        while l < len(w) and abs(w[l] - w[k]) < tol: l += 1
        Vs = V[:, k:l]
        s = sum(abs(np.trace(Vs.conj().T @ U @ Vs)) ** 2 for U in Us) / len(Us)
        out.append((round(float(w[k]), 4), l - k, round(float(s), 3))); k = l
    return out

def peierls(edges, faces_oriented, flux, n_edges):
    """Champ de jauge A sur les arêtes orientées (a<b) tel que la circulation sur chaque face orientée = flux.
    faces_oriented : listes de sommets dans l'ordre (sens sortant)."""
    eid = {e: k for k, e in enumerate(edges)}
    M = np.zeros((len(faces_oriented), n_edges))
    for f, cyc in enumerate(faces_oriented):
        for a, b in zip(cyc, cyc[1:] + cyc[:1]):
            if (a, b) in eid: M[f, eid[(a, b)]] += 1
            else: M[f, eid[(b, a)]] -= 1
    A = np.linalg.lstsq(M, flux, rcond=None)[0]
    assert np.allclose(M @ A, flux, atol=1e-9), "flux non réalisable"
    return A

def orient_out(pts, idx, center=np.zeros(3)):
    P = np.array([pts[i] for i in idx], float); c = P.mean(0); nrm = c - center
    # tri angulaire autour de la normale sortante (méthode : angle dans une base orthonormée du plan)
    e1 = P[0] - c; e1 /= np.linalg.norm(e1); e2 = np.cross(nrm / np.linalg.norm(nrm), e1)
    ang = np.arctan2((P - c) @ e2, (P - c) @ e1)
    return [idx[k] for k in np.argsort(ang)]

def ham_from_A(n, edges, A, amp=None):
    H = np.zeros((n, n), complex)
    for k, (a, b) in enumerate(edges):
        t = 1.0 if amp is None else amp[k]
        H[a, b] += t * np.exp(1j * A[k]); H[b, a] += t * np.exp(-1j * A[k])
    return H

def nullspace(M, tol=1e-9):
    if M.shape[0] < M.shape[1]:
        M = np.vstack([M, np.zeros((M.shape[1] - M.shape[0], M.shape[1]))])
    u, s, vt = np.linalg.svd(M, full_matrices=False)
    r = np.sum(s > tol * max(1.0, s[0] if len(s) else 1.0))
    return vt[r:].conj()
