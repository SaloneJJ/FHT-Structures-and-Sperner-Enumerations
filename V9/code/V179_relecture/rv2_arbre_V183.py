# Relecture V183 (1) : arbre à deux niveaux, configurations C0–C4, construction indépendante.
# Poids génériques DIFFÉRENTS de ceux de V183, jauge aléatoire appliquée au twist, bases aléatoires dans les niveaux propres.
import numpy as np, itertools as it
from rv_common import ROTS, HX, HY
rng = np.random.default_rng(2026)

leaf = [np.array(p, float) - 1.5 for p in it.product(range(4), repeat=3)]
L1 = [np.array(s, float) for s in it.product([-1., 1.], repeat=3)]
par = [int(np.argmin([np.linalg.norm(p - q) for q in L1])) for p in leaf]
W = dict(fr=1.0, co=0.83, rad=0.17, vert=0.61, on1=0.29, hop1=0.93, root=0.57, onr=0.41, tr=1.3)

def ks_sign(xint, d):
    # signe de Kawamoto–Smit : flux π sur toutes les plaquettes
    return [1, (-1) ** xint[0], (-1) ** (xint[0] + xint[1])][d]

def leaves(twist):
    n = 64; H = np.zeros((n, n), complex); key = {tuple(p): i for i, p in enumerate(leaf)}
    for i, p in enumerate(leaf):
        x = (p + 1.5).astype(int)
        for d in range(3):
            q = p.copy(); q[d] += 1
            j = key.get(tuple(q))
            if j is None: continue
            t = W["fr"] if par[i] == par[j] else W["co"]
            s = ks_sign(x, d) if twist else 1
            H[i, j] = H[j, i] = s * t
        H[i, i] = W["rad"] * p @ p
    if twist:   # jauge aléatoire (ne change ni les flux ni l'existence de relevés)
        D = np.diag(np.exp(2j * np.pi * rng.random(n))); H = D @ H @ D.conj().T
    return H

def cube_L1(twist, t):
    H = np.zeros((8, 8), complex)
    for a, b in it.combinations(range(8), 2):
        dv = L1[b] - L1[a]
        if np.count_nonzero(dv) == 1:
            d = int(np.flatnonzero(dv)[0]); lo = ((np.minimum(L1[a], L1[b]) + 1) / 2).astype(int)
            s = ks_sign(lo, d) if twist else 1
            H[a, b] = H[b, a] = s * t
    return H

def eigblock(Hs, which="low"):
    w, V = np.linalg.eigh(Hs); d = 1
    while d < len(w) and abs(w[d] - w[0]) < 1e-7: d += 1
    Q = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))[0]
    return V[:, :d] @ Q, w[0], d

def build(conf):
    twl = conf in ("C1", "C2", "C3"); Hl = leaves(twl)
    if conf in ("C0", "C1", "C4"):
        n = 73; H = np.zeros((n, n), complex); H[:64, :64] = Hl
        for i in range(64): H[i, 64 + par[i]] = H[64 + par[i], i] = W["vert"]
        H1 = cube_L1(conf == "C4", W["hop1"]); H[64:72, 64:72] = H1 + W["on1"] * np.eye(8)
        blocks = [1] * 72
        if conf == "C4":
            Vs, _, d = eigblock(H1)                          # racine = niveau propre bas du cube tordu (dim 4)
            Hn = np.zeros((72 + d, 72 + d), complex); Hn[:72, :72] = H[:72, :72]
            Hn[64:72, 72:] = W["root"] * Vs; Hn[72:, 64:72] = W["root"] * Vs.conj().T; Hn[72:, 72:] = W["onr"] * np.eye(d)
            return Hn, blocks + [d]
        H[64:72, 72] = H[72, 64:72] = W["root"]; H[72, 72] = W["onr"]
        return H, blocks + [1]
    # C2, C3 : niveau 1 = blocs isométriques sur le niveau propre −√3 du cube de frères tordu
    Ws = []
    for c in range(8):
        kids = [i for i in range(64) if par[i] == c]
        Hc = Hl[np.ix_(kids, kids)].copy(); np.fill_diagonal(Hc, 0)
        Vc, e0, d = eigblock(Hc); assert d == 4 and np.isclose(e0, -np.sqrt(3) * W["fr"])
        Wc = np.zeros((64, 4), complex); Wc[kids] = Vc; Ws.append(Wc)
    n = 96; H = np.zeros((n, n), complex); H[:64, :64] = Hl
    for c in range(8):
        sl = slice(64 + 4 * c, 68 + 4 * c)
        H[:64, sl] = W["vert"] * Ws[c]; H[sl, :64] = W["vert"] * Ws[c].conj().T; H[sl, sl] = W["on1"] * np.eye(4)
    for a, b in it.combinations(range(8), 2):
        dv = L1[b] - L1[a]
        if np.count_nonzero(dv) == 1:
            T = W["tr"] * Ws[a].conj().T @ Hl @ Ws[b]
            if conf == "C3":
                d = int(np.flatnonzero(dv)[0]); lo = ((np.minimum(L1[a], L1[b]) + 1) / 2).astype(int); T = T * ks_sign(lo, d)
            H[64 + 4 * a:68 + 4 * a, 64 + 4 * b:68 + 4 * b] = T; H[64 + 4 * b:68 + 4 * b, 64 + 4 * a:68 + 4 * a] = T.conj().T
    Vs, _, d = eigblock(H[64:, 64:])
    Hn = np.zeros((n + d, n + d), complex); Hn[:n, :n] = H
    Hn[64:96, n:] = W["root"] * Vs; Hn[n:, 64:96] = W["root"] * Vs.conj().T; Hn[n:, n:] = W["onr"] * np.eye(d)
    return Hn, [1] * 64 + [4] * 8 + [d]

def lift_blocks(H, blocks, R):
    """U bloc-permutation (sites -> images par R), solution de H U = U H (système linéaire, noyau par SVD)."""
    pos = leaf + L1 + [np.zeros(3)]
    off = np.cumsum([0] + blocks); nb = len(blocks); n = H.shape[0]
    img = []
    for p in pos:
        q = R @ p; img.append(int(np.argmin([np.linalg.norm(q - r) if len(r) == 3 else 9 for r in pos])))
    cols = []; var = []
    for b in range(nb):
        for i in range(blocks[b]):
            for j in range(blocks[b]):
                r, c = off[img[b]] + i, off[b] + j; var.append((r, c))
                # vec(HE - EH) pour E = e_r e_c^T
                M = np.zeros((n, n), complex); M[:, c] += H[:, r]; M[r, :] -= H[c, :]
                cols.append(M.ravel())
    A = np.array(cols).T
    s, vt = np.linalg.svd(A, full_matrices=False)[1:]
    null = vt[s < 1e-8 * s[0]].conj()
    if len(null) == 0: return None, 0
    U = np.zeros((n, n), complex)
    for v, (r, c) in zip(null[0], var): U[r, c] = v
    U *= np.sqrt(n / np.real(np.trace(U.conj().T @ U)))
    return U, len(null)

for conf in ["C0", "C1", "C2", "C3", "C4"]:
    H, blocks = build(conf)
    res = [lift_blocks(H, blocks, R) for R in ROTS]
    nfail = sum(r[0] is None for r in res)
    if nfail:
        print(" %s : pas de relevé pour %d/24 rotations (dim H = %d, bloc racine %d)" % (conf, nfail, len(H), blocks[-1])); continue
    Ux = res[[k for k, R in enumerate(ROTS) if np.array_equal(R, HX)][0]][0]
    Uy = res[[k for k, R in enumerate(ROTS) if np.array_equal(R, HY)][0]][0]
    unit = all(np.allclose(U.conj().T @ U, np.eye(len(U)), atol=1e-7) for U, _ in res)
    C = Ux @ Uy @ Ux.conj().T @ Uy.conj().T; lam = C[0, 0]
    w = np.linalg.eigvalsh(H); degs = []; k = 0
    while k < len(w) and len(degs) < 6:
        l = k
        while l < len(w) and abs(w[l] - w[k]) < 1e-7: l += 1
        degs.append(l - k); k = l
    print(" %s : relevés 24/24, noyau max %d, unitaires %s, [Ux,Uy] = %+.4f%+.4fi scalaire %s, dégénérescences basses %s (bloc racine %d)"
          % (conf, max(r[1] for r in res), unit, lam.real, lam.imag, np.allclose(C, lam * np.eye(len(C)), atol=1e-7), degs, blocks[-1]))
