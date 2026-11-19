"""Decomposition of real twisted group algebras R^eta[G] (A1, Def. def:twisted):
basis u_g, u_g u_h = (-1)^{e(g)e(h)} u_{gh}, e(g) = (1 - eta(g))/2.
Numerical Wedderburn decomposition: centre, eigen-blocks of a generic central element,
type of each block (R / C / H) by the signature of the regular trace form.
Usage: python3 twisted_decomp.py   (prints the checks and the S4 x S3 cases)"""
import itertools, numpy as np
from scipy.linalg import null_space

def sign(p):
    s, seen = 1, set()
    for i in range(len(p)):
        if i in seen: continue
        j, L = i, 0
        while j not in seen:
            seen.add(j); j = p[j]; L += 1
        if L % 2 == 0: s = -s
    return s

def compose(p, q):            # (p*q)(i) = p(q(i))
    return tuple(p[q[i]] for i in range(len(q)))

def group_closure(gens, n):
    e = tuple(range(n)); G = {e}; frontier = [e]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = compose(s, g)
                if h not in G: G.add(h); new.append(h)
        frontier = new
    return sorted(G)

def decompose(G, eta, tol=1e-8, seed=0):
    n = len(G); idx = {g: i for i, g in enumerate(G)}
    e = [0 if eta(g) == 1 else 1 for g in G]
    # left multiplication matrices L[g] : u_h -> u_g u_h
    L = np.zeros((n, n, n))
    for a, g in enumerate(G):
        for b, h in enumerate(G):
            L[a][idx[compose(g, h)], b] = (-1) ** (e[a] * e[b])
    R = np.zeros((n, n, n))            # right multiplication by u_h
    for b, h in enumerate(G):
        for a, g in enumerate(G):
            R[b][idx[compose(g, h)], a] = (-1) ** (e[a] * e[b])
    # centre: x with L_g x = R_g x for all g
    M = np.vstack([L[a] - R[a] for a in range(n)])
    MtM = M.T @ M                                    # n x n, same kernel as M
    wv, Vv = np.linalg.eigh(MtM)
    Z = Vv[:, wv < 1e-8 * max(1.0, wv.max())]
    rng = np.random.default_rng(seed)
    z = Z @ rng.standard_normal(Z.shape[1])
    Lz = np.tensordot(z, L, axes=1)          # left mult by z
    w, V = np.linalg.eig(Lz)
    # group eigenvalues
    blocks, used = [], np.zeros(n, bool)
    order = np.argsort(w.real + 1e-3 * w.imag)
    for i in order:
        if used[i]: continue
        lam = w[i]
        same = np.where(np.abs(w - lam) < 1e-6)[0]
        conj = np.where(np.abs(w - np.conj(lam)) < 1e-6)[0] if abs(lam.imag) > 1e-6 else np.array([], int)
        grp = np.union1d(same, conj); used[grp] = True
        blocks.append((lam, len(grp), abs(lam.imag) > 1e-6))
    out = []
    for lam, d, cplx in blocks:
        if cplx:
            m = int(round(np.sqrt(d / 2))); out.append(f"M{m}(C)" if m > 1 else "C")
            continue
        # real block: basis = real kernel of (Lz - lam)
        K = null_space(Lz - lam.real * np.eye(n), rcond=1e-9)
        # trace form B(x,y) = Tr(L_x L_y)
        Lx = [np.tensordot(K[:, k], L, axes=1) for k in range(K.shape[1])]
        B = np.array([[np.trace(a @ b) for b in Lx] for a in Lx])
        ev = np.linalg.eigvalsh((B + B.T) / 2)
        sig = int(np.sum(ev > tol) - np.sum(ev < -tol))
        dim = K.shape[1]
        if sig > 0:
            m = int(round(np.sqrt(dim))); out.append(f"M{m}(R)" if m > 1 else "R")
        else:
            m = int(round(np.sqrt(dim / 4))); out.append(f"M{m}(H)" if m > 1 else "H")
    return Z.shape[1], sorted(out, key=lambda s: (len(s), s))

def total_dim(parts):
    t = 0
    for p in parts:
        m = int(p[1:p.index('(')]) if p[0] == 'M' else 1
        k = p[p.index('(') + 1] if p[0] == 'M' else p[0]
        t += m * m * {'R': 1, 'C': 2, 'H': 4}[k]
    return t

if __name__ == '__main__':
    # S4 on 4 points (check: A_F = C + H + M3(C))
    S4 = group_closure([(1, 0, 2, 3), (1, 2, 3, 0)], 4)
    print("S4, sgn :", decompose(S4, sign))
    S5 = group_closure([(1, 0, 2, 3, 4), (1, 2, 3, 4, 0)], 5)
    print("S5, sgn :", decompose(S5, sign))
    # S4 x S3 on 4+3 points
    G = group_closure([(1, 0, 2, 3, 4, 5, 6), (1, 2, 3, 0, 4, 5, 6),
                       (0, 1, 2, 3, 5, 4, 6), (0, 1, 2, 3, 5, 6, 4)], 7)
    s4 = lambda g: sign(g[:4]); s3 = lambda g: sign(tuple(x - 4 for x in g[4:]))
    cases = {"sgn_V = sgn4 x sgn3": lambda g: sign(g),
             "sgn4 x 1": s4, "1 x sgn3": s3, "1 (untwisted)": lambda g: 1}
    print("|S4 x S3| =", len(G))
    for name, eta in cases.items():
        zc, parts = decompose(G, eta)
        has_AF = all(any(p == q for p in parts) for q in ("C", "H", "M3(C)"))
        print(f"S4xS3, eta = {name}: dim Z = {zc}; {' + '.join(parts)}; dim = {total_dim(parts)}; "
              f"sector A_F: {'YES' if has_AF else 'no'}")

def s4xc2():
    G = group_closure([(1, 0, 2, 3, 4, 5), (1, 2, 3, 0, 4, 5), (0, 1, 2, 3, 5, 4)], 6)
    s4 = lambda g: sign(g[:4]); c2 = lambda g: 1 if g[4] == 4 else -1
    for name, eta in {"sgn_V = sgn4 x sgnC2": lambda g: sign(g), "sgn4 x 1": s4,
                      "1 x sgnC2": c2}.items():
        zc, parts = decompose(G, eta)
        print(f"S4xC2, eta = {name}: {' + '.join(parts)}; dim = {total_dim(parts)}")
