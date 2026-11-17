# Verification (Fable) : signe du carre de la dualite d'Alexander combinatoire
# pour la paire de squelettes (sk_k, sk_{n-3-k}) du simplexe a n sommets.
# D_K : H_k(K) -> H_{n-3-k}(K^v), c |-> projection harmonique de comp(w), dw = c.
import itertools, numpy as np
def subsets(n, m): return list(itertools.combinations(range(n), m))
def bmat(n, m):  # bord C_{m-1}(simplexes a m sommets) -> simplexes a m-1 sommets
    rows = subsets(n, m-1); cols = subsets(n, m); ri = {s:i for i,s in enumerate(rows)}
    B = np.zeros((len(rows), len(cols)))
    for j,s in enumerate(cols):
        for t in range(m):
            f = s[:t]+s[t+1:]
            B[ri[f], j] = (-1)**t
    return B
def sgn(perm):
    p = list(perm); s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s
def comp(n, m):  # sous-ensembles de taille m -> complementaires (taille n-m), avec signe
    cols = subsets(n, m); rows = subsets(n, n-m); ri = {s:i for i,s in enumerate(rows)}
    C = np.zeros((len(rows), len(cols)))
    for j,s in enumerate(cols):
        c = tuple(x for x in range(n) if x not in s)
        C[ri[c], j] = sgn(s+c)
    return C
def nullspace(A, tol=1e-9):
    u,s,vt = np.linalg.svd(A); r = (s>tol).sum(); return vt[r:].T
def D(n, k):
    # K = sk_k : faces de taille <= k+1 ; cycles Z_k = ker d sur les (k+1)-ensembles
    Zk = nullspace(bmat(n, k+1)) if k >= 0 else np.ones((1,1))
    d_fill = bmat(n, k+2)                       # C_{k+1}(Delta) -> C_k
    kd = n-3-k                                   # dimension duale ; faces de taille kd+1 = n-k-2
    Zd = nullspace(bmat(n, kd+1)) if kd >= 0 else np.ones((1,1))
    Cm = comp(n, k+2)
    W = np.linalg.lstsq(d_fill, Zk, rcond=None)[0]   # remplissages
    img = Cm @ W
    P = Zd @ Zd.T                                # projecteur orthogonal sur les cycles duaux
    return Zk, Zd, Zd.T @ (P @ img)              # matrice dans les bases orthonormees
for n in (4,5,6,7,8):
    for k in range(0, n-2):
        kd = n-3-k
        if kd < 0: continue
        Zk, Zd, M1 = D(n, k)
        Zd2, Zk2, M2 = D(n, kd)
        # changement de base orthonormee : Zd2 vs Zd, Zk2 vs Zk
        T1 = Zd2.T @ Zd; T2 = Zk.T @ Zk2
        comp2 = T2 @ M2 @ T1 @ M1
        lam = np.trace(comp2)/comp2.shape[0]
        err = np.abs(comp2 - lam*np.eye(comp2.shape[0])).max()
        print(f"n={n} k={k} dual={kd} dim={comp2.shape[0]:3d}  D'D = {lam:+.4f} Id (ecart {err:.1e})  decalage={kd-k:+d}  (k+1)(n-2-k) pair: {((k+1)*(n-2-k))%2==0}")
