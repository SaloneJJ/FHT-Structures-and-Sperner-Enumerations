# Verification (Fable) : algebre reelle engendree par les cinq generateurs de l'octet a flux pi
import itertools, numpy as np
sites = list(itertools.product((0,1), repeat=3)); idx = {s:i for i,s in enumerate(sites)}
def hop(d, phase):
    M = np.zeros((8,8))
    for s in sites:
        t = list(s); t[d] ^= 1
        M[idx[tuple(t)], idx[s]] = phase(s)
    return M
a1 = hop(0, lambda s: 1); a2 = hop(1, lambda s: (-1)**s[0]); a3 = hop(2, lambda s: (-1)**(s[0]+s[1]))
g5 = np.diag([(-1)**sum(s) for s in sites]).astype(float)
g4 = g5 @ a1 @ a2 @ a3
G = [a1, a2, a3, g4, g5]
ok = all(np.allclose(X, X.T) and np.allclose(X@X, np.eye(8)) for X in G) and all(np.allclose(G[i]@G[j]+G[j]@G[i], 0) for i in range(5) for j in range(i))
# flux par face : produit des quatre sauts autour d'une plaquette
flux = [np.diag(a@b@a@b)[0] for a,b in ((a1,a2),(a1,a3),(a2,a3))]
w = a1@a2@a3@g4@g5
def span_dim(gens, field):
    basis = [np.eye(8, dtype=complex)]
    def rank(B):
        if field == 'R': V = np.array([np.concatenate([b.real.ravel(), b.imag.ravel()]) for b in B])
        else: V = np.array([b.ravel() for b in B])
        return np.linalg.matrix_rank(V, tol=1e-9)
    r = 1
    while True:
        new = basis + [b@g for b in basis for g in gens]
        # reduire
        V = np.array([np.concatenate([b.real.ravel(), b.imag.ravel()]) for b in new])
        u,s,vt = np.linalg.svd(V, full_matrices=False); k = (s>1e-9).sum()
        basis = [ (vt[i,:64] + 1j*vt[i,64:]).reshape(8,8) for i in range(k)]
        if k == r: return k, basis
        r = k
dR, BR = span_dim([X.astype(complex) for X in G], 'R')
dI, BI = span_dim([1j*X for X in G], 'R')
# commutant reel de l'algebre engendree par les Gamma dans M8(R)
A = []
for X in G:
    A.append(np.kron(np.eye(8), X) - np.kron(X.T, np.eye(8)))
u,s,vt = np.linalg.svd(np.vstack(A)); comm = (s<1e-9).sum() + (64-len(s))
# carres des elements imaginaires du commutant (structure H ?)
print("cinq involutions reelles symetriques anticommutantes :", ok, "| flux par face :", flux, "| produit =", round(w[0,0]), "x Id :", np.allclose(w, w[0,0]*np.eye(8)))
print("dim_R algebre engendree par Gamma_A :", dR, " (M2(H) = 16)")
print("dim_R algebre engendree par i Gamma_A :", dI, " (M4(C) = 32)")
print("dim_R commutant dans M8(R) :", comm, " (H = 4)")
