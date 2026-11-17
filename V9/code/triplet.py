import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__)) + '/'
import numpy as np, itertools
np.set_printoptions(linewidth=200)
TOL = 1e-9

# ---------- S4 ----------
G = list(itertools.permutations(range(4)))
idx = {g: k for k, g in enumerate(G)}
def mul(g, h):  # (g h)(x) = g(h(x))
    return tuple(g[h[x]] for x in range(4))
def inv(g):
    r = [0]*4
    for x in range(4): r[g[x]] = x
    return tuple(r)
def sgn(g):
    s = 1
    for i in range(4):
        for j in range(i+1, 4):
            if g[i] > g[j]: s = -s
    return s
def permmat(p, n):
    M = np.zeros((n, n))
    for x in range(n): M[p[x], x] = 1
    return M
def onb_sumzero(n):
    A = np.eye(n)[:, :n-1] - np.eye(n)[:, [n-1]]
    Q, _ = np.linalg.qr(A)
    return Q
Q4, Q3 = onb_sumzero(4), onb_sumzero(3)
pairings = [frozenset([frozenset([0,1]),frozenset([2,3])]),
            frozenset([frozenset([0,2]),frozenset([1,3])]),
            frozenset([frozenset([0,3]),frozenset([1,2])])]
def s3image(g):
    out = []
    for P in pairings:
        img = frozenset(frozenset(g[x] for x in b) for b in P)
        out.append(pairings.index(img))
    return tuple(out)
rho = {
 '1':   {g: np.eye(1) for g in G},
 's':   {g: np.eye(1)*sgn(g) for g in G},
 '2':   {g: Q3.T @ permmat(s3image(g), 3) @ Q3 for g in G},
 '3':   {g: Q4.T @ permmat(g, 4) @ Q4 for g in G},
}
rho['3p'] = {g: sgn(g)*rho['3'][g] for g in G}
dim = {k: rho[k][G[0]].shape[0] for k in rho}
IRR = ['1', 's', '2', '3', '3p']

# checks: homomorphisms, orthogonality, characters
for k in IRR:
    for g in G:
        assert np.allclose(rho[k][g] @ rho[k][g].T, np.eye(dim[k]))
        for h in G:
            assert np.allclose(rho[k][mul(g, h)], rho[k][g] @ rho[k][h])
chi = {k: np.array([np.trace(rho[k][g]) for g in G]) for k in IRR}
gram = np.array([[chi[a] @ chi[b] / 24 for b in IRR] for a in IRR])
assert np.allclose(gram, np.eye(5))
t12 = (1, 0, 2, 3)
assert np.isclose(chi['3'][idx[t12]], 1.0)   # convention : chi_3((12)) = +1
print("irreps OK, dims", dim)

# ---------- group algebra, A_F ----------
def ebar(g): return 0 if sgn(g) == 1 else 1
def u(g):  # u_g = i^{ebar(g)} g  as a 24-vector of complex coefficients
    v = np.zeros(24, complex); v[idx[g]] = 1j**ebar(g); return v
AF_basis = [u(g) for g in G]                       # real basis of A_F
def star(c):                                        # (sum c_g g)^* = sum conj(c_g) g^{-1}
    v = np.zeros(24, complex)
    for g in G: v[idx[inv(g)]] = np.conj(c[idx[g]])
    return v
def sigma(c):
    return np.array([np.conj(c[idx[g]])*sgn(g) for g in G])
for c in AF_basis: assert np.allclose(sigma(c), c)
def rep(k, c):
    return sum(c[idx[g]]*rho[k][g] for g in G)

# ---------- Hilbert space ----------
Nplus = [('1','1'), ('1','s'), ('2','1'), ('2','3'), ('3','1'), ('3','s')]
names_plus = ['nubarR', 'ebarR', 'L', 'Q', 'ubarR', 'dbarR']
comps = [(i, j, +1, n) for (i, j), n in zip(Nplus, names_plus)]
names_minus = ['nuR', 'eR', 'Lbar', 'Qbar', 'uR', 'dR']
comps += [(j, i, -1, n) for (i, j), n in zip(Nplus, names_minus)]
off, o = [], 0
for (i, j, s, n) in comps:
    off.append(o); o += dim[i]*dim[j]
NH = o
assert NH == 32
def blk(c): return slice(off[c], off[c] + dim[comps[c][0]]*dim[comps[c][1]])
gamma = np.zeros((NH, NH))
for c, (i, j, s, n) in enumerate(comps): gamma[blk(c), blk(c)] = s*np.eye(dim[i]*dim[j])
def left(cvec):
    M = np.zeros((NH, NH), complex)
    for c, (i, j, s, n) in enumerate(comps):
        M[blk(c), blk(c)] = np.kron(rep(i, cvec), np.eye(dim[j]))
    return M
# J : X in (i,j,s) |-> X^* in (j,i,-s) ;  J v = P conj(v)
P = np.zeros((NH, NH))
partner = {}
for c, (i, j, s, n) in enumerate(comps):
    for c2, (i2, j2, s2, n2) in enumerate(comps):
        if (i2, j2, s2) == (j, i, -s): partner[c] = c2
for c, (i, j, s, n) in enumerate(comps):
    c2 = partner[c]
    for r in range(dim[i]):
        for q in range(dim[j]):
            P[off[c2] + q*dim[i] + r, off[c] + r*dim[j] + q] = 1
def Jop(X): return P @ np.conj(X) @ P            # J X J^{-1}
def right(cvec): return Jop(left(star(cvec)))       # b° = J b* J^{-1}

# ---------- axioms ----------
res = {}
res['J2'] = np.abs(P @ P - np.eye(NH)).max()                       # J^2 = P conj(P) = P^2
res['Jgamma'] = np.abs(P @ gamma + gamma @ P).max()                  # J gamma = - gamma J
res['gamma2'] = np.abs(gamma @ gamma - np.eye(NH)).max()
res['P_orth'] = np.abs(P @ P.T - np.eye(NH)).max()                   # J antiunitary
L = [left(c) for c in AF_basis]
R = [right(c) for c in AF_basis]
res['left_is_rep_star'] = max(np.abs(left(star(c)) - left(c).conj().T).max() for c in AF_basis)
res['gamma_even'] = max(np.abs(gamma @ x - x @ gamma).max() for x in L)
res['order0'] = max(np.abs(a @ b - b @ a).max() for a in L for b in R)
res['u_unitary'] = max(np.abs(x @ x.conj().T - np.eye(NH)).max() for x in L)
# multiplicativity of left action on A_F (pi(u_g u_h) = pi(u_g) pi(u_h))
def gprod(c1, c2):
    v = np.zeros(24, complex)
    for g in G:
        for h in G: v[idx[mul(g, h)]] += c1[idx[g]]*c2[idx[h]]
    return v
res['left_mult'] = max(np.abs(left(gprod(a, b)) - left(a) @ left(b)).max()
                       for a in AF_basis[:6] for b in AF_basis)
res['right_antimult'] = max(np.abs(right(gprod(a, b)) - right(b) @ right(a)).max()
                            for a in AF_basis[:6] for b in AF_basis)
for k, v in res.items(): print(f"{k:18s} {v:.2e}")
# faithfulness / real dimension of pi(A_F)
Mreal = np.array([np.concatenate([x.real.ravel(), x.imag.ravel()]) for x in L])
print("dim_R pi(A_F) =", np.linalg.matrix_rank(Mreal, tol=1e-8))
# content : left multiplicities by parity
def content(sign, side):
    out = []
    for k in IRR:
        tot = 0
        for g in G:
            if side == 'L': Rg = left(np.eye(24)[idx[g]].astype(complex))
            else:           Rg = Jop(left(np.eye(24)[idx[inv(g)]].astype(complex)))
            proj = (np.eye(NH) + sign*gamma)/2
            tot += np.trace(proj @ Rg).real * chi[k][idx[g]]
        out.append(int(round(tot/24)))
    return out
print("gauche  + :", content(+1, 'L'), " - :", content(-1, 'L'))
print("droite  + :", content(+1, 'R'), " - :", content(-1, 'R'))
np.save(HERE + 'P.npy', P)
