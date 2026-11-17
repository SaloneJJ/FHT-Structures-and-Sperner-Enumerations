import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__)) + '/'
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(HERE + 'triplet.py').read())
import numpy as np

def gel(g): return np.eye(24)[idx[g]].astype(complex)
def central(k):
    v = np.zeros(24, complex)
    for g in G: v[idx[g]] = dim[k]/24*chi[k][idx[inv(g)]]
    return v
Plus = [c for c in range(12) if comps[c][2] == +1]
Minus = [c for c in range(12) if comps[c][2] == -1]
ip = np.concatenate([np.arange(blk(c).start, blk(c).stop) for c in Plus])
im = np.concatenate([np.arange(blk(c).start, blk(c).stop) for c in Minus])

# unknowns : odd D, entries (p in H+, q in H-) and (q, p)
unk = [(int(a), int(b)) for a in im for b in ip] + [(int(b), int(a)) for a in im for b in ip]
nU = len(unk)
def Dfrom(z):
    D = np.zeros((NH, NH), complex)
    for k, (r, c) in enumerate(unk): D[r, c] = z[k]
    return D

def lin_constraints(fun):
    """rows of complex-linear map z -> fun(D(z)) (flattened)"""
    cols = []
    for k in range(nU):
        z = np.zeros(nU, complex); z[k] = 1
        cols.append(fun(Dfrom(z)).ravel())
    return np.array(cols).T

def nullspace(M, tol=1e-9):
    if M.shape[0] == 0: return np.eye(M.shape[1], dtype=M.dtype)
    u_, s, vh = np.linalg.svd(M, full_matrices=False)
    r = (s > tol*max(1, s[0])).sum()
    return vh[r:].conj().T, s

gens = [(1, 0, 2, 3), (1, 2, 3, 0)]   # (12), (1234) : generate S4
A_g = [left(gel(g)) for g in gens]
B_g = [right(gel(g)) for g in gens]
def order1(D):
    return np.concatenate([((D@a - a@D)@b - b@(D@a - a@D)).ravel() for a in A_g for b in B_g])

def real_dim(basis_complex, extra_real=None):
    """basis : list of complex matrices B_k ; D = sum z_k B_k. Impose D = D^*, J D J^-1 = D."""
    n = len(basis_complex)
    rows = []
    for k in range(n):
        for part in (1, 1j):
            B = part*basis_complex[k]
            h = B - B.conj().T
            j = B - Jop(B)
            rows.append(np.concatenate([h.ravel(), j.ravel()] + ([f(B).ravel() for f in extra_real] if extra_real else [])))
    M = np.array(rows).T
    Mr = np.vstack([M.real, M.imag])
    ns, s = nullspace(Mr)
    sols = []
    for v in ns.T:
        z = v[0::2] + 1j*v[1::2]
        sols.append(sum(z[k]*basis_complex[k] for k in range(n)))
    return sols

def complex_basis(M):
    ns, s = nullspace(M)
    return [Dfrom(v) for v in ns.T], s

print("unknowns (odd D, complex) :", nU)
# V0 / V1
allB = [Dfrom(np.eye(nU)[k]) for k in range(nU)]
def rd_simple(Bs, J=True):
    n = len(Bs); rows = []
    for k in range(n):
        for part in (1, 1j):
            B = part*Bs[k]
            r = [(B - B.conj().T).ravel()]
            if J: r.append((B - Jop(B)).ravel())
            rows.append(np.concatenate(r))
    M = np.array(rows).T; Mr = np.vstack([M.real, M.imag])
    return Mr.shape[1] - np.linalg.matrix_rank(Mr, tol=1e-9)
print("V0 impair autoadjoint            dim_R =", rd_simple(allB, J=False))
print("V1 + J'D = DJ'                   dim_R =", rd_simple(allB, J=True))

# V2 : order one
M1 = lin_constraints(order1)
B1, s1 = complex_basis(M1)
print("ordre un (impair, complexe)      dim_C =", len(B1), " gap sv:", s1[len(s1)-len(B1)-1] if len(B1) < len(s1) else None)
S2 = real_dim(B1)
print("V2 + autoadjoint + J'            dim_R =", len(S2))
# check full order one on all 24x24 group elements and A_F basis
Lall = [left(gel(g)) for g in G]; Rall = [right(gel(g)) for g in G]
err = 0
for D in S2:
    err = max(err, np.abs(D - D.conj().T).max(), np.abs(D - Jop(D)).max(), np.abs(D@gamma + gamma@D).max())
    for a in Lall:
        C = D@a - a@D
        for b in Rall: err = max(err, np.abs(C@b - b@C).max())
print("   verif complete ordre un / J / impair / autoadj :", f"{err:.1e}")

# block structure of V2
def block_rank(sols, c1, c2):
    vecs = [np.concatenate([D[blk(c2), blk(c1)].real.ravel(), D[blk(c2), blk(c1)].imag.ravel()]) for D in sols]
    return np.linalg.matrix_rank(np.array(vecs), tol=1e-8)
print("   blocs H+ -> H- non nuls (dim_R de la projection) :")
tot = 0
for c1 in Plus:
    for c2 in Minus:
        r = block_rank(S2, c1, c2)
        if r: print(f"     {comps[c1][3]:7s} -> {comps[c2][3]:5s} : {r}"); tot += r
print("   somme :", tot)

# V3 : photon / C_F condition, for the three transposition classes and both signs
e1, es, e2 = central('1'), central('s'), central('2')
results = {}
for tau in [(1, 0, 2, 3), (2, 1, 0, 3), (3, 1, 2, 0), (0, 1, 3, 2)]:
    for eps in (+1, -1):
        aF = 1j*(e1 - es) + eps*gprod(e2, u(tau))
        assert np.allclose(sigma(aF), aF)
        X = left(aF)
        S3 = real_dim(B1, extra_real=[lambda B, X=X: B@X - X@B])
        results[(tau, eps)] = len(S3)
print("V3 + [D, pi(a_CF)] = 0           dim_R par choix (tau, signe) :", results)
tau = (1, 0, 2, 3); aF = 1j*(e1 - es) + gprod(e2, u(tau)); X = left(aF)
S3 = real_dim(B1, extra_real=[lambda B: B@X - X@B])
for c1 in Plus:
    for c2 in Minus:
        r = block_rank(S3, c1, c2)
        if r: print(f"     {comps[c1][3]:7s} -> {comps[c2][3]:5s} : {r}")

# V4 : commute with left rho(S4)  (Theorem V19 §3)
rhoL = [left(gel(g)) for g in G]
S4L = real_dim(B1, extra_real=[lambda B, Y=Y: B@Y - Y@B for Y in [left(gel(g)) for g in gens]])
print("V4 + commute avec rho(S4) (gauche) dim_R =", len(S4L),
      " max|[D,pi(A_F)]| =", f"{max([0]+[np.abs(D@x - x@D).max() for D in S4L for x in L]):.1e}")
for c1 in Plus:
    for c2 in Minus:
        r = block_rank(S4L, c1, c2)
        if r: print(f"     {comps[c1][3]:7s} -> {comps[c2][3]:5s} : {r}")
# V5 : commute with the diagonal S4 = Ad(pi(u_g) J pi(u_g) J^-1)
diag = [left(u(g)) @ Jop(left(u(g))) for g in G]
for g, U in zip(G, diag):
    assert np.abs(U.imag).max() < 1e-12
S5 = real_dim(B1, extra_real=[lambda B, Y=Y: B@Y - Y@B for Y in [diag[idx[g]] for g in gens]])
print("V5 + commute avec S4 diagonal (jauge) dim_R =", len(S5))
for c1 in Plus:
    for c2 in Minus:
        r = block_rank(S5, c1, c2)
        if r: print(f"     {comps[c1][3]:7s} -> {comps[c2][3]:5s} : {r}")
S35 = real_dim(B1, extra_real=[lambda B: B@X - X@B] + [lambda B, Y=Y: B@Y - Y@B for Y in [diag[idx[g]] for g in gens]])
print("V3+V5 dim_R =", len(S35))

# commutants : moduli of the identification
def skew_commutant(ops, Jcond=False):
    n = NH; rows = []
    basis = []
    for r in range(n):
        for c in range(n):
            E = np.zeros((n, n), complex); E[r, c] = 1; basis += [E, 1j*E]
    cols = []
    for E in basis:
        parts = [(E + E.conj().T).ravel()] + [(E@o - o@E).ravel() for o in ops]
        if Jcond: parts.append((E - Jop(E)).ravel())
        cols.append(np.concatenate(parts))
    M = np.array(cols).T; Mr = np.vstack([M.real, M.imag])
    return Mr.shape[1] - np.linalg.matrix_rank(Mr, tol=1e-9)
c_rho_gamma = skew_commutant(rhoL[:4] + [left(gel(g)) for g in [(1,2,3,0),(1,0,2,3)]] + [gamma.astype(complex)])
c_full = skew_commutant([left(gel(g)) for g in gens] + [right(gel(g)) for g in gens] + [gamma.astype(complex)], Jcond=True)
c_noJ = skew_commutant([left(gel(g)) for g in gens] + [right(gel(g)) for g in gens] + [gamma.astype(complex)])
print("dim Lie commutant (rho, gamma)            :", c_rho_gamma)
print("dim Lie commutant (gauche, droite, gamma) :", c_noJ, " avec J :", c_full)
