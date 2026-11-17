#!/usr/bin/env python3
"""D20 (c) -- espaces REELS des D_F admissibles et partie finie de l'action de masse.
Criteres : FHT_28/04_paysage/CRITERES_D20_fonctionnelle_a_posteriori.md (ecrits avant calcul).
H_F = formes harmoniques reelles (noyau de d + d* par degre, base orthonormee) : produit scalaire canonique, Aut agit orthogonalement.
J' : intertwiner de D17 (A s(g) rho_imp(sigma g sigma) = rho_pair(g) A), rendu orthogonal par decomposition polaire ; J' = [[0, U], [U^T, 0]] o conjugaison.
D = [[0, Q^+], [Q, 0]] (hermitien, impair), Q = X + iY (impair <- pair) ; contraintes : J'D = DJ' ; ordre un pour les generateurs purement gauches ;
photon (option) : [D, rho(g)] = 0 pour le stabilisateur du sommet p.
Action de masse (partie finie, T1-T3 des criteres) : kappa = 48 Tr (Q^+Q)^2 / (Tr Q^+Q)^2 (= 1 ssi rang 48 et spectre degenere), minimisee sur l'espace admissible ;
et la valeur exacte chaleur sum_i exp(-t (eps m)^2 s_i^2) (s normalise : sum s_i^2 = 48), eps m = 0,1, t = 10, 20, 40.
Usage : python3 D20_DF_reel.py K4 K4:photon X_comp X_comp:photon X_chi K26 ... | tee D20_c.log"""
import sys, os, json, itertools, time, random
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.optimize import minimize
H_ = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(H_, 'D19_photon.py')).read()
exec(src.split("pv = {")[0].replace("def analyse(name, X, pvert=None):", "def _analyse_D19(name, X, pvert=None):"))
rng = np.random.default_rng(20261115)

def harmonic(C):
    """bases orthonormees des formes harmoniques reelles par degre (d >= 1 ; homologie reduite = non reduite pour d >= 1)"""
    B = {}
    for d in range(1, C.dmax + 1):
        r, c, v = [], [], []
        for j, s in enumerate(C.S[d]):
            for t in range(len(s)): r.append(C.idx[d - 1][s[:t] + s[t + 1:]]); c.append(j); v.append((-1.0) ** t)
        B[d] = sp.csr_matrix((v, (r, c)), shape=(len(C.S[d - 1]), len(C.S[d])))
    Z = {}
    for d in range(1, C.dmax + 1):
        L = B[d].T @ B[d] + (B[d + 1] @ B[d + 1].T if d + 1 in B else 0)
        L = sp.csc_matrix(L, dtype=float); n = L.shape[0]
        if n <= 3000:
            w, U = np.linalg.eigh(L.toarray())
        else:
            w, U = sla.eigsh(L, k=min(n - 2, 80), sigma=-1e-2, which='LM')
        z = w < 1e-8
        if z.sum():
            assert z.sum() < len(w) - 1 and w[~z].min() > 1e-4, ('noyau mal separe', d, z.sum())
            Q_, _ = np.linalg.qr(U[:, z]); Z[d] = Q_
    return Z

def perm_matrix(C, d, g):
    r, c, v = [], [], []
    for j, s in enumerate(C.S[d]):
        t = [g[i] for i in s]; inv = sum(1 for a in range(len(t)) for b in range(a + 1, len(t)) if t[a] > t[b])
        r.append(C.idx[d][tuple(sorted(t))]); c.append(j); v.append(-1.0 if inv % 2 else 1.0)
    return sp.csr_matrix((v, (r, c)), shape=(len(C.S[d]),) * 2)

def rho(C, Z, g, degs):
    blocks = [Z[d].T @ (perm_matrix(C, d, g) @ Z[d]) for d in degs]
    n = sum(b.shape[0] for b in blocks); M = np.zeros((n, n)); o = 0
    for b in blocks: k = b.shape[0]; M[o:o + k, o:o + k] = b; o += k
    return M

def reduce_space(B, constraints, chunk=256, tol=1e-9):
    """B : base orthonormee (lignes) de l'espace courant ; chaque contrainte f(E) -> (k, m) reelle lineaire ; on garde les combinaisons annulees"""
    for f in constraints:
        if B.shape[0] == 0: break
        R = np.concatenate([f(B[a:a + chunk]) for a in range(0, B.shape[0], chunk)], axis=0)
        G = R @ R.T; w, V = np.linalg.eigh(G); sc = max(w.max(), 1.0)
        B = V[:, w < tol * sc].T @ B
        B, _ = np.linalg.qr(B.T); B = B.T
    return B

def analyse(name, X, pvert=None):
    t0 = time.time(); C = AF9.Cx(name, X); nv = len(C.V); V = C.V
    nonL = [C.ix[v] for v in V if not v.startswith('L_')]; Lv = [C.ix[v] for v in V if v.startswith('L_')]
    if pvert: Lv = [C.ix[pvert]] + [x for x in Lv if x != C.ix[pvert]]
    gens, order, levels = aut_sgs_base(C, nonL + Lv)
    left = [g for g in gens if all(g[u] == u for u in nonL)]
    D4 = [g for g in left if g[C.ix[pvert]] == C.ix[pvert]] if pvert else []
    sig = [C.ix[AF9.sw(v)] for v in V]; gidx = [C.ix[v] for v in ('g0', 'g1', 'g2', 'g3')]
    def sgn_g(g):
        pos = {p: k for k, p in enumerate(gidx)}; perm = [pos[g[p]] for p in gidx]
        return -1 if sum(1 for a in range(4) for b in range(a + 1, 4) if perm[a] > perm[b]) % 2 else 1
    Z = harmonic(C); betti = {d: z.shape[1] for d, z in Z.items()}
    ev = [d for d in sorted(Z) if d % 2 == 0]; od = [d for d in sorted(Z) if d % 2 == 1]
    nE = sum(betti[d] for d in ev); nO = sum(betti[d] for d in od)
    print(f'== {name} : |Aut| = {order} ; {len(gens)} generateurs dont {len(left)} gauches ; photon {len(D4)} ; Betti reels {betti} ({time.time()-t0:.0f} s)', flush=True)
    # J'
    pairs = []
    for g in gens:
        h = tuple(sig[g[sig[i]]] for i in range(nv))
        pairs.append((rho(C, Z, g, ev), sgn_g(g) * rho(C, Z, h, od)))
    orth = max(np.abs(M @ M.T - np.eye(nE)).max() for M, _ in pairs)
    def mapA(E):
        out = []
        for M, N in pairs:
            A = E.reshape(-1, nE, nO); out.append((A @ N - np.einsum('ij,kjl->kil', M, A)).reshape(E.shape[0], -1))
        return np.concatenate(out, axis=1)
    BA = reduce_space(np.eye(nE * nO), [lambda E, M=M, N=N: (E.reshape(-1, nE, nO) @ N - np.einsum('ij,kjl->kil', M, E.reshape(-1, nE, nO))).reshape(E.shape[0], -1) for M, N in pairs])
    A = np.tensordot(rng.standard_normal(BA.shape[0]), BA, axes=1).reshape(nE, nO)
    u, s, vt = np.linalg.svd(A); U = u @ vt                          # polaire : intertwiner orthogonal
    resJ = max(np.abs(U @ N - M @ U).max() for M, N in pairs)
    print(f'   rho orthogonale (ecart {orth:.1e}) ; dim des intertwiners J\' = {BA.shape[0]} ; U orthogonal, residu d\'entrelacement {resJ:.1e}', flush=True)
    Me = {g: rho(C, Z, g, ev) for g in left}; Mo = {g: rho(C, Z, g, od) for g in left}
    Jl = np.block([[np.zeros((nE, nE)), U], [U.T, np.zeros((nO, nO))]])
    R = {g: np.block([[Me[g], np.zeros((nE, nO))], [np.zeros((nO, nE)), Mo[g]]]) for g in left}
    Ro = {g: Jl @ R[g] @ Jl.T for g in left}                          # J rho(g) J^-1 (rho reelle)
    oz = max(np.abs(R[g] @ Ro[h] - Ro[h] @ R[g]).max() for g in left for h in left)
    n = nE + nO; k = nO * nE
    def Dbatch(E):                                                    # parametres (X, Y) -> D complexes
        Xs = E[:, :k].reshape(-1, nO, nE); Ys = E[:, k:].reshape(-1, nO, nE); Q = Xs + 1j * Ys
        D = np.zeros((E.shape[0], n, n), complex); D[:, nE:, :nE] = Q; D[:, :nE, nE:] = np.conj(np.transpose(Q, (0, 2, 1))); return D
    def cplx(o, E): o = o.reshape(E.shape[0], -1); return np.concatenate([o.real, o.imag], axis=1)
    cons = [lambda E: cplx(np.einsum('ij,kjl->kil', Jl, np.conj(Dbatch(E))) - np.einsum('kij,jl->kil', Dbatch(E), Jl), E)]
    for g in D4: cons.append(lambda E, g=g: cplx(np.einsum('kij,jl->kil', Dbatch(E), R[g]) - np.einsum('ij,kjl->kil', R[g], Dbatch(E)), E))
    for g in left:
        for h in left:
            def f(E, g=g, h=h):
                D = Dbatch(E); Cg = np.einsum('kij,jl->kil', D, R[g]) - np.einsum('ij,kjl->kil', R[g], D)
                return cplx(np.einsum('kij,jl->kil', Cg, Ro[h]) - np.einsum('ij,kjl->kil', Ro[h], Cg), E)
            cons.append(f)
    BD = reduce_space(np.eye(2 * k), cons, chunk=128)
    dim = BD.shape[0]
    print(f'   ordre zero (ecart {oz:.1e}) ; dimension REELLE des D_F admissibles (hermitiens) = {dim}  ({time.time()-t0:.0f} s)', flush=True)
    res = dict(name=name, betti=betti, dim_reelle=dim)
    if dim == 0: return res
    Qb = BD[:, :k].reshape(dim, nO, nE) + 1j * BD[:, k:].reshape(dim, nO, nE)
    ranks = [np.linalg.matrix_rank(np.tensordot(rng.standard_normal(dim), Qb, axes=1), tol=1e-8) for _ in range(3)]
    def kappa_and_grad(c):
        Q = np.tensordot(c, Qb, axes=1); M = Q.conj().T @ Q; s2 = np.trace(M).real; s4 = np.sum(np.abs(M) ** 2)
        f = nE * s4 / s2 ** 2
        # derivee : d s4 = 4 Re Tr(M dQ^+ Q) ; d s2 = 2 Re Tr(Q^+ dQ)
        G4 = 4 * np.real(np.einsum('kij,ij->k', Qb.conj(), Q @ M)); G2 = 2 * np.real(np.einsum('kij,ij->k', Qb.conj(), Q))
        return f, nE * (G4 / s2 ** 2 - 2 * s4 * G2 / s2 ** 3)
    best = None
    for trial in range(20):
        r = minimize(kappa_and_grad, rng.standard_normal(dim), jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-12))
        if best is None or r.fun < best.fun: best = r
    Q = np.tensordot(best.x, Qb, axes=1); s = np.linalg.svd(Q, compute_uv=False); s = s * np.sqrt(nE / np.sum(s ** 2))
    heat = {t: float(np.sum(np.exp(-t * 0.01 * s ** 2)) + (nE - np.sum(s > 1e-9)) * 0) for t in (10, 20, 40)}
    res.update(rang_generique=ranks, kappa_min=float(best.fun), valeurs_singulieres_min=np.round(np.sort(s)[::-1], 4).tolist(), chaleur_somme_paires=heat)
    print(f'   rang generique de Q = {ranks} (max {nE}) ; kappa_min = {best.fun:.6f} (1 = rang plein degenere) ; valeurs singulieres (normalisees, somme des carres {nE}) :',
          np.round(np.sort(s)[::-1], 3).tolist(), flush=True)
    print(f'   action de masse (chaleur, paires) sum exp(-t (eps m)^2 s^2) : {heat}  ({time.time()-t0:.0f} s)', flush=True)
    return res

pv = {'K4': 'L_p0', 'K26': 'L_p0', 'X_comp': 'L_l2'}
out = []
for arg in sys.argv[1:]:
    name, _, opt = arg.partition(':')
    out.append(analyse(name + (' (photon)' if opt else ''), cx[name](), pv[name] if opt else None))
json.dump(out, open(os.path.join(H_, 'D20_c_' + '_'.join(a.replace(':', '-').replace("'", 'p') for a in sys.argv[1:]) + '.json'), 'w'), indent=1, default=float)
