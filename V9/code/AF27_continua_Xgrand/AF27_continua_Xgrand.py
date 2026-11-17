#!/usr/bin/env python3
"""AF27 -- continus portes par X_grand (principe et filtre : 08_descente_Xgrand/CONTINUA_XGRAND.md, ecrits avant calcul).
(1) 2O dans SU(2) ; spin 3/2 = Sym^3 : irreductible, quaternionique ; GL(2,3) dans U(2) par isoclinisme (2T u i.(2O \ 2T)), verifie contre GL(2,3) sur F3 ;
    indicateur tordu par det.  (2) filtre de protection des candidats.  (3) S^4 = HP^1 : sous-groupes H de 2O qui coupent 4_spin en droites quaternioniques.
Usage : python3 AF27_continua_Xgrand.py | tee AF27.log   (< 1 min, numpy seul)
"""
import itertools
import numpy as np
from collections import Counter

def q2m(q):
    a, b, c, d = q
    return np.array([[a + 1j * b, c + 1j * d], [-c + 1j * d, a - 1j * b]])
def group_2O():
    S = set()
    for i in range(4):
        for s in (1, -1): v = [0.] * 4; v[i] = s; S.add(tuple(v))
    for s in itertools.product((.5, -.5), repeat=4): S.add(s)
    r = 1 / np.sqrt(2)
    for i, j in itertools.combinations(range(4), 2):
        for a in (r, -r):
            for b in (r, -r): v = [0.] * 4; v[i] = a; v[j] = b; S.add(tuple(np.round(v, 12)))
    return [q2m(np.array(x)) for x in S]
O2 = group_2O(); assert len(O2) == 48

def sym_power(M, n):
    """action de M (2x2) sur les polynomes homogenes de degre n en (x, y) : base x^(n-k) y^k"""
    from math import comb
    out = np.zeros((n + 1, n + 1), complex)
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    # (x, y) -> (a x + c y, b x + d y) : substitution
    for k in range(n + 1):
        # (a x + c y)^(n-k) (b x + d y)^k
        p = np.zeros(n + 1, complex)
        for i in range(n - k + 1):
            for j in range(k + 1):
                p[i + j] += comb(n - k, i) * comb(k, j) * (c ** i) * (a ** (n - k - i)) * (d ** j) * (b ** (k - j))
        out[:, k] = p
    # normaliser en base orthonormee (poids binomiaux) pour avoir une representation unitaire
    from math import sqrt
    W = np.diag([sqrt(comb(n, k)) for k in range(n + 1)])
    return np.linalg.inv(W) @ out @ W

def key(M): return tuple(np.round(M, 6).ravel())
def order(M):
    P = M.copy(); k = 1
    while not np.allclose(P, np.eye(len(M))): P = P @ M; k += 1
    return k
def char_norm(rep): return np.mean([abs(np.trace(R)) ** 2 for R in rep])
def fs(rep, eta=None):
    eta = eta if eta is not None else [1] * len(rep)
    return np.mean([e * np.trace(R @ R) for e, R in zip(eta, rep)]).real

# ---------------------------------------------------------------- (1)
print('(1) 2O : ordres', dict(sorted(Counter(order(g) for g in O2).items())))
# 2T = elements a coordonnees quaternioniques dans {0, +-1} ou {+-1/2}
def is2T(M):
    q = np.array([M[0, 0].real, M[0, 0].imag, M[0, 1].real, M[0, 1].imag])
    return np.allclose(np.sort(np.abs(q)), [0, 0, 0, 1]) or np.allclose(np.abs(q), .5)
T2 = [g for g in O2 if is2T(g)]; assert len(T2) == 24
GL = [g if is2T(g) else 1j * g for g in O2]
keys = {key(g) for g in GL}; closed = all(key(a @ b) in keys for a in GL for b in GL)
def gl23_orders():
    els = [np.array(m).reshape(2, 2) for m in itertools.product(range(3), repeat=4)]
    els = [m for m in els if round(np.linalg.det(m)) % 3 != 0]
    def ordF(m):
        P = m.copy(); k = 1
        while not np.array_equal(P % 3, np.eye(2, dtype=int)): P = (P @ m) % 3; k += 1
        return k
    return dict(sorted(Counter(ordF(m) for m in els).items()))
print('    GL(2,3) dans U(2) par isoclinisme : ferme', closed, '; ordres', dict(sorted(Counter(order(g) for g in GL).items())),
      '; GL(2,3) sur F3 :', gl23_orders())
detchar = [1 if is2T(g) else -1 for g in O2]          # caractere det de GL(2,3) (= signe du revetement S4 -> +-1, porte par le coset)
for name, G, eta in (('2O', O2, None), ('GL(2,3)', GL, detchar)):
    rho = [sym_power(g, 3) for g in G]
    unit = all(np.allclose(R @ R.conj().T, np.eye(4)) for R in rho)
    print(f'    {name} : Sym^3 (spin 3/2) unitaire {unit}, <chi,chi> = {char_norm(rho):.3f}, indicateur FS = {fs(rho):+.3f}',
          f', indicateur tordu par det = {fs(rho, eta):+.3f}' if eta else '')

# ---------------------------------------------------------------- (2) filtre de protection
def ncomp(rep):   # <chi, chi> : nombre de composantes (sans multiplicite) pour une rep. complexe
    return round(char_norm(rep), 6)
# S^4 = HP^1 : harmonique 1 = 5 (spin 2) ; sous O (image de 2O) : Sym^4
print('(2) S^4 (spin 2 = Sym^4) sous 2O : <chi,chi> =', ncomp([sym_power(g, 4) for g in O2]), '(1 = irreductible)')
print('    S^7 spinoriel (8 reel = 4_spin quaternionique) : irreductible sur R car FS = -1 et <chi,chi> = 1 (voir (1))')
# S5 et son 4 standard
def S5_std():
    P = []
    B = np.linalg.qr(np.vstack([np.eye(5)[:4] - np.eye(5)[1:5]]).T)[0]       # base orthonormee de la somme nulle (5 x 4)
    for p in itertools.permutations(range(5)):
        M = np.eye(5)[list(p)]; P.append(B.T @ M @ B)
    return P
S5 = S5_std(); assert all(np.allclose(R @ R.T, np.eye(4)) for R in S5)
print('    S5 : 4 standard <chi,chi> =', ncomp(S5))
def on_sym2_traceless(R):
    n = len(R); basis = []
    for i in range(n):
        for j in range(i, n):
            E = np.zeros((n, n)); E[i, j] = E[j, i] = 1; basis.append(E)
    Bm = np.array([b.ravel() for b in basis]).T
    T = np.array([(R @ b @ R.T).ravel() for b in basis]).T
    A = np.linalg.lstsq(Bm, T, rcond=None)[0]
    return A - np.eye(len(basis)) * 0
def on_antisym(R):
    n = len(R); basis = []
    for i in range(n):
        for j in range(i + 1, n):
            E = np.zeros((n, n)); E[i, j] = 1; E[j, i] = -1; basis.append(E)
    Bm = np.array([b.ravel() for b in basis]).T
    return np.linalg.lstsq(Bm, np.array([(R @ b @ R.T).ravel() for b in basis]).T, rcond=None)[0]
sym = [on_sym2_traceless(R) for R in S5]; anti = [on_antisym(R) for R in S5]
print('    S3 = bord du 4-simplexe : harmonique 2 (Sym^2 de 4, avec la trace) <chi,chi> =', ncomp(sym), '(1 + 4 + 5 attendu : 3)')
print('    CP^3 : su(4) = Lambda^2 (6) + i Sym^2_0 (9) sous S5 : <chi,chi> Lambda^2 =', ncomp(anti), '; Sym^2 avec trace =', ncomp(sym), '=> 15 = 6 + 4 + 5')
print('    Sp(2) et U(4) (actions L x R) : 4 (x) 4* irreductible pour G x G des que 4 l est pour G (produit tensoriel exterieur)')

# ---------------------------------------------------------------- (3) brisure M2(H) -> H + H : sous-groupes H de 2O
def closure(gens):
    els = {key(np.eye(2)): np.eye(2, dtype=complex)}; frontier = list(els.values())
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = a @ g; k = key(b)
                if k not in els: els[k] = b; new.append(b)
        frontier = new
    return els
subs = {}
for a, b in itertools.combinations_with_replacement(range(48), 2):
    H = closure([O2[a], O2[b]]); subs[frozenset(H)] = list(H.values())
print(f'(3) sous-groupes de 2O : {len(subs)} (engendres par deux elements)')
rows = Counter()
for Hk, H in subs.items():
    rho = [sym_power(h, 3) for h in H]; n = ncomp(rho)
    # composantes : on diagonalise le commutant par projection sur les classes ; ici on note <chi,chi> et l'indicateur FS de la restriction
    ind = fs(rho)
    sig = (len(H), dict(sorted(Counter(order(h) for h in H).items())).__repr__(), n, round(ind, 3))
    rows[sig] += 1
print('    (ordre, ordres des elements, <chi,chi> de 4_spin|H, FS de 4_spin|H) : nombre de sous-groupes conjugues')
for sig, c in sorted(rows.items()): print('     ', sig, ':', c)

# ---------------------------------------------------------------- (3b) droites quaternioniques invariantes (points fixes dans S^4 = HP^1)
rho48 = [sym_power(g, 3) for g in O2]
assert all(np.allclose(R @ R.conj().T, np.eye(4)) for R in rho48), 'Sym^3 non unitaire'
# structure quaternionique : J = C o conj, C rho(g)^* = rho(g) C, C C^* = -1
A = np.vstack([np.kron(np.eye(4), R.conj().T) - np.kron(R, np.eye(4)) for R in rho48])     # vec(C R^* - R C) = 0 (lignes)
_, s, Vh = np.linalg.svd(A); C = Vh[-1].conj().reshape(4, 4); C = C / np.sqrt(abs(np.trace(C @ C.conj().T)) / 4)
assert np.allclose(C @ C.conj(), -np.eye(4)), 'pas quaternionique'
def invariant_H_lines(H):
    rho = [sym_power(h, 3) for h in H]
    Am = np.vstack([np.kron(R, np.eye(4)) - np.kron(np.eye(4), R.T) for R in rho])          # commutant : M R = R M
    _, s, Vh = np.linalg.svd(Am); B = Vh[np.sum(s > 1e-9):].conj(); comm = [b.reshape(4, 4) for b in B]
    rng = np.random.default_rng(1); M = sum(rng.normal() * c for c in comm); M = M + M.conj().T   # element hermitien generique du commutant
    w, U = np.linalg.eigh(M); groups = []
    for x in np.unique(np.round(w, 6)):
        groups.append(U[:, np.isclose(w, x, atol=1e-5)])
    # sous-espaces propres ; J-stables ? (J W = C conj(W))
    out = []
    for G_ in groups:
        JW = C @ G_.conj(); P = G_ @ G_.conj().T
        out.append((G_.shape[1], bool(np.allclose(P @ JW, JW))))
    return len(comm), out
print('(3b) sous-groupes de brisure : dim du commutant ; sous-espaces propres (dim, stable par J)')
seen = set()
for Hk, H in sorted(subs.items(), key=lambda kv: len(kv[1])):
    sig = (len(H), tuple(sorted(Counter(order(h) for h in H).items())))
    if sig in seen or len(H) not in (8, 12, 16, 24): continue
    seen.add(sig); print('    ordre', len(H), dict(sig[1]), ':', invariant_H_lines(H))
