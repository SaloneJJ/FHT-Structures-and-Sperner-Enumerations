"""Relecture Ec4 (V149 §3.2) : dimension reelle de {B hermitien, impair (B gamma_F = -gamma_F B), J_F B J_F^-1 = B,
[B, pi(a)] = 0 et [B, J_F pi(b) J_F^-1] = 0 pour a, b dans A_F}. Base reelle COMPLETE de A_F (24 elements, x2 avec J).
Methode : (1) reduction exacte par un element diagonal generique de l'algebre (B_ij = 0 si etiquettes differentes) ;
(2) contraintes lineaires reelles restantes, SVD."""
import numpy as np, itertools
from rel_F import *
basis = real_basis(); gens = basis + [Jconj(UF, b) for b in basis]
rng = np.random.default_rng(5)
diag_gens = [g for g in gens if np.allclose(g, np.diag(np.diag(g)))]
Dg = sum(rng.normal() * g for g in diag_gens)
lab = np.round(np.diag(Dg), 10)
pairs = [(i, j) for i in range(dF) for j in range(dF) if lab[i] == lab[j]]
print(f'{len(diag_gens)} generateurs diagonaux ; {len(set(lab))} etiquettes ; inconnues complexes apres reduction : {len(pairs)}')
nu = len(pairs)
def mat(v):   # v : vecteur reel de longueur 2 nu -> B
    B = np.zeros((dF, dF), complex)
    for t, (i, j) in enumerate(pairs): B[i, j] = v[2 * t] + 1j * v[2 * t + 1]
    return B
def nullspace(conds):
    Gm = np.zeros((2 * nu, 2 * nu))
    for f in conds:          # matrice de Gram sum_f A_f^T A_f (evite de stocker A)
        cols = []
        for t in range(2 * nu):
            e = np.zeros(2 * nu); e[t] = 1; c = f(mat(e)).ravel(); cols.append(np.concatenate([c.real, c.imag]))
        A = np.array(cols).T; Gm += A.T @ A
    w, v = np.linalg.eigh(Gm); return v[:, w < 1e-9 * w.max()].T
comm = [lambda B, X=X: B @ X - X @ B for X in gens]
herm = lambda B: B - B.conj().T; odd = lambda B: B @ GF + GF @ B; jinv = lambda B: Jconj(UF, B) - B
ns_c = nullspace(comm); print('commutant (dimension REELLE, = 2 x complexe) :', len(ns_c), ' attendu 2 x (10 x 9 + 36) = 252')
ns_all = nullspace(comm + [herm, odd, jinv]); print('Ec4 : hermitien + impair + J_F-invariant + commutant : dimension reelle', len(ns_all))
ns_even = nullspace(comm + [herm, jinv]); print('temoin sans la condition impair :', len(ns_even), '(auteurs : 66)')
ns_noJ = nullspace(comm + [herm, odd]); print('sans la condition J_F :', len(ns_noJ), '(attendu 18 : C complexe 3x3 quelconque)')
supp = np.zeros(dF, bool)
for v in ns_all:
    B = mat(v); supp |= (np.abs(B).sum(0) > 1e-9) | (np.abs(B).sum(1) > 1e-9)
lbl = [(('I', 'II')[i // dI], (i % dI) // (nB * nG), ((i % dI) // nG) % nB, i % nG) for i in np.where(supp)[0]]
print('support :', int(supp.sum()), 'etats :', lbl)
# controle : un B de l'espace a-t-il la forme Majorana [[0, C^+],[C, 0]] avec C symetrique ?
v = ns_all[0]; B = mat(v); I_ = [idx(0, 0, g, 0) for g in range(3)]; II_ = [idx(0, 0, g, 1) for g in range(3)]
C = B[np.ix_(II_, I_)]; print('C symetrique :', np.allclose(C, C.T))
# generateurs aleatoires (methode des auteurs) : combien en faut-il ?
for nrand in (1, 2, 4):
    gr = []
    for _ in range(nrand):
        a = random_alg(rng); ja = Jconj(UF, a); gr += [a + a.conj().T, 1j * (a - a.conj().T), ja + ja.conj().T, 1j * (ja - ja.conj().T)]
    print(f'  avec {nrand} element(s) aleatoire(s) (parties hermitienne ET antihermitienne) : dimension', len(nullspace([lambda B, X=X: B @ X - X @ B for X in gr] + [herm, odd, jinv])))
gr = []
for _ in range(3):
    a = random_alg(rng); ja = Jconj(UF, a); gr += [a + a.conj().T, ja + ja.conj().T]
print('  avec 3 elements aleatoires, parties hermitiennes seulement (a + a^+) :', len(nullspace([lambda B, X=X: B @ X - X @ B for X in gr] + [herm, odd, jinv])),
      '(les parties hermitiennes ne separent pas lambda de conj(lambda))')

print('\n--- test equitable des generateurs aleatoires (reduction par les sous-espaces propres d une combinaison des generateurs fournis, sans base diagonale) ---')
def dim_with(gr):
    X = sum(rng.normal() * g for g in gr); w, V = np.linalg.eigh(X)
    grp = []; cur = [0]
    for i in range(1, dF):
        if abs(w[i] - w[cur[0]]) < 1e-8: cur.append(i)
        else: grp.append(cur); cur = [i]
    grp.append(cur)
    unk = [(a, b) for g in grp for a in g for b in g]; n_ = len(unk)
    def mk(v):
        B = np.zeros((dF, dF), complex)
        for t, (a, b) in enumerate(unk): B[a, b] = v[2 * t] + 1j * v[2 * t + 1]
        return V @ B @ V.conj().T
    Gm = np.zeros((2 * n_, 2 * n_))
    for f in [lambda B, Y=Y: B @ Y - Y @ B for Y in gr] + [herm, odd, jinv]:
        cols = []
        for t in range(2 * n_):
            e = np.zeros(2 * n_); e[t] = 1; c = f(mk(e)).ravel(); cols.append(np.concatenate([c.real, c.imag]))
        A = np.array(cols).T; Gm += A.T @ A
    ww = np.linalg.eigvalsh(Gm); return int(np.sum(ww < 1e-9 * ww.max())), n_
for nrand, both in ((3, False), (6, False), (2, True), (4, True)):
    gr = []
    for _ in range(nrand):
        a = random_alg(rng); ja = Jconj(UF, a)
        gr += [a + a.conj().T, ja + ja.conj().T] + ([1j * (a - a.conj().T), 1j * (ja - ja.conj().T)] if both else [])
    print(f'  {nrand} elements aleatoires, {"parties hermitienne et antihermitienne" if both else "parties hermitiennes seules"} : dimension {dim_with(gr)[0]}')
