#!/usr/bin/env python3
"""AUDIT FABLE -- A2, section 4 (registre T2, T5) : controles combinatoires independants.
 (1) Proposition « Unique decomposition » : la matrice des multiplicites T du bimodule d'une generation s'ecrit a b^t + c d^t (a, b, c, d dans N^5)
     d'une seule facon (a l'echange des deux termes pres) : recherche exhaustive.
 (2) Theoreme « The one-generation triple » : dimension reelle 32 de l'espace des D admissibles (comptage des blocs par la regle de l'enonce).
 (3) Theoreme « Photon » : restrictions a D_4 = Stab(appariement) des representations de S_4 ; <chi_3, 1> = <chi_3, sgn> = 0 ; Res 2 = 1 + sgn ;
     comptage 84 = 4 x 18 + 12.
 (4) Proposition « Intersection form » : matrice 3 [[0,-2,-2],[2,0,2],[2,-2,0]], antisymetrique, de rang 2.
 (5) Lemme « Opposite faces » : H_1(K_4) = 3' = 3 (x) sgn (caracteres).
Usage : python3 AF8_bimodule_controles.py      (numpy)
"""
import itertools, numpy as np
irr = ['1', 's', '2', '3', "3'"]; dim = {'1': 1, 's': 1, '2': 2, '3': 3, "3'": 3}
Np = [('1', '1'), ('1', 's'), ('2', '1'), ('2', '3'), ('3', '1'), ('3', 's')]; Nm = [(j, i) for i, j in Np]
T = np.zeros((5, 5), int)
for i, j in Np + Nm: T[irr.index(i), irr.index(j)] += 1
ok = True
print('(1) T =\n', T, '; support :', int((T > 0).sum()), 'entrees ; dimension', sum(dim[i] * dim[j] for i, j in Np + Nm))
sols = set()
vecs = list(itertools.product(range(3), repeat=5))
for a in vecs:
    if not any(a): continue
    for b in vecs:
        if not any(b): continue
        R = T - np.outer(a, b)
        if (R < 0).any() or not R.any(): continue
        # R = c d^t ? rang 1 a coefficients entiers >= 0 : c = une colonne non nulle divisee par un diviseur commun, d determine
        cols = [k for k in range(5) if R[:, k].any()]
        c0 = R[:, cols[0]]
        for div in (1, 2):
            if (c0 % div).any(): continue
            c = c0 // div; i0 = int(np.nonzero(c)[0][0])
            if (R[i0, :] % c[i0]).any(): continue
            dd = R[i0, :] // c[i0]
            if (np.outer(c, dd) == R).all():
                sols.add(frozenset([(tuple(a), tuple(b)), (tuple(int(x) for x in c), tuple(int(x) for x in dd))]))
print('    decompositions T = a b^t + c d^t :', len(sols))
for s in sols: print('    ', sorted(s))
exp = frozenset([((1, 1, 1, 0, 0), (1, 0, 0, 1, 0)), ((1, 0, 0, 1, 0), (1, 1, 1, 0, 0))])
ok &= sols == {exp}; print('    unique et egale a A B^t + B A^t, A = 1 + s + 2, B = 1 + 3 :', sols == {exp})

print('(2) Espace des D admissibles sur une generation')
# blocs entre composantes de chiralites opposees (N+ <-> N-) avec i = k ou j = l ; on compte les paires non ordonnees {X, Y}, puis J identifie un bloc et son J-partenaire
tot = 0; detail = []; majorana = []
for (i, j) in Np:
    for (k, l) in Nm:
        if i != k and j != l: continue                      # bloc nul (regle de l'enonce)
        if i == k and j == l:                               # composante commune a N+ et N- : (1,1) = nu_R-bar -> nu_R, bloc de Majorana
            majorana.append(((i, j), dim[i] * dim[j])); continue
        d_c = dim[j] * dim[l] if i == k else dim[i] * dim[k]   # 1 (x) S, S dans Hom(V_j, V_l), ou S (x) 1, S dans Hom(V_i, V_k)
        detail.append(((i, j), (k, l), d_c))
print('    blocs N+ -> N- autorises (dimension complexe) :', [(x, y, d_) for x, y, d_ in detail], '; Majorana :', majorana)
# D autoadjoint : la moitie N+ -> N- determine D. J D = D J identifie le bloc x -> y au bloc J y -> J x (J echange (i, j) et (j, i)) : orbites de taille 2
pairs = {}
for (x, y, d_c) in detail:
    jx, jy = (y[1], y[0]), (x[1], x[0])
    pairs[frozenset([(x, y), (jx, jy)])] = d_c
assert all(len(k) == 2 for k in pairs)
dimC = sum(pairs.values()); dimR = 2 * dimC + 2 * sum(d_ for _, d_ in majorana)      # Majorana (1x1, fixe par l'involution) : un nombre complexe
print(f'    {len(pairs)} orbites de blocs (dimension complexe {dimC}) + Majorana : dimension reelle = 2 x {dimC} + 2 = {dimR} (article : 2(2+2+2+2+1+3+3) + 2 = 32)'); ok &= dimR == 32

print('(3) Restrictions a D_4 = Stab(appariement {01|23})')
S4 = list(itertools.permutations(range(4)))
def sgn(p): return (-1) ** sum(1 for a in range(4) for b in range(a + 1, 4) if p[a] > p[b])
def chi3(p): return sum(1 for a in range(4) if p[a] == a) - 1
pair = [frozenset({frozenset({0, 1}), frozenset({2, 3})}), frozenset({frozenset({0, 2}), frozenset({1, 3})}), frozenset({frozenset({0, 3}), frozenset({1, 2})})]
def act(p, pr): return frozenset(frozenset(p[x] for x in blk) for blk in pr)
def chi2(p): return sum(1 for pr in pair if act(p, pr) == pr) - 1
D4 = [p for p in S4 if act(p, pair[0]) == pair[0]]
ip = lambda f, g: sum(f(p) * g(p) for p in D4) / len(D4)
one = lambda p: 1
print(f'    |D_4| = {len(D4)} ; <chi_3, 1> = {ip(chi3, one)}, <chi_3, sgn> = {ip(chi3, sgn)}, <chi_2, 1> = {ip(chi2, one)}, <chi_2, sgn> = {ip(chi2, sgn)}, <1, sgn> = {ip(one, sgn)}')
ok &= len(D4) == 8 and ip(chi3, one) == 0 and ip(chi3, sgn) == 0 and ip(chi2, one) == 1 and ip(chi2, sgn) == 1 and ip(one, sgn) == 0
print('    parametres : 4 Yukawa dans M_3(C) (4 x 18) + Majorana symetrique 3x3 complexe (12) =', 4 * 18 + 12, '(article : 84)')
# sous-groupes : seul D_4 (3 conjugues) donne exactement cela ? on verifie au moins que A_4 et S_3 ne conviennent pas
A4 = [p for p in S4 if sgn(p) == 1]; S3 = [p for p in S4 if p[3] == 3]
for name, G in (('A_4', A4), ('S_3', S3), ('V_4', [p for p in A4 if chi3(p) in (3, -1)])):
    ipG = lambda f, g: sum(f(p) * g(p) for p in G) / len(G)
    print(f'    {name} : <chi_2,1> = {ipG(chi2, one)}, <chi_2,sgn> = {ipG(chi2, sgn)}, <chi_3,1> = {ipG(chi3, one)}, <chi_3,sgn> = {ipG(chi3, sgn)}, <1,sgn> = {ipG(one, sgn)}')

print('(4) Forme d\'intersection')
r = {'C': {'1': 1, 's': 1}, 'H': {'2': 2}, '3': {'3': 1, "3'": 1}}
M = np.array([[3 * (sum(r[p].get(i, 0) * r[q].get(j, 0) for i, j in Np) - sum(r[p].get(i, 0) * r[q].get(j, 0) for i, j in Nm)) for q in 'CH3'] for p in 'CH3'])
print('   ', M.tolist(), '; antisymetrique :', (M == -M.T).all(), '; rang', np.linalg.matrix_rank(M))
ok &= (M == 3 * np.array([[0, -2, -2], [2, 0, 2], [2, -2, 0]])).all() and np.linalg.matrix_rank(M) == 2

print('(5) H_1(K_4) comme S_4-module : caractere = chi(C_1) - chi(C_0) + 1 (K_4 connexe)')
edges = list(itertools.combinations(range(4), 2))
def chiC1(p):
    t = 0
    for (a, b) in edges:
        if {p[a], p[b]} == {a, b}: t += 1 if p[a] == a else -1
    return t
chiH1 = lambda p: chiC1(p) - (chi3(p) + 1) + 1
same = all(chiH1(p) == sgn(p) * chi3(p) for p in S4); print("    chi(H_1(K_4)) = sgn x chi_3 (= 3') :", same); ok &= same
print('VERDICT :', 'conforme' if ok else 'ECART')
