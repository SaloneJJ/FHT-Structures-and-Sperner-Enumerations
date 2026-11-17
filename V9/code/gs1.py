#!/usr/bin/env python3
"""GS1 -- grande symetrie M4(H)+M8(C) = (M2(H)+M4(C)) (x) M2(R) par une piece de doublement D. Criteres : CRITERES_GS1.md.
Usage : python3 gs1.py | tee GS1.log"""
import sys, os, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..', 'FHT_28', '06_code', 'AF20_algebres_tordues')))
from twisted_alg import closure, decompose, sign
from fht_tools import aut_order, betti_z2

def tree_aut(n, E):
    Es = {frozenset(e) for e in E}; c = 0
    for p in itertools.permutations(range(n)):
        if all(frozenset((p[a], p[b])) in Es for a, b in E): c += 1
    return c
T8 = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)]      # chemin a0..a6 + pendant en a2
T7 = [(0,1),(1,2),(2,3),(3,4),(4,5),(2,6)]            # chemin a0..a5 + pendant en a2
print('G-a  |Aut(T8)| =', tree_aut(8, T8), ' |Aut(T7)| =', tree_aut(7, T7))
def D_of(T, n, k=3, mark=0):
    F = [[f'c{i}_{a}', f'c{i}_{b}'] for i in range(k) for a, b in T]
    for i in range(mark): F.append([f'c{i}_{0}', f'm{i}'])            # marquage : un pendant de plus sur la copie i
    V = sorted({x for f in F for x in f}); return V, F
for nm, T, n in (('D8', T8, 8), ('D7', T7, 7)):
    V, F = D_of(T, n); print(f'     {nm} : {len(V)} sommets, Betti Z/2 {betti_z2(F)}, |Aut| = {aut_order(V, F)}')
    # signe de la permutation des sommets qui echange les copies 0 et 1
    perm = list(range(3 * n)); 
    for a in range(n): perm[a], perm[n + a] = n + a, a
    print(f'     {nm} : signe (permutation des sommets) d une transposition des copies = {sign(tuple(perm)):+d}')
for mk in (1, 2):
    V, F = D_of(T8, 8, 3, mk); print(f'G-e  D8 avec {mk} copie(s) marquee(s) : Betti {betti_z2(F)}, |Aut| = {aut_order(V, F)}')

# ---- groupes en permutations
F3 = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
def glperm(M):
    return tuple(F3.index(((M[0][0]*x + M[0][1]*y) % 3, (M[1][0]*x + M[1][1]*y) % 3)) for x, y in F3)
mats = [((a, b), (c, d)) for a in range(3) for b in range(3) for c in range(3) for d in range(3) if (a*d - b*c) % 3]
detof = {glperm(M): (1 if (M[0][0]*M[1][1] - M[0][1]*M[1][0]) % 3 == 1 else -1) for M in mats}
GL = sorted(detof); assert len(GL) == 48
S5 = closure([(1,0,2,3,4), (1,2,3,4,0)]); S3 = closure([(1,0,2), (1,2,0)]); S2 = closure([(1,0)])
def prod(G, H, eG, eH):
    n, m = len(G[0]), len(H[0]); P = []; ev = {}
    for g in G:
        for h in H:
            p = tuple(list(g) + [n + x for x in h]); P.append(p); ev[p] = eG(g) * eH(h)
    return sorted(P), (lambda p: ev[p])
one = lambda p: 1; sg = lambda p: sign(p); dt = lambda p: detof[p]
print('\nControles (A1) :')
decompose(GL, dt, 'GL(2,3), det'); decompose(S5, sg, 'S5, sgn'); decompose(S3, one, 'S3, trivial'); decompose(S3, sg, 'S3, sgn')
print('\nG-c  piece de doublement paire (eta trivial sur S3) :')
G, e = prod(GL, S3, dt, one); fs = decompose(G, e, 'GL(2,3) x S3, det x 1')
print('     -> facteur sur 4spin x 2 :', [f for f in fs if f in ('M4(H)', 'M8(R)', 'M8(C)')])
G, e = prod(S5, S3, sg, one); fc = decompose(G, e, 'S5 x S3, sgn x 1')
print('     -> facteur sur 4 x 2 :', [f for f in fc if f in ('M4(H)', 'M8(R)', 'M8(C)')])
print('\nG-d  contre-epreuve, piece de doublement impaire (eta = sgn sur S3) :')
G, e = prod(GL, S3, dt, sg); fs2 = decompose(G, e, 'GL(2,3) x S3, det x sgn')
print('     -> facteur sur 4spin x 2 :', [f for f in fs2 if f in ('M4(H)', 'M8(R)', 'M8(C)')])
G, e = prod(S5, S3, sg, sg); fc2 = decompose(G, e, 'S5 x S3, sgn x sgn')
print('     -> facteur sur 4 x 2 :', [f for f in fc2 if f in ('M4(H)', 'M8(R)', 'M8(C)')])
print('\nG-e  une copie marquee (S2, eta trivial) :')
decompose(S2, one, 'S2, trivial')
G, e = prod(GL, S2, dt, one); decompose(G, e, 'GL(2,3) x S2, det x 1')
G, e = prod(S5, S2, sg, one); decompose(G, e, 'S5 x S2, sgn x 1')
