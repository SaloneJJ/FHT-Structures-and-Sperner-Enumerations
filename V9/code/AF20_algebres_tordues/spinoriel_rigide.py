"""Rigidification (point 2) appliquee au point 3 : on recolle a X le long de l'arbre acyclique Z = c * {8 vecteurs}
un gadget contractile Y = c * Gamma ; Gamma : pour chaque base ordonnee (v, w), un sommet b relie a v (bras de longueur 1), a w (longueur 2) et a v + w (longueur 3) : code la structure lineaire.
Lemme : Y n X = Z acyclique, Y contractile => H~(X u Y) = H~(X) (Mayer-Vietoris)."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fht_tools import betti_z2, aut_order, maximal
vecs=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
neg=lambda v: ((-v[0])%3,(-v[1])%3)
line=lambda v: min(vecs.index(v), vecs.index(neg(v)))
X=[[f'v{i}', f'l{line(v)}'] for i,v in enumerate(vecs)]+[[f'v{i}','c'] for i in range(8)]
G=[]
for i,v in enumerate(vecs):
    for j,w in enumerate(vecs):
        if (v[0]*w[1]-v[1]*w[0])%3!=0:
            b,m,n1,n2=f'b{i}_{j}',f'm{i}_{j}',f'n{i}_{j}',f'o{i}_{j}'; u=vecs.index(((v[0]+w[0])%3,(v[1]+w[1])%3)); G+=[[b,f'v{i}'],[b,m],[m,f'v{j}'],[b,n1],[n1,n2],[n2,f'v{u}']]
Y=[e+['c'] for e in G]
Xp=maximal([frozenset(f) for f in X+Y])
V=sorted(set(x for f in Xp for x in f))
print('X : Betti', betti_z2(X), '| X\' : sommets', len(V), 'facettes', len(Xp), 'Betti', betti_z2([list(f) for f in Xp]))
print('|Aut(X\')| =', aut_order(V,[list(f) for f in Xp]), '(attendu 48 = |GL(2,3)|)')
