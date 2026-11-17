"""S_spin rigidifie au plus court : gadget = les 8 droites affines de AG(2,3) ne passant pas par 0 (triplets v+w+u=0),
chacune un sommet t relie a ses 3 points, cone sur c (Y = c * Gamma, Y n X = etoile de c, acyclique)."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fht_tools import betti_z2, aut_order, maximal
vecs=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
neg=lambda v: ((-v[0])%3,(-v[1])%3)
line=lambda v: min(vecs.index(v), vecs.index(neg(v)))
X=[[f'v{i}', f'l{line(v)}'] for i,v in enumerate(vecs)]+[[f'v{i}','c'] for i in range(8)]
T=sorted({tuple(sorted(t)) for t in itertools.combinations(range(8),3) if all((sum(vecs[i][k] for i in t))%3==0 for k in range(2))})
print('triplets a somme nulle :', len(T))
Y=[[f't{n}', f'v{i}', 'c'] for n,t in enumerate(T) for i in t]
S=maximal([frozenset(f) for f in X+Y]); V=sorted(set(x for f in S for x in f)); F=[list(f) for f in S]
print('sommets', len(V), 'facettes', len(F), 'Betti Z/2', betti_z2(F), '|Aut|', aut_order(V,F))
import json; json.dump(F, open('S_spin_min.json','w'))
