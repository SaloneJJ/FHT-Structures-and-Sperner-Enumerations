"""Point 3 : un petit complexe ou le centre z = -I de GL(2,3) agit par -1 sur l'homologie.
X = (8 vecteurs non nuls de F_3^2) + (4 droites) + apex c ; aretes v--droite(v), c--v.  H1 = fonctions impaires sur les vecteurs (dim 4)."""
from twisted_alg import *
import numpy as np, itertools, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fht_tools import betti_z2, aut_order
vecs=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
def mp(M): return tuple(vecs.index(((M[0][0]*a+M[0][1]*b)%3,(M[1][0]*a+M[1][1]*b)%3)) for a,b in vecs)
GL=closure([mp([[1,1],[0,1]]),mp([[0,1],[2,0]]),mp([[2,0],[0,1]])])
line=lambda v: min(vecs.index(v), vecs.index(((-v[0])%3,(-v[1])%3)))
lines=sorted(set(line(v) for v in vecs))
F=[[f'v{i}', f'l{line(v)}'] for i,v in enumerate(vecs)]+[[f'v{i}','c'] for i in range(8)]
V=sorted(set(x for f in F for x in f))
print('sommets',len(V),'aretes',len(F),'Betti Z/2',betti_z2(F),'|Aut|',aut_order(V,F))
# caractere de H1 (fonctions impaires) sur GL(2,3)
C,cl,idx,reps,X=char_table(GL); k=len(C); sizes=np.array([len(c) for c in C])
neg=lambda i: vecs.index(((-vecs[i][0])%3,(-vecs[i][1])%3))
def chi_odd(p): return (sum(1 for i in range(8) if p[i]==i) - sum(1 for i in range(8) if p[i]==neg(i)))/2
chi=np.array([chi_odd(r) for r in reps])
# H1 du graphe = (aretes - sommets + 1) en caracteres : verifie par la formule directe
mult=[np.round(np.sum(sizes*chi*np.conj(X[j]))/len(GL),6) for j in range(k)]
e0=[l for l in range(k) if reps[l]==tuple(range(8))][0]
zc=[l for l in range(k) if reps[l]==mp([[2,0],[0,2]])][0]
print('decomposition de H1 en irreductibles de GL(2,3) (degre, valeur du centre, multiplicite) :')
for j in range(k):
    if abs(mult[j])>1e-6: print('   degre',int(round(X[j][e0].real)),' z ->',np.round(X[j][zc].real/X[j][e0].real,3),' mult',mult[j])
