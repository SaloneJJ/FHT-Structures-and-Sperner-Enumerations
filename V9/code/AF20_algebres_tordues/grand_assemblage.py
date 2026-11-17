"""Assemblage de la grande symetrie, gabarit de K26 : X = (S^L * B^R * G0) u (B^L * S^R * K4),
S = S_spin_min (21 sommets, Aut = GL(2,3), H1 = spinoriel 4), B = B0 (6 sommets, Aut = S5, H1 = 4 de S5) ; sur chaque cote S et B partagent un sommet (c = w)."""
import sys, os, json, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fht_tools import betti_z2, aut_order, maximal
S=json.load(open('S_spin_min.json'))
B=[['z','x0','x1','x2','x3'],['w','z'],['w','x0'],['w','x1'],['w','x2'],['w','x3']]
ren=lambda f,side: [side+('w' if x=='c' else x) for x in f]          # c (apex de S) identifie a w (apex de B)
print('B0 : Betti', betti_z2(B), '|Aut|', aut_order(sorted({x for f in B for x in f}),B))
J=[[*ren(a,'L_'),*ren(b,'R_')] for a in S for b in B]
print('joint S*B : Betti', betti_z2(J))
G0=['g0','g1','g2','g3']; K4=[list(e) for e in itertools.combinations(G0,2)]
X=[ren(a,'L_')+ren(b,'R_')+[g] for a in S for b in B for g in G0]+[ren(b,'L_')+ren(a,'R_')+e for b in B for a in S for e in K4]
X=[list(f) for f in maximal([frozenset(f) for f in X])]; V=sorted({x for f in X for x in f})
t=time.time(); bt=betti_z2(X); print('X_grand : sommets', len(V), 'facettes', len(X), 'Betti', bt, f'({time.time()-t:.0f} s)', flush=True)
t=time.time(); print('|Aut(X_grand)| =', aut_order(V,X), '; attendu (48*120)^2*24 =', (48*120)**2*24, f'({time.time()-t:.0f} s)')
json.dump(X, open('X_grand.json','w'))
