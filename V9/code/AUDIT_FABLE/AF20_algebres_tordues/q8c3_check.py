from burnside_signe_v2 import *
vecs=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
def mp(M): return tuple(vecs.index(((M[0][0]*a+M[0][1]*b)%3,(M[1][0]*a+M[1][1]*b)%3)) for a,b in vecs)
q8=[mp([[0,1],[2,0]]),mp([[1,1],[1,2]])]
ext=lambda p: tuple(p)+(8,9,10); c3=tuple(range(8))+(9,10,8)
G=closure([ext(q8[0]),ext(q8[1]),c3])
C,cl,idx,reps,X=char_table(G); k=len(C)
gens=[[int(x) for x in induced(G,sorted(H),eps,reps)] for H in subgroups(G) for eps in homs_pm(H)]
B=lattice_basis(gens)
e=tuple(range(11)); e0=[l for l in range(k) if reps[l]==e][0]
for a in range(k):
    chi=X[a]
    if np.allclose(chi.imag,0,atol=1e-6): continue
    b=[bb for bb in range(k) if np.allclose(X[bb],np.conj(chi),atol=1e-6)][0]
    v=np.round((chi+X[b]).real).astype(int)
    if v[e0]==4: print('V dim 4 :', [(m, in_lattice(B, m*v)) for m in (1,2,3,4)]); break
