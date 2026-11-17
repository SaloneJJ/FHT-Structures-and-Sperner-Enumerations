import numpy as np, itertools, time, scipy.sparse.linalg as sla
import AF23_dirac_natif_S3 as A
out=[]
def P(*a):
    s=' '.join(str(x) for x in a); print(s, flush=True); out.append(s)
# (a) localisation des plus bas modes (N = 3) : poids pres des aretes de la racine (>= 2 coordonnees au bord de la facette) et par niveau
N=3; H,cells=A.build(N,'S3'); m=2**N
E,U=sla.eigsh(H.tocsc(),k=24,sigma=1e-6,which='LM')
o=np.argsort(np.abs(E)); lev=np.array([c[0] for c in cells])
def nb(c): return sum(1 for x in c[2] if x in (0,2**c[0]-1))
fin=lev==N; edge2=np.array([nb(c)>=2 for c in cells]); face1=np.array([nb(c)==1 for c in cells]); bulk=np.array([nb(c)==0 for c in cells])
vol=[np.mean(edge2[fin]),np.mean(face1[fin]),np.mean(bulk[fin])]
P('(a) N = 3 ; fractions de volume (niveau fin) : pres des aretes de la racine %.3f, pres des faces %.3f, interieur %.3f'%tuple(vol))
for k in o[:16]:
    w=(np.abs(U[:,k])**2).reshape(-1,8).sum(1); wf=w[fin].sum()
    P('   E = %+.5f : niveau fin %.2f ; dans le niveau fin : aretes %.3f, faces %.3f, interieur %.3f'%(E[k],wf,w[fin&edge2].sum()/wf,w[fin&face1].sum()/wf,w[fin&bulk].sum()/wf))
# (b) mu = 0
A.MU=0.0
for N in (2,3):
    H,cells=A.build(N,'S3')
    if N==2: Ev=np.linalg.eigvalsh(H.toarray())
    else: Ev=sla.eigsh(H.tocsc(),k=48,sigma=1e-6,which='LM',return_eigenvectors=False)
    a=np.sort(np.abs(Ev))[:40]
    P('(b) mu = 0, N = %d : plus petites |E| :'%N, np.round(a[:28],5).tolist())
A.MU=1.0
# (c) masse de gout invariante sous O ?
Xg,Yg,G4,G5,FG=A.Xg,A.Yg,A.G4,A.G5,A.FG
basis=[]
gens=[G5,G4]+list(Yg)
# commutant : resoudre [M, g] = 0 pour g dans gens (M 8x8 complexe)
rows=[np.kron(np.eye(8),g.T)-np.kron(g,np.eye(8)) for g in gens]
Am=np.vstack(rows); _,s,Vh=np.linalg.svd(Am); ns=Vh[np.sum(s>1e-9):].conj()
C=[v.reshape(8,8) for v in ns]
P('(c) commutant de {G5, G4, Y_1, Y_2, Y_3} : dimension', len(C), '; FG y est :', np.linalg.matrix_rank(np.array([c.ravel() for c in C]+[FG.ravel()]))==len(C))
# rotations de 90 degres autour de chaque axe c : permutation (a,b)->(b,-a) des deux autres indices
rots=[]
for c in range(3):
    a,b=[x for x in range(3) if x!=c]
    cmap=[0,1,2]; cmap[a]=b; cmap[b]=a; sig=[1,1,1]; sig[b]=-1     # e_a -> e_b, e_b -> -e_a
    rots.append(A.transport(None,cmap,sig))
rows=[]
for V in rots:
    rows.append(np.kron(V,V.conj())-np.eye(64))       # vec(V M V^+) - vec(M)  (vectorisation par lignes)
# M = sum x_i C_i : contraintes
Cm=np.array([c.ravel() for c in C]).T
Bm=np.vstack([R@Cm for R in rows]); _,s2,Vh2=np.linalg.svd(Bm); k=Cm.shape[1]-np.sum(s2>1e-9)
P('(c) invariants sous O (rotations de 90 degres) dans ce commutant : dimension', k, '(1 = scalaires seulement)')
FGr=[np.allclose(V@FG@V.conj().T,FG) for V in rots]; P('    FG invariant sous les rotations de 90 degres :', FGr)
# (d) parite des degres d'aretes : 288-cellules (2O) et ses subdivisions ; bord du 4-cube
from scipy.spatial import ConvexHull
r=1/np.sqrt(2); Pp=set()
for i in range(4):
    for s in (1,-1): v=[0]*4; v[i]=s; Pp.add(tuple(v))
for s in itertools.product((0.5,-0.5),repeat=4): Pp.add(s)
for i,j in itertools.combinations(range(4),2):
    for a in (r,-r):
        for b in (r,-r): v=[0.0]*4; v[i]=a; v[j]=b; Pp.add(tuple(v))
X=np.array(sorted(Pp)); T0=[tuple(sorted(int(x) for x in f)) for f in ConvexHull(X).simplices]
def edge_degrees(tets):
    from collections import Counter
    c=Counter()
    for t in tets:
        for e in itertools.combinations(t,2): c[e]+=1
    return c
deg=edge_degrees(T0); from collections import Counter
P('(d) 288-cellules (2O) : degres des aretes (nombre de tetraedres autour) :', dict(Counter(deg.values())))
