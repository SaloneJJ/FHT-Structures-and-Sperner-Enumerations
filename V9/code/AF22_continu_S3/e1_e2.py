import itertools, numpy as np, json, time
from scipy.spatial import ConvexHull
import scipy.sparse as sp, scipy.sparse.linalg as sla
r=1/np.sqrt(2); P=set()
for i in range(4):
    for s in (1,-1): v=[0]*4; v[i]=s; P.add(tuple(v))
for s in itertools.product((0.5,-0.5),repeat=4): P.add(s)
for i,j in itertools.combinations(range(4),2):
    for a in (r,-r):
        for b in (r,-r): v=[0.0]*4; v[i]=a; v[j]=b; P.add(tuple(v))
X=np.array(sorted(P)); H=ConvexHull(X); T0=[tuple(sorted(int(x) for x in f)) for f in H.simplices]
def bary(tets):
    """subdivision barycentrique : sommets = simplexes ; tetraedres = drapeaux complets"""
    S=set()
    for t in tets:
        for k in range(1,5): S.update(itertools.combinations(t,k))
    idx={s:i for i,s in enumerate(sorted(S))}
    new=[]
    for t in tets:
        for perm in itertools.permutations(t):
            new.append(tuple(sorted(idx[tuple(sorted(perm[:k]))] for k in range(1,5))))
    return new, len(idx)
def graph_lap(simplices, n):
    E=set()
    for t in simplices:
        for a,b in itertools.combinations(t,2): E.add((min(a,b),max(a,b)))
    E=np.array(sorted(E)); I=np.r_[E[:,0],E[:,1]]; J=np.r_[E[:,1],E[:,0]]
    A=sp.csr_matrix((np.ones(len(I)),(I,J)),shape=(n,n)); deg=np.array(A.sum(1)).ravel()
    return sp.diags(deg)-A, A, deg
def analyse(name, L, A, deg, n, k=40):
    w=np.sort(sla.eigsh(L.tocsc(),k=k,sigma=-1e-3,which='LM',return_eigenvectors=False)) if n>3000 else np.sort(np.linalg.eigvalsh(L.toarray()))[:k]
    w=w[1:]
    # paquets par egalite stricte (1e-6 relatif)
    packs=[]; cur=[w[0]]
    for x in w[1:]:
        if abs(x-cur[-1])<=1e-6*abs(x): cur.append(x)
        else: packs.append(cur); cur=[x]
    packs.append(cur)
    W=0.5*(sp.identity(n)+sp.diags(1/deg)@A); rng=np.random.default_rng(0); st=rng.choice(n,size=min(50,n),replace=False)
    v=np.zeros((n,len(st))); v[st,np.arange(len(st))]=1; ts=np.unique(np.round(np.logspace(0,np.log10(min(3000,10*n)),30)).astype(int)); Pt=[]; tt=0
    for t in ts:
        while tt<t: v=W@v; tt+=1
        Pt.append(np.mean(v[st,np.arange(len(st))]))
    ds=-2*np.gradient(np.log(Pt),np.log(ts))
    print(f'{name}: n={n}  premiers paquets (taille, valeur) :', [(len(p), round(float(p[0]),5)) for p in packs[:8]], flush=True)
    print('   d_s :', [(int(a),round(float(b),2)) for a,b in zip(ts,ds)][::3], flush=True)
    return packs
out={}
# E1 : barycentrique
T=T0; n=48
for lev in (1,2):
    T,n=bary(T); L,A,deg=graph_lap(T,n); p=analyse(f'E1 bary niveau {lev}', L,A,deg,n)
    out[f'E1_{lev}']=[(len(q),float(q[0])) for q in p[:8]]
# E2 : bord du 4-cube, grille de cote m = 2^n
for lev in (2,3,4,5):
    m=2**lev; pts=[p for p in itertools.product(range(m+1),repeat=4) if any(c in (0,m) for c in p)]
    idx={p:i for i,p in enumerate(pts)}; I=[];J=[]
    for p in pts:
        for d in range(4):
            q=list(p); q[d]+=1; q=tuple(q)
            if q in idx:
                # arete sur le bord : les deux extremites sur une meme face (une coordonnee commune egale a 0 ou m, autre que d)
                if any((p[e] in (0,m)) and e!=d for e in range(4)): I.append(idx[p]); J.append(idx[q])
    n=len(pts); A=sp.csr_matrix((np.ones(2*len(I)),(I+J,J+I)),shape=(n,n)); deg=np.array(A.sum(1)).ravel(); L=sp.diags(deg)-A
    p=analyse(f'E2 bord du 4-cube, cote {m}', L, A, deg, n)
    out[f'E2_{m}']=[(len(q),float(q[0])) for q in p[:8]]
json.dump(out, open('e1_e2_resultats.json','w'), indent=1)
