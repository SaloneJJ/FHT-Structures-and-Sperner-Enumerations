"""Raffinements hierarchiques de l'orbite 2O dans S^3 : T0 = bord de l'enveloppe convexe, T_n = subdivisions arretieres ; tests (a)-(d)."""
import itertools, numpy as np, json, time, sys
from scipy.spatial import ConvexHull
import scipy.sparse as sp, scipy.sparse.linalg as sla
r=1/np.sqrt(2)
P=set()
for i in range(4):
    for s in (1,-1):
        v=[0]*4; v[i]=s; P.add(tuple(v))
for s in itertools.product((0.5,-0.5),repeat=4): P.add(s)
for i,j in itertools.combinations(range(4),2):
    for a in (r,-r):
        for b in (r,-r):
            v=[0.0]*4; v[i]=a; v[j]=b; P.add(tuple(v))
X=np.array(sorted(P)); print('points', len(X))
H=ConvexHull(X); tets=[tuple(sorted(f)) for f in H.simplices]
print('T0 : facettes (tetraedres)', len(tets))
def subdivide(tets):
    mid={}; new=[]
    def m(a,b):
        k=(min(a,b),max(a,b))
        if k not in mid: mid[k]=nv[0]; nv[0]+=1
        return mid[k]
    for (a,b,c,d) in tets:
        ab,ac,ad,bc,bd,cd=m(a,b),m(a,c),m(a,d),m(b,c),m(b,d),m(c,d)
        new+=[(a,ab,ac,ad),(b,ab,bc,bd),(c,ac,bc,cd),(d,ad,bd,cd)]
        new+=[(ac,bd,ab,ad),(ac,bd,ad,cd),(ac,bd,cd,bc),(ac,bd,bc,ab)]
    return [tuple(sorted(t)) for t in new]
def betti_Q(tets):
    S=set()
    for t in tets:
        for k in range(1,5): S.update(itertools.combinations(t,k))
    byd={}
    for s in S: byd.setdefault(len(s)-1,[]).append(s)
    byd[-1]=[()]
    idx={k:{s:i for i,s in enumerate(sorted(v))} for k,v in byd.items()}
    rk={}
    for k in range(0,4):
        I,J,Vv=[],[],[]
        for j,s in enumerate(sorted(byd[k])):
            for i in range(len(s)): I.append(idx[k-1][s[:i]+s[i+1:]]); J.append(j); Vv.append((-1)**i)
        M=np.zeros((len(byd[k-1]),len(byd[k]))); M[I,J]=Vv; rk[k]=np.linalg.matrix_rank(M)
    rk[4]=0
    return {k:len(byd[k])-rk.get(k,0)-rk[k+1] for k in range(-1,4) if len(byd[k])-rk.get(k,0)-rk[k+1]}
def laplacian(tets,n):
    E=set()
    for t in tets:
        for a,b in itertools.combinations(t,2): E.add((a,b))
    I=[a for a,b in E]+[b for a,b in E]; J=[b for a,b in E]+[a for a,b in E]
    A=sp.csr_matrix((np.ones(len(I)),(I,J)),shape=(n,n)); deg=np.array(A.sum(1)).ravel()
    return sp.diags(deg)-A, A, deg, len(E)
res={}
levels=[tets]; nv=[len(X)]; nvs=[len(X)]
for lev in range(1,4):
    levels.append(subdivide(levels[-1])); nvs.append(nv[0])
for lev,T in enumerate(levels):
    n=nvs[lev]; t0=time.time()
    L,A,deg,ne=laplacian(T,n)
    out={'sommets':n,'aretes':ne,'tetraedres':len(T)}
    if lev<=1: out['betti_Q']={int(a):int(b) for a,b in betti_Q(T).items()}
    k=min(60,n-2)
    w=np.sort(sla.eigsh(L.tocsc(),k=k,sigma=-1e-3,which='LM',return_eigenvectors=False)) if n>2000 else np.sort(np.linalg.eigvalsh(L.toarray()))[:k]
    out['bas_spectre']=[round(float(x),5) for x in w]
    # dimension spectrale : marche paresseuse W = (I + D^-1 A)/2
    Dinv=sp.diags(1/deg); W=0.5*(sp.identity(n)+Dinv@A)
    rng=np.random.default_rng(0); starts=rng.choice(n,size=min(50,n),replace=False)
    vec=np.zeros((n,len(starts))); vec[starts,np.arange(len(starts))]=1.0
    ts=np.unique(np.round(np.logspace(0,np.log10(min(4000,20*n)),40)).astype(int)); Pt={}; cur=vec.copy(); tt=0
    for t in ts:
        while tt<t: cur=W@cur; tt+=1
        Pt[int(t)]=float(np.mean(cur[starts,np.arange(len(starts))]))
    tl=np.array(sorted(Pt)); pl=np.array([Pt[t] for t in tl])
    ds=-2*np.gradient(np.log(pl),np.log(tl))
    out['ds']=[(int(a),round(float(b),3)) for a,b in zip(tl,ds)]
    out['duree_s']=round(time.time()-t0,1)
    res[lev]=out
    print(lev, n, ne, len(T), out.get('betti_Q'), 'bas spectre', out['bas_spectre'][:16], flush=True)
    print('   d_s(t) :', out['ds'][::3], flush=True)
json.dump(res, open("s3_2O_resultats.json","w"), indent=1, default=lambda o: int(o) if hasattr(o,"item") else str(o))
