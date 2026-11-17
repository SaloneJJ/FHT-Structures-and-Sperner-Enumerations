import sys, os; _H = os.path.dirname(os.path.abspath(__file__)); [sys.path.insert(0, os.path.join(_H, "..", p)) for p in ("AF23_dirac_S3", "V352")]; os.chdir(_H)
import numpy as np, scipy.sparse as sp, itertools, AF23_dirac_natif_S3 as A, sys
o = {(0,0):1,(1,1):1,(2,0):1,(3,1):1,(0,1):-1,(1,0):-1,(2,1):-1,(3,0):-1}
fac=[(i,s) for i in range(4) for s in (0,1)]
def Hor(N, mu):
    A.MU=mu; G4=A.G4; H,cells=A.build(N,'S3'); A.G4=0*G4; H0,_=A.build(N,'S3'); A.G4=G4
    sg=np.repeat([o[c[1]] for c in cells],8).astype(float)
    return (H0+sp.diags(sg)@(H-H0)).toarray(), cells
X,Y=A.Xg,A.Yg; G5=A.G5; Am=(X[0]-X[1])/np.sqrt(2); Bm=(X[0]+X[1]-2*X[2])/np.sqrt(6)
gens=[G5,Am,Bm,Y[0],Y[1],Y[2]]; names=['G5','A','B','Y1','Y2','Y3']
mons=[]
for sub in itertools.product((0,1),repeat=6):
    M=np.eye(8,dtype=complex)
    for b,g in zip(sub,gens):
        if b: M=M@g
    mons.append(('.'.join(n for b,n in zip(sub,names) if b) or '1',M))
pairs=[(a,b) for a in range(8) for b in range(a+1,8) if fac[a][0]!=fac[b][0]]
def apply(H,t,cf):
    Sg=np.ones((8,8))
    for (a,b),x in zip(pairs,t): Sg[a,b]=Sg[b,a]=x
    return H*Sg[cf][:,cf]
def search(H,cells):
    nc=len(cells); cf=np.array([fac.index(c[1]) for c in cells]); found=[]
    Hb=H.reshape(nc,8,nc,8).transpose(0,2,1,3)
    same=cf[:,None]==cf[None,:]; nz=np.abs(Hb).sum((2,3))>1e-12
    for name,U in mons:
        R=np.einsum('ij,uvjk,lk->uvil',U,Hb.conj(),U.conj())
        if not np.allclose(R[same&nz],-Hb[same&nz],atol=1e-10): continue
        # signes requis entre facettes
        req={}; ok=True
        for u,v in zip(*np.nonzero(~same&nz)):
            r=R[u,v]; h=-Hb[u,v]; k=np.argmax(np.abs(h)); e=(r.ravel()[k]/h.ravel()[k])
            if not (np.isclose(abs(e),1) and np.allclose(r,e*h,atol=1e-10) and np.isclose(e.imag,0)): ok=False; break
            key=(cf[u],cf[v]); e=int(np.sign(e.real))
            if req.setdefault(key,e)!=e: ok=False; break
        if not ok: continue
        eps={0:1}; st=[0]
        while st and ok:
            a=st.pop()
            for (x,y),e in req.items():
                if x!=a: continue
                w=eps[a]*e
                if y in eps:
                    if eps[y]!=w: ok=False
                else: eps[y]=w; st.append(y)
        if ok: found.append((name,int(np.round((U@U.conj())[0,0].real))))
    return found
N=int(sys.argv[1]) if len(sys.argv)>1 else 1
H,cells=Hor(N,1.0); cf=np.repeat([fac.index(c[1]) for c in cells],8)
confs=[('AF23 (SVD)',[1]*len(pairs))]+[(f'aleatoire {s}',np.random.default_rng(s).choice([-1,1],len(pairs))) for s in range(4)]
for name,t in confs:
    Ht=apply(H,t,cf); f=search(Ht,cells)
    E=np.sort(np.linalg.eigvalsh(Ht)); dbl=np.allclose(E[0::2],E[1::2],atol=1e-9)
    print(f'N={N} {name:12s} : J trouves (U, J^2) = {f} ; spectre symetrique E <-> -E : {np.allclose(E,-E[::-1],atol=1e-9)}',flush=True)
