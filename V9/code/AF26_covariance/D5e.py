import sys, os; _H = os.path.dirname(os.path.abspath(__file__)); [sys.path.insert(0, os.path.join(_H, "..", p)) for p in ("AF23_dirac_S3", "V352")]; os.chdir(_H)
import numpy as np, sys
exec(open('D5c.py').read().split("def search")[0])
def search_cells(H,cells,sign=-1):
    nc=len(cells); Hb=H.reshape(nc,8,nc,8).transpose(0,2,1,3); nz=np.abs(Hb).sum((2,3))>1e-12; found=[]
    diag=np.eye(nc,dtype=bool)
    for name,U in mons:
        R=np.einsum('ij,uvjk,lk->uvil',U,Hb,U.conj())
        if not np.allclose(R[diag&nz],sign*Hb[diag&nz],atol=1e-10): continue
        req={}; ok=True
        for u,v in zip(*np.nonzero(~diag&nz)):
            r=R[u,v]; h=sign*Hb[u,v]; k=np.argmax(np.abs(h)); e=r.ravel()[k]/h.ravel()[k]
            if not (np.isclose(abs(e),1) and np.isclose(e.imag,0,atol=1e-9) and np.allclose(r,e.real*h,atol=1e-10)): ok=False; break
            req.setdefault(u,[]).append((v,int(np.sign(e.real))))
        if not ok: continue
        eps={0:1}; st=[0]
        while st and ok:
            a=st.pop()
            for b,e in req.get(a,[]):
                w=eps[a]*e
                if b in eps:
                    if eps[b]!=w: ok=False; break
                else: eps[b]=w; st.append(b)
        if ok and len(eps)==nc: found.append((name,int(np.round((U@U)[0,0].real))))
    return found
N=int(sys.argv[1]); H,cells=Hor(N,1.0); cf=np.repeat([fac.index(c[1]) for c in cells],8)
for name,t in [('AF23 (SVD)',[1]*len(pairs))]+[(f'aleatoire {s}',np.random.default_rng(s).choice([-1,1],len(pairs))) for s in range(4)]:
    Ht=apply(H,t,cf); E=np.sort(np.linalg.eigvalsh(Ht))
    print(f'N={N} {name:12s} : Gamma (unitaire, GH=-HG) = {search_cells(Ht,cells,-1)} ; sym. unitaire (UH=HU) = {search_cells(Ht,cells,+1)} ; E<->-E : {np.allclose(E,-E[::-1],atol=1e-9)}',flush=True)
