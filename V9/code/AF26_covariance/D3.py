import sys, os; _H = os.path.dirname(os.path.abspath(__file__)); [sys.path.insert(0, os.path.join(_H, "..", p)) for p in ("AF23_dirac_S3", "V352", "AF24_holonomies")]; sys.path.insert(0, _H)
import numpy as np, scipy.sparse as sp, time, AF23_dirac_natif_S3 as A
o = {(0,0):1,(1,1):1,(2,0):1,(3,1):1,(0,1):-1,(1,0):-1,(2,1):-1,(3,0):-1}
fac=[(i,s) for i in range(4) for s in (0,1)]
def Hor(N, mu):
    A.MU=mu; G4=A.G4; H,cells=A.build(N,'S3'); A.G4=0*G4; H0,_=A.build(N,'S3'); A.G4=G4
    sg=np.repeat([o[c[1]] for c in cells],8).astype(float)
    return (H0+sp.diags(sg)@(H-H0)).tocsr(), cells
N=2; t=time.time()
H1,cells=Hor(N,1.0); H0,_=Hor(N,0.0); HF=(H1-H0).tocsr()
res={}
for name,U in [('A-',[1]*8),('B~',[1,1,0,0,0,0,1,1]),('C~',[0,1,1,1,1,1,1,1]),('D~',[0,0,1,1,1,1,1,1])]:
    eps=np.repeat([1-2*U[fac.index(c[1])] for c in cells],8).astype(float)
    H=(H0+sp.diags(eps)@HF).toarray(); assert np.allclose(H,H.conj().T)
    E=np.sort(np.abs(np.linalg.eigvalsh(H))); res[name]=E
    print(name,'n_leger =',int(np.sum(E<0.3586)),'; 12 plus petites |E| :',np.round(E[:12],4).tolist(),f'({time.time()-t:.0f} s)',flush=True)
print('controle B~ = A :',np.max(np.abs(res['A-'][:40]-res['B~'][:40])))
