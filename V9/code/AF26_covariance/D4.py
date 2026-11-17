import sys, os; _H = os.path.dirname(os.path.abspath(__file__)); [sys.path.insert(0, os.path.join(_H, "..", p)) for p in ("AF23_dirac_S3", "V352", "AF24_holonomies")]; sys.path.insert(0, _H)
import numpy as np, scipy.sparse as sp, AF23_dirac_natif_S3 as A, time
o = {(0,0):1,(1,1):1,(2,0):1,(3,1):1,(0,1):-1,(1,0):-1,(2,1):-1,(3,0):-1}
fac=[(i,s) for i in range(4) for s in (0,1)]
def Hor(N, mu):
    A.MU=mu; G4=A.G4; H,cells=A.build(N,'S3'); A.G4=0*G4; H0,_=A.build(N,'S3'); A.G4=G4
    sg=np.repeat([o[c[1]] for c in cells],8).astype(float)
    return (H0+sp.diags(sg)@(H-H0)).toarray(), cells
N=2; t0=time.time()
Hs={mu:Hor(N,mu) for mu in (1.0,0.0)}; cells=Hs[1.0][1]; cf=np.repeat([fac.index(c[1]) for c in cells],8)
pairs=[(a,b) for a in range(8) for b in range(a+1,8) if fac[a][0]!=fac[b][0]]
def apply(H,t):
    S=np.ones((8,8))
    for (a,b),x in zip(pairs,t): S[a,b]=S[b,a]=x
    return H*S[cf][:,cf]
def run(name,t):
    out=[]
    for mu in (1.0,0.0):
        E=np.sort(np.abs(np.linalg.eigvalsh(apply(Hs[mu][0],t)))); out.append(E)
    print(f'{name:22s} mu=1 : |E|min {out[0][0]:.4f}, n_leger {int(np.sum(out[0]<0.3586))}, 6 plus petites {np.round(out[0][:6],4).tolist()} ; mu=0 : |E|min {out[1][0]:.4f}  ({time.time()-t0:.0f} s)',flush=True)
    return out
ref=run('AF23 (signes SVD)',[1]*len(pairs))
g=[-1 if 0 in p else 1 for p in pairs]; ctl=run('jauge facette (0,0)',g)
print('   controle de jauge : ecart',max(np.abs(ref[0]-ctl[0]).max(),np.abs(ref[1]-ctl[1]).max()))
for seed in range(4):
    rng=np.random.default_rng(seed); run(f'aleatoire {seed}',rng.choice([-1,1],len(pairs)))
