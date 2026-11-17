import sys, os; _H = os.path.dirname(os.path.abspath(__file__)); [sys.path.insert(0, os.path.join(_H, "..", p)) for p in ("AF23_dirac_S3", "V352", "AF24_holonomies")]; sys.path.insert(0, _H)
import numpy as np, scipy.sparse as sp, itertools, AF23_dirac_natif_S3 as A
o = {(0,0):1,(1,1):1,(2,0):1,(3,1):1,(0,1):-1,(1,0):-1,(2,1):-1,(3,0):-1}
def Hor(N, mu, kill=None):
    A.MU=mu; G4=A.G4; H,cells=A.build(N,'S3'); A.G4=0*G4; H0,_=A.build(N,'S3'); A.G4=G4
    sg=np.repeat([o[c[1]] for c in cells],8).astype(float)
    return (H0+sp.diags(sg)@(H-H0)).tocsr(), cells
def sym(N, cells, sigma, refl):
    """g : x_a -> x_sigma(a), puis x_b -> m - x_b pour b dans refl ; renvoie P (unitaire) et le signe chi par facette"""
    m=2**N; idx={c:k for k,c in enumerate(cells)}; nc=len(cells)
    P=np.zeros((8*nc,8*nc),complex); chi={}
    for (n,f,co) in cells:
        mm=2**n; i,s=f; tA=[a for a in range(4) if a!=i]
        i2=sigma[i]; s2=(1-s) if i2 in refl else s; g=(i2,s2); tB=[a for a in range(4) if a!=i2]
        pos={}
        for k,a in enumerate(tA):
            b=sigma[a]; pos[b]=(mm-1-co[k]) if b in refl else co[k]
        co2=tuple(pos[b] for b in tB)
        pi=[tB.index(sigma[a]) for a in tA]; sg=[-1 if sigma[a] in refl else 1 for a in tA]
        V=A.transport(None,pi,sg)
        u=idx[(n,f,co)]; v=idx[(n,g,co2)]; P[8*v:8*v+8,8*u:8*u+8]=V
        chi[f]=1 if np.allclose(V@A.FG@V.conj().T,A.FG) else -1
    return P,chi
N=1
for mu in (0.0,1.0):
    H,cells=Hor(N,mu); Hd=H.toarray()
    for name,sigma,refl in [('1<->2 et x0 -> -x0',[0,2,1,3],{0}),('1<->2 seul',[0,2,1,3],set()),('x0 -> -x0 seul',[0,1,2,3],{0}),('cycle 1->2->3->1',[0,2,3,1],set())]:
        P,chi=sym(N,cells,sigma,refl)
        assert np.allclose(P@P.conj().T,np.eye(len(P)))
        R=P@Hd@P.conj().T
        print(f'mu={mu} {name:22s} : |P H P^+ - H| = {np.abs(R-Hd).max():.3e}')
print('--- recherche des signes par facette')
fac=[(i,s) for i in range(4) for s in (0,1)]
H,cells=Hor(N,0.0); Hd=H.toarray(); cf=np.repeat([fac.index(c[1]) for c in cells],8)
for name,sigma,refl in [('1<->2 et x0 -> -x0',[0,2,1,3],{0}),('cycle 1->2->3->1',[0,2,3,1],set()),('1<->2 seul',[0,2,1,3],set())]:
    P,chi=sym(N,cells,sigma,refl); best=None
    for signs in itertools.product((1,-1),repeat=7):
        sv=np.array((1,)+signs)[cf]; Ps=P*sv[None,:]
        d=np.abs(Ps@Hd@Ps.conj().T-Hd).max()
        if best is None or d<best[0]: best=(d,signs)
    print(f'{name:22s} : meilleur ecart {best[0]:.3e} ; chi par facette {chi}')
