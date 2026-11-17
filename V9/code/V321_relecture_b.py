import numpy as np
rng=np.random.default_rng(5)
Bz=np.array([[1,-1,0,0],[1,1,-2,0],[1,1,1,-3]],float);Bz=(Bz.T/np.linalg.norm(Bz,axis=1)).T
def rho(g):
    P=np.zeros((4,4))
    for i in range(4):P[g[i],i]=1
    return Bz@P@Bz.T
C=rho((0,2,3,1));Z=rho((1,0,2,3))
w,Ue=np.linalg.eig(C);Ue=Ue[:,np.argsort(np.angle(w))]
def s23(Uv):
    U=Ue.conj().T@Uv; P2=np.abs(U)**2
    k1=int(np.argmin([abs(P2[:,k].max()-2/3) for k in range(3)])); e=int(np.argmax(P2[:,k1]))
    rest=[k for k in range(3) if k!=k1]; k3=min(rest,key=lambda k:P2[e,k])
    mu,tau=[i for i in range(3) if i!=e]
    return P2[mu,k3]/(1-P2[e,k3]) if True else None
out=[]
for _ in range(400):
    A=rng.normal(size=(3,3));M=A+A.T;M=(M+Z.T@M@Z)/2
    d=rng.normal(size=3)+1j*rng.normal(size=3)
    Yn=Ue@np.diag(d)@Ue.conj().T          # complexe, C3-invariant
    mI=Yn.T@np.linalg.inv(M)@Yn; m=M*0.0+mI
    ev,V=np.linalg.eigh(m.conj().T@m); out.append(s23(V))
out=np.array(out);print("type I pur, Y_nu complexe C3-inv, M_R reel : sin2th23 in [%.3f, %.3f]"%(out.min(),out.max()))
