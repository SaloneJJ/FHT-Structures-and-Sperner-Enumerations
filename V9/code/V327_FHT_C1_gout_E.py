# V327 (FHT_25) — C1 : le goût E double-t-il les masses ? Réseau BBH (3 chaînes SSH, flux pi), Gamma5 = cube des frères,
# H_B = H_impair ⊗ 1 + c Gamma5 ⊗ gamma_F, insertion eps_D Gamma5 ⊗ D_F (F jouet de dimension 4, gamma_F = diag(1,1,-1,-1)). numpy seul.
import numpy as np
sx=np.array([[0,1],[1,0]],complex); sz=np.diag([1.,-1]).astype(complex); I2=np.eye(2)
def k3(a,b,c): return np.kron(np.kron(a,b),c)
def h(tf,tc,k): return np.array([[0,tf+tc*np.exp(-1j*k)],[tf+tc*np.exp(1j*k),0]])
def H(tf,tc,k): return k3(h(tf,tc,k[0]),I2,I2)+k3(sz,h(tf,tc,k[1]),I2)+k3(sz,sz,h(tf,tc,k[2]))
G5=(k3(sx,I2,I2)+k3(sz,sx,I2)+k3(sz,sz,sx))/np.sqrt(3)
tf,tc=0.4,1.0; rng=np.random.default_rng(327)
ks=[rng.uniform(-np.pi,np.pi,3) for _ in range(15)]
# commutant de {H(k)}
Mrows=[np.kron(H(tf,tc,k),np.eye(8))-np.kron(np.eye(8),H(tf,tc,k).T) for k in ks]
U,s,Vh=np.linalg.svd(np.vstack(Mrows)); null=Vh[s.size-np.sum(s<1e-9):] if np.sum(s<1e-9) else None
dimC=int(np.sum(s<1e-9)); print("dim commutant {H(k)} :",dimC)
Cs=[v.reshape(8,8) for v in Vh[-dimC:]]
print("commutant commute avec Gamma5 :",all(np.allclose(C@G5,G5@C,atol=1e-8) for C in Cs))
# F jouet : gamma_F, D_F impair, valeurs singulieres distinctes
gF=np.diag([1,1,-1,-1]).astype(complex); Y=np.array([[0.3,0.1],[0.05,0.7]]); DF=np.zeros((4,4),complex); DF[:2,2:]=Y; DF[2:,:2]=Y.T
eps=0.05; worst=1e9
for k in ks:
    Hk=H(tf,tc,k); Hev=(Hk+G5@Hk@G5)/2; Hod=Hk-Hev
    Ht=np.kron(Hod,np.eye(4))+np.kron(Hev,gF)+eps*np.kron(G5,DF)
    ev=np.sort(np.linalg.eigvalsh(Ht)); gaps=np.diff(ev)
    # chaque niveau doublement degenere : ev[0]=ev[1], ev[2]=ev[3], ...
    worst=min(worst,np.max(np.abs(ev[0::2]-ev[1::2]))) if False else worst
    pairs=np.max(np.abs(ev[0::2]-ev[1::2])); worst=min(worst,1e9); 
    if 'mx' not in dir(): mx=0
    mx=max(mx,pairs)
print("écart max entre niveaux appariés (doublement) :",mx)
# au point de Dirac K : masses du bord
K=np.array([np.pi]*3); Hk=H(tf,tc,K); Hev=(Hk+G5@Hk@G5)/2; Hod=Hk-Hev
mu=np.linalg.svd(Y,compute_uv=False)
for tfx in (tc,):   # t_f = t_c : masse de paroi nulle, on lit eps*mu
    Hk=H(tfx,tc,K); Hev=(Hk+G5@Hk@G5)/2; Hod=Hk-Hev
    Ht=np.kron(Hod,np.eye(4))+np.kron(Hev,gF)+eps*np.kron(G5,DF)
    ev=np.sort(np.abs(np.linalg.eigvalsh(Ht)))
    print("t_f=t_c, k=K : |E| =",np.round(ev,5)," attendu eps*mu =",np.round(eps*mu,5))
