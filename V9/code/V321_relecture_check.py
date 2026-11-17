import numpy as np
rng=np.random.default_rng(321)
v1=np.array([3,-1,-1,-1])/np.sqrt(12); v2=np.array([0,1,1,-2])/np.sqrt(6); v3=np.array([0,1,-1,0])/np.sqrt(2)
B=np.array([v1,v2,v3]).T
def perm(p):
    P=np.zeros((4,4))
    for i in range(4): P[p[i],i]=1
    return B.T@P@B
R=perm((0,2,3,1))            # 3-cycle fixant g0
sgnR=1                       # pair : 3 et 3' coincident sur C3
# T1a : seule hypothese H=YY^T reel symetrique, C3-invariant (e_R quelconque)
worst=0
for _ in range(5000):
    A=rng.normal(size=(3,3)); H=A@A.T
    H=sum(np.linalg.matrix_power(R,k)@H@np.linalg.matrix_power(R,k).T for k in range(3))/3
    ev=np.sort(np.linalg.eigvalsh(H)); g=min(np.diff(ev)); worst=max(worst,g)
print("T1a H reel sym C3-inv : max sur 5000 de l'ecart minimal =",worst)
# T1b : Y in Hom(3_R,3_L) equivariant sous C3 avec e_R en 3' (S4) -> sur C3 identique ; e_R = trois singulets reels triviaux:
for _ in range(3):
    A=rng.normal(size=(3,3))
    Y=sum(np.linalg.matrix_power(R,k)@A for k in range(3))/3   # e_R triviaux : rho_L Y = Y
    print("T1b e_R singulets triviaux, Y reel : val. sing.",np.round(np.linalg.svd(Y,compute_uv=False),4))
# T1c : Jarlskog avec tout reel
def Jinv(Ha,Hb): return np.linalg.det(Ha@Hb-Hb@Ha)
A=rng.normal(size=(3,3));Hu=A@A.T; A=rng.normal(size=(3,3));Hd=A@A.T
print("T1c det[Hu,Hd] reels =",Jinv(Hu,Hd))
A=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));Hc=A@A.conj().T
print("     det[Hc,Hd] (un complexe) =",Jinv(Hc,Hd))
# T1d : TM1 + delta in {0,pi} -> sin^2 th23 via la relation de somme
s13=0.022; s=np.sqrt(s13)
for cd in (1,-1):
    cot2=-cd*2*np.sqrt(2)*s*np.sqrt(1-3*s13)/(1-5*s13)
    t=0.5*np.arctan2(1,cot2); print("T1d cos d=%+d : sin^2 th23 = %.3f"%(cd,np.sin(t)**2))
# T2 : poids de w=(e0-ej)/sqrt2 sur v1^(k) (direction fixe par Stab(g_k))
for j in (1,2,3):
    w=np.zeros(4); w[0]=1; w[j]=-1; w/=np.sqrt(2)
    out=[]
    for k in range(4):
        u=-np.ones(4); u[k]=3; u/=np.sqrt(12); out.append(round(np.dot(u,w)**2,4))
    print(f"T2 arete (0 {j}) : poids sur v1^(k), k=0..3 :",out)
# T2 GJ realiste en base de chaine, v1=tau
r=1/206.77
E=1.0; F=np.sqrt(r)*E   # m_e ~ F^2/E
Y=np.array([[16.8,0,0],[0,E,F],[0,F,0]])
ev,U=np.linalg.eigh(Y@Y.T); U=U[:,np.argsort(ev)]   # colonnes e, mu, tau en base (v1,v2,v3)
for j in (1,3):
    w=np.zeros(4); w[0]=1; w[j]=-1; w/=np.sqrt(2); wc=B.T@w
    print(f"T2 GJ v1=tau, (0 {j}) : |U_a1|^2 (e,mu,tau) =",np.round(np.abs(U.T@wc)**2,4))
# echappatoire : e ~ v1 dans le bloc GJ, tau sur le plan ; e = cos t v1 - sin t u, u dans le plan
t=np.arctan(np.sqrt(r)); w=np.zeros(4); w[0]=1; w[1]=-1; w/=np.sqrt(2); wc=B.T@w
vals=[]
for phi in np.linspace(0,np.pi,181):
    u=np.array([0,np.cos(phi),np.sin(phi)]); e=np.array([np.cos(t),0,0])-np.sin(t)*u
    vals.append(np.dot(e,wc)**2)
vals=np.array(vals); s12=1-vals/(1-s13)
print("Echapp. e~v1 (GJ) : |U_e1|^2 in [%.4f, %.4f] ; sin^2 th12 in [%.3f, %.3f] (TM1 : 0.318)"%(vals.min(),vals.max(),s12.min(),s12.max()))
# etiquetage different : colonne fixee = nu3 avec v1 = tau, (0 3)
w=np.zeros(4); w[0]=1; w[3]=-1; w/=np.sqrt(2); wc=B.T@w
print("Colonne fixee lue comme nu3, (0 3), GJ : |U_a3|^2 (e,mu,tau)=",np.round(np.abs(U.T@wc)**2,4))
