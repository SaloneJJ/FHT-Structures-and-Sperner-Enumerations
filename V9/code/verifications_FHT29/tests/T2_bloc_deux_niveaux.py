# T2 (Fable, 16/11/2026) -- Les trois generations sont-elles les moments du niveau parent ?
# Reseau BBH uniforme (trois chaines SSH, flux pi) sur un bloc periodique 4^3 = huit octets. numpy seul, quelques secondes.
# Criteres : CRITERES_T1_T2.md (definition 1).
import itertools, numpy as np
def Hbloc(tf,tc,L=4):
    S=list(itertools.product(range(L),repeat=3)); ix={s:i for i,s in enumerate(S)}; H=np.zeros((L**3,L**3))
    for s in S:
        for d in range(3):
            t=list(s); t[d]=(t[d]+1)%L
            amp=(tf if s[d]%2==0 else tc)*(-1)**sum(s[:d])
            H[ix[tuple(t)],ix[s]]+=amp; H[ix[s],ix[tuple(t)]]+=amp
    return H,S,ix
sx=np.array([[0,1],[1,0]],complex); sz=np.diag([1.,-1]).astype(complex); I2=np.eye(2); k3=lambda a,b,c:np.kron(np.kron(a,b),c)
Xg=[k3(sx,I2,I2),k3(sz,sx,I2),k3(sz,sz,sx)]
for tf,tc in ((1.0,1.0),(0.7,1.0),(0.4,1.0)):
    H,S,ix=Hbloc(tf,tc); ev=np.linalg.eigvalsh(H)
    lev={}
    for e in np.round(np.abs(ev),6): lev[e]=lev.get(e,0)+1
    print(f"t_f={tf}, t_c={tc} : niveaux |E| (multiplicite) :", "  ".join(f"{e:.4f} ({m})" for e,m in sorted(lev.items())))
    pred={}
    for k in itertools.product((0,np.pi),repeat=3):
        w=sum(1 for x in k if x==0); m=np.sqrt(w*(tf+tc)**2+(3-w)*(tf-tc)**2); pred.setdefault(w,[]).append(round(m,6))
    print("   formule m_w = sqrt(w (t_f+t_c)^2 + (3-w)(t_f-t_c)^2), w = nombre de directions a k_d = 0 :", {w:(len(v),v[0]) for w,v in sorted(pred.items())})
# etiquettes du niveau parent au point uniforme : translations d'un octet et moments
H,S,ix=Hbloc(1.0,1.0)
T=[]
for d in range(3):
    M=np.zeros((64,64))
    for s in S:
        t=list(s); t[d]=(t[d]+2)%4; M[ix[tuple(t)],ix[s]]=1   # deux pas : les signes echelonnes sont invariants
    T.append(M)
print("translations d'un octet : commutent a H :",[bool(np.allclose(M@H,H@M)) for M in T])
E,V=np.linalg.eigh(H); out={}
for i in range(64):
    v=V[:,i]; lab=tuple(int(round(v@(M@v))) for M in T)   # valeurs propres +-1 = e^{i k_d}, si l'etat est propre
    pass
# projecteurs par etiquette
for lab in itertools.product((1,-1),repeat=3):
    P=np.eye(64)
    for d in range(3): P=P@(np.eye(64)+lab[d]*T[d])/2
    Hl=P@H@P; e=np.linalg.eigvalsh(Hl); e=np.round(np.abs(e[np.abs(e)>1e-9]),4); w=sum(1 for x in lab if x==1)
    out.setdefault(w,[]).append((lab,sorted(set(e)) if len(e) else [0.0], int(round(np.trace(P)))))
for w in sorted(out):
    print(f"   w = {w} : {len(out[w])} moment(s), |E| = {out[w][0][1]}, dimension par moment {out[w][0][2]}")
# echange k -> k + (pi,pi,pi) : multiplication par (-1)^(n1+n2+n3) sur les positions d'octet
Sh=np.diag([(-1.0)**(s[0]//2+s[1]//2+s[2]//2) for s in S])
print("echange k -> k+(pi,pi,pi) : envoie H(t_f,t_c) sur H(t_f,-t_c) :",bool(np.allclose(Sh@H@Sh,Hbloc(1.0,-1.0)[0])),"| commute a H :",bool(np.allclose(Sh@H,H@Sh)))
print("conclusion : un seul moment sans masse (w = 0) ; triplet w = 1 a |E| = 2 t ; triplet w = 2 a |E| = 2 sqrt(2) t ; pas d'appariement.")
