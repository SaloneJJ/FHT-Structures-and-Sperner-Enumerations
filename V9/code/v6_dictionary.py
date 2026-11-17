# Verification (Fable) : dictionnaire entre les six generateurs de l'octet et les pieces de K26 ; test du pont G0 = diagonales
import itertools, numpy as np
sites=list(itertools.product((0,1),repeat=3)); idx={s:i for i,s in enumerate(sites)}
def hop(d,ph):
    M=np.zeros((8,8))
    for s in sites:
        t=list(s); t[d]^=1; M[idx[tuple(t)],idx[s]]=ph(s)
    return M
c=[hop(0,lambda s:1),hop(1,lambda s:(-1)**s[0]),hop(2,lambda s:(-1)**(s[0]+s[1]))]
Pi=[np.diag([(-1.0)**s[d] for s in sites]) for d in range(3)]
ch=[c[d]@Pi[d] for d in range(3)]                       # c chapeau (antisymetriques, carre -1)
G5=(c[0]+c[1]+c[2])/np.sqrt(3); A=(c[0]-c[1])/np.sqrt(2); B=(c[0]+c[1]-2*c[2])/np.sqrt(6)
al=[1j*x for x in ch]; G4=G5@al[0]@al[1]@al[2]; F=1j*A@B
six=[G5.astype(complex),A.astype(complex),B.astype(complex)]+al
acomm=lambda X,Y:np.allclose(X@Y+Y@X,0); comm=lambda X,Y:np.allclose(X@Y-Y@X,0)
print("six generateurs {G5,A,B,alpha_d} : involutions hermitiennes anticommutantes :", all(np.allclose(X@X,np.eye(8)) and np.allclose(X,X.conj().T) for X in six) and all(acomm(six[i],six[j]) for i in range(6) for j in range(i)))
print("F = iAB : F^2=1", np.allclose(F@F,np.eye(8)), "| hermitien", np.allclose(F,F.conj().T), "| imaginaire pur", np.allclose(F.real,0),
      "| F = i(c1c2+c2c3+c3c1)/sqrt3 :", np.allclose(F,1j*(c[0]@c[1]+c[1]@c[2]+c[2]@c[0])/np.sqrt(3)),
      "| F G5 = i c1c2c3 :", np.allclose(F@G5,1j*c[0]@c[1]@c[2]))
print("F commute a G5, alpha_d, G4 :", all(comm(F,X) for X in [G5]+al+[G4]), "| anticommute a A, B :", acomm(F,A) and acomm(F,B))
# releves des rotations
H=c[0]+c[1]+c[2]; lifts=[]
for perm in itertools.permutations(range(3)):
    for sg in itertools.product((1,-1),repeat=3):
        R=np.zeros((3,3))
        for i in range(3): R[i,perm[i]]=sg[i]
        if round(np.linalg.det(R))!=1: continue
        Pg=np.zeros((8,8))
        for s in sites:
            v=R@(2*np.array(s)-1); Pg[idx[tuple(int((x+1)//2) for x in v)],idx[s]]=1
        for signs in itertools.product((1,-1),repeat=7):
            W=np.diag((1,)+signs)@Pg
            if np.allclose(W@H@W.T,H): lifts.append((R,W)); break
def cls(R):
    t=round(np.trace(R)); 
    return {3:'e',1:'C4',0:'C3',-1:('C2' if abs(np.diag(R)).sum()==3 else "C2'")}[t]
reps={}
for R,W in lifts: reps.setdefault(cls(R),W)
order=['e',"C2'",'C2','C3','C4']   # = e, transpositions, doubles transpositions, 3-cycles, 4-cycles de S4
tab={'1':[1,1,1,1,1],'s':[1,-1,1,1,-1],'2':[2,0,2,-1,0],'3':[3,1,-1,0,-1],"3'":[3,-1,-1,0,1]}
def char(basis):
    out=[]
    for k in order:
        W=reps[k]; out.append(int(round(sum(np.trace(b.conj().T@W@b@W.T)/np.trace(b.conj().T@b) for b in basis).real)))
    return out
def name(ch_):
    # decomposition en irreductibles de S4
    sizes=[1,6,3,8,6]; res=[]
    for n_,t in tab.items():
        m=round(sum(s*a*b for s,a,b in zip(sizes,ch_,t))/24)
        if m: res.append(f"{m}x{n_}" if m>1 else n_)
    return " + ".join(res)
for nm,bs in (("G5",[G5]),("{A,B}",[A,B]),("{alpha_d}",al),("G4",[G4]),("F",[F]),("Cl(A,B)={1,A,B,AB}",[np.eye(8),A,B,A@B]),
              ("algebre de saveur {1,G4A,G4B,AB}",[np.eye(8),G4@A,G4@B,A@B]),
              ("Cl(A,B).span{alpha_d,G4} (16)",[a@v for a in (np.eye(8),A,B,A@B) for v in al+[G4]])):
    ch_=char(bs); print(f"  {nm:34s}: caractere {ch_}  =  {name(ch_)}")
print("  une generation de K26, (1+2+s)x(1+3) : ", [ (tab['1'][i]+tab['2'][i]+tab['s'][i])*(tab['1'][i]+tab['3'][i]) for i in range(5)])
V=np.array([ (a@v).ravel() for a in (np.eye(8),A,B,A@B) for v in al+[G4]]); print("  rang des 16 operateurs :", np.linalg.matrix_rank(V))
# test du pont : operateurs sur les grandes diagonales
w3=1j*c[0]@c[1]@c[2]; trip=[w3@Pi[1]@Pi[2], w3@Pi[2]@Pi[0], w3@Pi[0]@Pi[1]]
print("somme sur les diagonales = F G5 :", np.allclose(w3,F@G5), "| caractere du triplet :", name(char(trip)))
print("le triplet commute-t-il a G5 ?", [comm(t,G5) for t in trip], "| a alpha_d ?", [[comm(t,a) for a in al] for t in trip][0])
A_=[]
for X in [G5]+al: A_.append(np.kron(np.eye(8),X)-np.kron(X.T,np.eye(8)))
s_=np.linalg.svd(np.vstack(A_),compute_uv=False); print("dimension (complexe) du commutant de {G5, alpha_d} :", int((s_<1e-9).sum()))
