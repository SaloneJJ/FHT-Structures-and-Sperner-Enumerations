# Verification (Fable) : cloture des deux points ouverts (bimodule de l'octet ; generations = sommets)
import itertools, numpy as np
sites=list(itertools.product((0,1),repeat=3)); idx={s:i for i,s in enumerate(sites)}
def hop(d,ph):
    M=np.zeros((8,8))
    for s in sites:
        t=list(s); t[d]^=1; M[idx[tuple(t)],idx[s]]=ph(s)
    return M
a=[hop(0,lambda s:1),hop(1,lambda s:(-1)**s[0]),hop(2,lambda s:(-1)**(s[0]+s[1]))]
eps=np.diag([(-1.0)**sum(s) for s in sites]); P=eps@a[0]@a[1]@a[2]
Pi=[np.diag([(-1.0)**s[d] for s in sites]) for d in range(3)]
def alg_dim(gens):
    basis=[np.eye(8,dtype=complex)]; r=1
    while True:
        new=basis+[b@g for b in basis for g in gens]
        V=np.array([np.concatenate([b.real.ravel(),b.imag.ravel()]) for b in new])
        u,s,vt=np.linalg.svd(V,full_matrices=False); k=int((s>1e-9).sum())
        basis=[(vt[i,:64]+1j*vt[i,64:]).reshape(8,8) for i in range(k)]
        if k==r: return k
        r=k
def real_structure(gens):
    # J = U kappa commutant aux generateurs : U conj(G) = G U ; renvoie la dimension des solutions et le signe de J^2
    rows=[np.kron(np.eye(8),G)-np.kron(G.conj().T,np.eye(8)).T.conj().T*0 for G in gens]  # placeholder
    A=[]
    for G in gens:
        # vec(G U - U conj(G)) = (I (x) G - conj(G)^T (x) I) vec(U)   (vec par colonnes)
        A.append(np.kron(np.eye(8),G)-np.kron(G.conj().T,np.eye(8)))
    u,s,vt=np.linalg.svd(np.vstack(A)); null=vt[(s>1e-9).sum():].conj()
    if len(null)==0: return 0,None
    U=null[0].reshape(8,8,order='F'); J2=U@U.conj(); c=J2[0,0]/abs(J2[0,0]) if abs(J2[0,0])>1e-9 else np.trace(J2)/8
    return len(null), round((np.trace(J2)/np.trace(U@U.conj().T)).real,6)
five=[x.astype(complex) for x in a+[eps,P]]
six=[x.astype(complex) for x in a]+[1j*a[d]@Pi[d] for d in range(3)]
ok6=all(np.allclose(X,X.conj().T) and np.allclose(X@X,np.eye(8)) for X in six) and all(np.allclose(six[i]@six[j]+six[j]@six[i],0) for i in range(6) for j in range(i))
print("six involutions hermitiennes anticommutantes :",ok6)
for name,g in (("Gamma (six)",six),("i Gamma (six)",[1j*x for x in six]),("Gamma (cinq)",five),("i Gamma (cinq)",[1j*x for x in five])):
    n,j2=real_structure(g); print(f"  {name:15s}: dim_R = {alg_dim(g):3d} ; structures antilineaires commutantes : {n} ; J^2 = {j2}")
# releves projectifs des 24 rotations du cube sur l'octet a flux pi
H=a[0]+a[1]+a[2]; lifts=[]
for perm in itertools.permutations(range(3)):
    for sg in itertools.product((1,-1),repeat=3):
        R=np.zeros((3,3)); 
        for i in range(3): R[i,perm[i]]=sg[i]
        if round(np.linalg.det(R))!=1: continue
        Pg=np.zeros((8,8))
        for s in sites:
            v=R@(2*np.array(s)-1); Pg[idx[tuple(int((x+1)//2) for x in v)],idx[s]]=1
        found=None
        for signs in itertools.product((1,-1),repeat=7):
            D=np.diag((1,)+signs); W=D@Pg
            if np.allclose(W@H@W.T,H): found=W; break
        lifts.append(found)
print("releves trouves :",sum(x is not None for x in lifts),"/ 24")
L=[x.astype(complex) for x in lifts]
print("  dim_R algebre des rotations :",alg_dim(L)," ; avec les cinq Gamma :",alg_dim(L+five)," ; commutent-elles ? ",all(np.allclose(x@g,g@x) for x in L for g in five))
# S4 permutant quatre generateurs : caractere du relevement d'un 3-cycle sur S (dim 4)
t=lambda i,j:(five[i]-five[j])/np.sqrt(2)
c3=t(0,1)@t(1,2); print("  trace du relevement d'un 3-cycle sur C^8 :",round(np.trace(c3).real,6),"-> sur S :",round(np.trace(c3).real/2,6),"(4 irreductible : +1 ; 2+2' : -2)")
# generations : invariants du stabilisateur d'un sommet dans H_k(sk_k Delta^{n-1})
def subsets(n,m): return list(itertools.combinations(range(n),m))
def bmat(n,m):
    rows=subsets(n,m-1); cols=subsets(n,m); ri={s:i for i,s in enumerate(rows)}; B=np.zeros((len(rows),len(cols)))
    for j,s in enumerate(cols):
        for q in range(m): B[ri[s[:q]+s[q+1:]],j]=(-1)**q
    return B
def psign(seq):
    s=1
    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            if seq[i]>seq[j]: s=-s
    return s
print("invariants de Stab(g_0) dans H_k(sk_k), par n (lignes) et k = 0..n-2 :")
for n in (4,5,6,7):
    out=[]
    for k in range(0,n-1):
        cols=subsets(n,k+1); ci={s:i for i,s in enumerate(cols)}
        u,s,vt=np.linalg.svd(bmat(n,k+1)); Z=vt[(s>1e-9).sum():].T if k>=0 else None
        Pz=Z@Z.T; avg=np.zeros((len(cols),len(cols))); cnt=0
        for p in itertools.permutations(range(1,n)):
            g=(0,)+p; M=np.zeros((len(cols),len(cols)))
            for j,sx in enumerate(cols):
                im=[g[v] for v in sx]; M[ci[tuple(sorted(im))],j]=psign(im)
            avg+=M; cnt+=1
            if cnt>=720: break
        avg/=cnt; out.append(int(round(np.trace(Pz@avg))))
    print("  n =",n,":",out)
