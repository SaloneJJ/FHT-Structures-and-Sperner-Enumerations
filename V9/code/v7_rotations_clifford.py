# Verification (Fable) : algebre engendree par les releves des 24 rotations et les cinq generateurs de l'octet ; type reel ou quaternionique
import itertools, numpy as np, os, sys
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'v5_closure.py')).read()
exec(src.split('print("six involutions')[0])
H=a[0]+a[1]+a[2]; lifts=[]
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
            if np.allclose(W@H@W.T,H): lifts.append(W.astype(complex)); break
G=lifts+five
print("dim_R de l'algebre engendree :",alg_dim(G))
A=[np.kron(np.eye(8),X)-np.kron(X.conj().T,np.eye(8)) for X in G]
u,s,vt=np.linalg.svd(np.vstack(A)); null=vt[(s>1e-9).sum():].conj()
U=null[0].reshape(8,8,order='F'); J2=U@U.conj()
print("structures antilineaires commutantes :",len(null),"; J^2 =",round((np.trace(J2)/np.trace(U@U.conj().T)).real,6),"(+1 : type reel, donc M8(R))")
