from twisted_alg import *
import numpy as np, itertools
# GL(2,3) sur les 8 vecteurs non nuls de F_3^2
vecs=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
def mat_perm(M): return tuple(vecs.index(((M[0][0]*a+M[0][1]*b)%3,(M[1][0]*a+M[1][1]*b)%3)) for a,b in vecs)
GL=closure([mat_perm([[1,1],[0,1]]),mat_perm([[0,1],[2,0]]),mat_perm([[2,0],[0,1]])]); print('GL(2,3)',len(GL))
def det_char(p):
    # determinant via image de (1,0),(0,1)
    a,c=vecs[p[vecs.index((1,0))]]; b,d=vecs[p[vecs.index((0,1))]]
    return 1 if (a*d-b*c)%3==1 else -1
decompose(GL, lambda p:1, 'GL(2,3), trivial'); decompose(GL, det_char, 'GL(2,3), det')
# 2O : quaternions unitaires (groupe octaedral binaire), action reguliere
def qmul(a,b):
    a0,a1,a2,a3=a; b0,b1,b2,b3=b
    return (a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2, a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0)
s=1/np.sqrt(2)
gens=[(s,s,0,0),(0.5,0.5,0.5,0.5)]
E=[(1.0,0.0,0.0,0.0)]; fr=list(E)
key=lambda q: tuple(round(x,6)+0.0 for x in q)
seen={key(E[0])}
while fr:
    x=fr.pop()
    for g in gens:
        y=qmul(x,g)
        if key(y) not in seen: seen.add(key(y)); E.append(y); fr.append(y)
print('2O', len(E)); K=[key(q) for q in E]; ix={k:i for i,k in enumerate(K)}
perms=[tuple(ix[key(qmul(g,x))] for x in E) for g in E]
O=closure(perms[:]) ; print(len(O))
# caractere lineaire non trivial de 2O : signe de l'image dans S4 (parite) = action sur les 6 axes... on utilise sign de la permutation reguliere ? (non) -> on prend le caractere via l'abelianise : C2
# image dans S4 : action par conjugaison sur les 4 paires d'axes {+-(1+-i+-j+-k)/2}... on se contente du caractere : g dans 2T (24 elements du sous-groupe d'indice 2) -> +1
T=set(k for k in K if all(abs(abs(x)-0.5)<1e-6 for x in k) or sum(abs(x)>1e-6 for x in k)==1)
def chi2T(p):
    g=E[[i for i in range(len(E)) if perms[i]==p][0]] if False else None
    return 1
idxp={p:i for i,p in enumerate(perms)}
eps2O=lambda p: 1 if K[idxp[p]] in T else -1
print('|2T|',len(T))
decompose(perms and sorted(perms), lambda p:1, '2O, trivial'); decompose(sorted(perms), eps2O, '2O, eps (2O/2T)')
