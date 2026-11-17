from twisted_alg import *
# B4 = C2 wr S4 : permutations signees de 4 objets, sur 8 points (i, i+4 = -i)
def sp(perm, signs):
    p=[0]*8
    for i in range(4):
        j=perm[i]; s=signs[i]
        p[i]= j if s==1 else j+4; p[i+4]= j+4 if s==1 else j
    return tuple(p)
B4=closure([sp((1,0,2,3),(1,1,1,1)), sp((1,2,3,0),(1,1,1,1)), sp((0,1,2,3),(-1,1,1,1))]); print(len(B4))
def sgnS4(p):  # signe de la permutation des 4 paires
    q=tuple(p[i]%4 for i in range(4)); return sign(q)
def nneg(p): return (-1)**sum(1 for i in range(4) if p[i]>=4)
for nm,eps in [('B4 trivial',lambda p:1),('B4 sgn(S4)',sgnS4),('B4 (-1)^neg',nneg),('B4 sgn*neg',lambda p:sgnS4(p)*nneg(p))]: decompose(B4,eps,nm)
