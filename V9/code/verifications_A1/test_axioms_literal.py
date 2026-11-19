"""Literal check of the axioms (A1)-(A5) (A5 component-wise) of Article I for v<=4: enumerate all families E of subsets of size >= 2 and all level functions E -> {1,...,v-1}; count the realizations up to isomorphism. Expected 1, 4, 23 for v = 2, 3, 4 (and 1 for v = 1). Self-contained version of the original V0_FHT_Test_axiomes.py."""
import sys, itertools
from itertools import product, combinations

class M:
    @staticmethod
    def popcount(x): return bin(x).count('1')
    @staticmethod
    def canon_py(v,E):
        best=None
        for p in itertools.permutations(range(v)):
            key=tuple(sorted(sum(1<<p[x] for x in range(v) if e>>x&1) for e in E))
            if best is None or key<best: best=key
        return best

def components(edges):
    left=list(edges); out=[]
    while left:
        reach=left[0]; grp=[left[0]]; rest=left[1:]; ch=True
        while ch:
            ch=False; keep=[]
            for e in rest:
                if e & reach: reach|=e; grp.append(e); ch=True
                else: keep.append(e)
            rest=keep
        out.append((reach,grp)); left=rest
    return out

def valid(v, E, lev):
    full=(1<<v)-1
    N=max(lev)
    comps={0:[(1<<x,[]) for x in range(v)]}
    for n in range(1,N+1):
        En=[e for e,l in zip(E,lev) if l==n]
        if not En: return False
        cs=components(En)
        for reach,grp in cs:
            # connected Sperner
            for a in grp:
                for b in grp:
                    if a!=b and (a&b)==a: return False
        comps[n]=cs
    # A3
    if len(comps[N])!=1 or comps[N][0][0]!=full: return False
    # A2 disjointness within level
    for n in range(1,N+1):
        vs=[r for r,_ in comps[n]]
        for i in range(len(vs)):
            for j in range(i+1,len(vs)):
                if vs[i]&vs[j]: return False
    # A4
    for n in range(1,N+1):
        for _,grp in comps[n]:
            for e in grp:
                for k in range(0,n):
                    for r,_ in comps[k]:
                        if e & r and not ((e & r)==r and e!=r): return False
    # A5 component-wise
    for n in range(1,N+1):
        for r,_ in comps[n]:
            if not any((r2 & r)==r2 and r2!=r for r2,_ in comps[n-1]): return False
    return True

for v in (2,3,4):
    subs=[sum(1<<x for x in c) for r in range(2,v+1) for c in combinations(range(v),r)]
    classes=set()
    for k in range(1,len(subs)+1):
        for E in combinations(subs,k):
            if M.popcount(0)>0: pass
            for lev in product(range(1,v),repeat=k):
                if valid(v,E,lev):
                    classes.add(M.canon_py(v,E)); break
    print("v=",v,"literal count =",len(classes))
