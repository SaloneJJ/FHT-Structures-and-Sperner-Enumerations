# Verification (Fable) : pieces A et S de A2 -- signe des sommets, caracteres de l'homologie reduite
import itertools, numpy as np
def closure(maxfaces):
    F=set()
    for m in maxfaces:
        m=tuple(sorted(m))
        for r in range(0,len(m)+1):
            for s in itertools.combinations(m,r): F.add(s)
    return F
def psign(seq):
    seq=list(seq); s=1
    for i in range(len(seq)):
        for j in range(i+1,len(seq)):
            if seq[i]>seq[j]: s=-s
    return s
def perm_sign(p):  # p: dict vertex->vertex
    keys=sorted(p); return psign([p[k] for k in keys])
def homology_chars(F, perms):
    byd={}
    for s in F: byd.setdefault(len(s)-1,[]).append(s)
    for d in byd: byd[d].sort()
    idx={d:{s:i for i,s in enumerate(byd[d])} for d in byd}
    top=max(byd)
    def bd(d):  # C_d -> C_{d-1}
        M=np.zeros((len(byd[d-1]),len(byd[d])))
        for j,s in enumerate(byd[d]):
            for t in range(len(s)): M[idx[d-1][s[:t]+s[t+1:]],j]=(-1)**t
        return M
    B={d:bd(d) for d in range(0,top+1)}
    def act(g,d):
        M=np.zeros((len(byd[d]),len(byd[d])))
        for j,s in enumerate(byd[d]):
            im=[g[v] for v in s]; M[idx[d][tuple(sorted(im))],j]=psign(im)
        return M
    def proj(M):
        if M.size==0: return np.zeros((M.shape[0],M.shape[0]))
        u,sv,vt=np.linalg.svd(M,full_matrices=False); r=(sv>1e-9).sum(); U=u[:,:r]; return U@U.T
    res={}
    for d in range(0,top+1):
        Pb=proj(B[d+1]) if d+1<=top else np.zeros((len(byd[d]),)*2)   # bords dans C_d
        Pc=proj(B[d].T)                                               # (ker d)^perp
        row=[]
        for g in perms:
            G=act(g,d); row.append(round(np.trace(G)-np.trace(G@Pb)-np.trace(G@Pc)))
        res[d]=row
    return res
# ---- piece A : Z={z,x0..x3}, X3={p0,p1,p2}, f, q
x=[0,1,2,3]; z=4; pair={5:((0,1),(2,3)),6:((0,2),(1,3)),7:((0,3),(1,2))}; f=8; q=9
maxA=[(z,0,1,2,3),(5,6,7),(q,)]+[(f,a,b) for a,b in itertools.combinations((5,6,7),2)]
for p,(ab,cd) in pair.items(): maxA+=[(p,z)+ab,(p,z)+cd]
FA=closure(maxA)
def liftA(s):
    g={z:z,f:f,q:q}
    for i in x: g[i]=s[i]
    for p,(ab,cd) in pair.items():
        im=tuple(sorted((s[ab[0]],s[ab[1]])))
        for p2,(ab2,cd2) in pair.items():
            if im in (ab2,cd2): g[p]=p2
    return g
reps={'e':(0,1,2,3),'(01)':(1,0,2,3),'(01)(23)':(1,0,3,2),'(012)':(1,2,0,3),'(0123)':(1,2,3,0)}
print("A : faces max", len(maxA), "sommets 10")
print("  sgn_V sur les 24 elements :", sorted(set(perm_sign(liftA(s)) for s in itertools.permutations(x))))
hc=homology_chars(FA,[liftA(s) for s in reps.values()])
print("  caracteres de H~_d(A) sur e,(01),(01)(23),(012),(0123) :", {d:hc[d] for d in hc if any(hc[d])})
print("  attendus : 1=[1,1,1,1,1]  2=[2,0,2,-1,0]  sgn=[1,-1,1,1,-1]")
# ---- piece S : F_3^2
vecs=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
def line(v): return frozenset({v,((2*v[0])%3,(2*v[1])%3)})
lines=sorted({line(v) for v in vecs}, key=sorted)
aff=[]
for a,b in vecs:  # droite ax+by=1
    aff.append(frozenset((u,w) for u in range(3) for w in range(3) if (a*u+b*w)%3==1))
names={}; 
for v in vecs: names[('v',v)]=len(names)
for l in lines: names[('l',l)]=len(names)
names[('c',)]=len(names)
for t in aff: names[('t',t)]=len(names)
c=names[('c',)]
maxS=[(names[('v',v)],names[('l',line(v))]) for v in vecs]+[(names[('v',v)],c) for v in vecs]+[(names[('t',t)],names[('v',v)],c) for t in aff for v in t]
FS=closure(maxS)
GL=[m for m in itertools.product(range(3),repeat=4) if (m[0]*m[3]-m[1]*m[2])%3!=0]
def liftS(m):
    ap=lambda v:((m[0]*v[0]+m[1]*v[1])%3,(m[2]*v[0]+m[3]*v[1])%3)
    g={c:c}
    for v in vecs: g[names[('v',v)]]=names[('v',ap(v))]
    for l in lines: g[names[('l',l)]]=names[('l',frozenset(ap(v) for v in l))]
    for t in aff: g[names[('t',t)]]=names[('t',frozenset(ap(v) for v in t))]
    return g
ok=all(perm_sign(liftS(m))==(1 if (m[0]*m[3]-m[1]*m[2])%3==1 else -1) for m in GL)
print("S : sommets", len(names), "| |GL(2,3)| =", len(GL), "| sgn_V = det sur les 48 elements :", ok)
hs=homology_chars(FS,[liftS((1,0,0,1)),liftS((2,0,0,2)),liftS((1,0,0,2))])
print("  caracteres de H~_d(S) sur I, -I, diag(1,-1) :", {d:hs[d] for d in hs if any(hs[d])}, " attendu H~_1 : [4,-4,0]")
