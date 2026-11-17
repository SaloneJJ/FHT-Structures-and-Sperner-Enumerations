import sys, time, itertools
from s4enum import *
from enum_pieces import lefschetz_orbit
SZ={'S4':1,'A4':2,'D4':3,'S3':4,'C4':6,"V4'":6,'C3':8}
def lam_target(want):
    # Lambda non reduit = 1 + sum_k (-1)^k tr(H~_k)
    out=[]
    for c in range(5):
        v=1
        for (k,r),m in want.items(): v+=(-1)**k*m*CHI[r][c]
        out.append(v)
    return out
def all_complexes(types, tl=60):
    perms,lab=build_W(types); n=len(lab)
    orbs=[[tuple(sorted(x)) for x in o] for o in simplex_orbits(perms,n)]
    orbs.sort(key=lambda o:(len(o[0]),o[0]))
    key={frozenset(s):i for i,o in enumerate(orbs) for s in o}
    faces=[{key[frozenset(o[0][:i]+o[0][i+1:])] for i in range(len(o[0]))} if len(o[0])>1 else set() for o in orbs]
    L=[lefschetz_orbit(o,perms) for o in orbs]
    res=[]; t0=time.time()
    def dfs(i,inc,lam):
        if time.time()-t0>tl: raise TimeoutError
        if i==len(orbs): res.append((tuple(lam),inc)); return
        if len(orbs[i][0])==1: dfs(i+1,inc+[i],[a+b for a,b in zip(lam,L[i])]); return
        dfs(i+1,inc,lam)
        if faces[i]<=set(inc): dfs(i+1,inc+[i],[a+b for a,b in zip(lam,L[i])])
    try: dfs(0,[],[0]*5); st='complet'
    except TimeoutError: st='incomplet'
    return perms, orbs, res, st
REPS3={'3':{'3':1},"3'":{"3'":1},'1+2':{'1':1,'2':1},'s+2':{'s':1,'2':1},'1+1+1':{'1':3},'s+s+s':{'s':3},'1+1+s':{'1':2,'s':1},'s+s+1':{'s':2,'1':1}}
TENS_S={'3':"3'","3'":'3','1+2':'s+2','s+2':'1+2','1+1+1':'s+s+s','s+s+s':'1+1+1','1+1+s':'s+s+1','s+s+1':'1+1+s'}
def want(rep, d): return {(d,r):m for r,m in REPS3[rep].items()}
def multisets(maxn):
    order=list(SZ); out=[]
    def rec(i,cur,tot):
        if i==len(order):
            if cur: out.append(cur); 
            return
        k=0
        while tot+k*SZ[order[i]]<=maxn: rec(i+1,cur+[order[i]]*k,tot+k*SZ[order[i]]); k+=1
    rec(0,[],0); return sorted(out,key=lambda w:sum(SZ[t] for t in w))
if __name__=='__main__':
    maxn=int(sys.argv[1]); t0=time.time()
    best={}
    for W in multisets(maxn):
        n=sum(SZ[t] for t in W)
        perms,orbs,res,st=all_complexes(W, tl=25)
        lamset={}
        for lam,inc in res: lamset.setdefault(lam,[]).append(inc)
        for rep in REPS3:
            for d in (0,1,2):
                key=(rep,d)
                if key in best: continue
                tC=tuple(lam_target(want(rep,d))); tCp=tuple(lam_target(want(TENS_S[rep],d+1)))
                if tC not in lamset or tCp not in lamset: continue
                okC=[inc for inc in lamset[tC] if homology_chars([x for j in inc for x in orbs[j]],perms)==want(rep,d)]
                okCp=[inc for inc in lamset[tCp] if homology_chars([x for j in inc for x in orbs[j]],perms)==want(TENS_S[rep],d+1)]
                if okC and okCp:
                    best[key]=(n,'+'.join(W),st, [orbs[j][0] for j in okC[0]], [orbs[j][0] for j in okCp[0]])
        if time.time()-t0>160: print('arret temps', n); break
    for k,v in sorted(best.items(), key=lambda kv:kv[1][0]):
        print(k, 'sommets', v[0], v[1], v[2], '\n    C  :', v[3], '\n    C\' :', v[4])
