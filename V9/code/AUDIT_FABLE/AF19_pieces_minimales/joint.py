import sys, time, itertools
from s4enum import *
from enum_pieces import lefschetz_orbit, TARGET, WANT
def enum_on(perms, Vsub, which, tl=120):
    Vsub=sorted(Vsub); n=len(Vsub)
    orbs=[]; seen=set()
    for k in range(1,n+1):
        for s in itertools.combinations(Vsub,k):
            fs=frozenset(s)
            if fs in seen: continue
            o={frozenset(p[i] for i in s) for p in perms}
            if not all(set(x)<=set(Vsub) for x in o): continue
            seen|=o; orbs.append(sorted(tuple(sorted(x)) for x in o))
    key={frozenset(s):i for i,o in enumerate(orbs) for s in o}
    faces=[{key[frozenset(o[0][:i]+o[0][i+1:])] for i in range(len(o[0]))} if len(o[0])>1 else set() for o in orbs]
    L=[lefschetz_orbit(o,perms) for o in orbs]; tg=[TARGET[w] for w in which]
    out=[]; t0=time.time(); cnt=[0]
    def dfs(i,inc,lam):
        if time.time()-t0>tl: raise TimeoutError
        if i==len(orbs):
            cnt[0]+=1
            for w,t in zip(which,tg):
                if lam==t:
                    S=[x for j in inc for x in orbs[j]]
                    if homology_chars(S,perms)==WANT[w]: out.append((w,frozenset(S)))
            return
        if len(orbs[i][0])==1: dfs(i+1,inc+[i],[a+b for a,b in zip(lam,L[i])]); return
        dfs(i+1,inc,lam)
        if faces[i]<=set(inc): dfs(i+1,inc+[i],[a+b for a,b in zip(lam,L[i])])
    try: dfs(0,[],[0]*5); st='complet'
    except TimeoutError: st='incomplet'
    return out, st, cnt[0]
if __name__=='__main__':
    types=sys.argv[1].split('+'); perms,lab=build_W(types); print(lab)
    VA=[int(x) for x in sys.argv[2].split(',')]; VB=[int(x) for x in sys.argv[3].split(',')]
    t=time.time(); As,sa,ca=enum_on(perms,VA,['A0']); print('A', len(As), sa, ca, round(time.time()-t,1), flush=True)
    t=time.time(); Bs,sb,cb=enum_on(perms,VB,['B0']); print('B', len(Bs), sb, cb, round(time.time()-t,1), flush=True)
    good=[]
    for _,A in As:
        for _,B in Bs:
            Z=A&B
            if not Z: continue
            if homology_chars(sorted(Z),perms)=={}: good.append((A,B))
    print('paires (A,B) avec A inter B acyclique :', len(good))
    for A,B in good[:3]:
        fa=[s for s in A if not any(set(s)<set(t) for t in A)]; fb=[s for s in B if not any(set(s)<set(t) for t in B)]
        print(' A facettes', sorted(fa,key=len), '\n B facettes', sorted(fb,key=len), '\n Z', sorted(A&B,key=len)[-3:])
