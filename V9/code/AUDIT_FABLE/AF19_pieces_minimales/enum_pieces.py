import sys, time, itertools
from s4enum import *
REPS = [S4.index(g) for g in [(0,1,2,3),(1,0,2,3),(1,0,3,2),(1,2,0,3),(1,2,3,0)]]
TYPES = {1:['S4'],2:['A4'],3:['D4'],4:['S3',"V4'"+'_?'][:1]+['C4'][:0],}
ORB = [('S4',1),('A4',2),('D4',3),('S3',4),('C4',6),("V4'",6),('C3',8)]
# NB : S4/C4 et S4/V4' ont 6 points ; V4' = sous-groupe de Klein non normal
def sign_sort(im):
    arr=sorted(range(len(im)),key=lambda a:im[a]); sg=1; arr=list(arr)
    for a in range(len(arr)):
        while arr[a]!=a: b=arr[a]; arr[a],arr[b]=arr[b],arr[a]; sg=-sg
    return sg
def lefschetz_orbit(orb, perms):
    """contribution de l'orbite (liste de simplexes) a Lambda(g) non reduit, pour g dans REPS"""
    out=[]
    for gi in REPS:
        g=perms[gi]; tot=0
        for s in orb:
            im=[g[v] for v in s]
            if set(im)==set(s): tot+=(-1)**(len(s)-1)*sign_sort(im)
        out.append(tot)
    return out
TARGET = {'B0':[-3,-1,1,0,1],'B1':[5,3,1,2,1],'A0':[1,1,1,4,1],'A1':[1,1,1,-2,1]}
WANT = {'B0':{(1,'1'):1,(1,'3'):1}, 'B1':{(0,'1'):1,(0,'3'):1},
        'A0':{(0,'1'):1,(1,'2'):1,(2,'s'):1}, 'A1':{(1,'1'):1,(2,'2'):1,(3,'s'):1}}
def search(types, which, limit_time=150):
    perms,lab=build_W(types); n=len(lab)
    orbs=[[tuple(sorted(x)) for x in o] for o in simplex_orbits(perms,n)]
    orbs.sort(key=lambda o:(len(o[0]),o[0]))
    key={frozenset(s):i for i,o in enumerate(orbs) for s in o}
    faces=[]
    for o in orbs:
        s=o[0]; faces.append({key[frozenset(s[:i]+s[i+1:])] for i in range(len(s))} if len(s)>1 else set())
    L=[lefschetz_orbit(o,perms) for o in orbs]
    nv=sum(1 for o in orbs if len(o[0])==1)
    found=[]; cnt=[0]; t0=time.time()
    tg=[TARGET[w] for w in which]
    def dfs(i, inc, lam):
        if time.time()-t0>limit_time: raise TimeoutError
        if i==len(orbs):
            cnt[0]+=1
            for w,t in zip(which,tg):
                if lam==t:
                    S=[x for j in inc for x in orbs[j]]
                    h=homology_chars(S,perms)
                    if h==WANT[w]: found.append((w,[orbs[j][0] for j in inc]))
            return
        o=orbs[i]
        if len(o[0])==1:   # sommets : tous presents
            dfs(i+1, inc+[i], [a+b for a,b in zip(lam,L[i])]); return
        dfs(i+1, inc, lam)
        if faces[i] <= set(inc):
            dfs(i+1, inc+[i], [a+b for a,b in zip(lam,L[i])])
    try: dfs(0,[], [0]*5); status='complet'
    except TimeoutError: status='INCOMPLET (temps)'
    return status, cnt[0], len(orbs), found
if __name__=='__main__':
    which=sys.argv[1].split(',')
    for spec in sys.argv[2:]:
        types=spec.split('+')
        t=time.time(); st,c,no,f=search(types, which)
        print(f'{spec:20s} orbites {no:4d} complexes {c:9d} {st} ; trouves {[(w, len(s)) for w,s in f][:5]} ({time.time()-t:.1f} s)', flush=True)
