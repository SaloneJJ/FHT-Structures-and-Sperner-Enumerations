import itertools, numpy as np, sys, time
S4 = list(itertools.permutations(range(4)))
def mul(a,b): return tuple(a[b[i]] for i in range(4))
def cls(g):
    seen=set(); L=[]
    for i in range(4):
        if i in seen: continue
        j=i; c=0
        while j not in seen: seen.add(j); j=g[j]; c+=1
        L.append(c)
    return tuple(sorted(L))
CL = {(1,1,1,1):0,(1,1,2):1,(2,2):2,(1,3):3,(4,):4}
CHI = {'1':[1,1,1,1,1],'s':[1,-1,1,1,-1],'2':[2,0,2,-1,0],'3':[3,1,-1,0,-1],"3'":[3,-1,-1,0,1]}
SIZE=[1,6,3,8,6]
# sous-groupes (representants des classes) pour les orbites transitives S4/H
def closure(gens):
    G={tuple(range(4))}; fr=list(G)
    while fr:
        x=fr.pop(); 
        for g in gens:
            y=mul(x,g)
            if y not in G: G.add(y); fr.append(y)
    return G
SUB = {'S4':closure([(1,0,2,3),(1,2,3,0)]), 'A4':closure([(1,2,0,3),(1,0,3,2)]), 'D4':closure([(1,2,3,0),(1,0,2,3)][:1]+[(2,1,0,3)]),
       'S3':closure([(1,0,2,3),(1,2,0,3)]), 'C4':closure([(1,2,3,0)]), "V4'":closure([(1,0,2,3),(0,1,3,2)]), 'C3':closure([(1,2,0,3)]),
       'C2t':closure([(1,0,2,3)]), 'C2d':closure([(1,0,3,2)]), '1':{tuple(range(4))}}
def cosets(H):
    reps=[]; seen=set()
    for g in S4:
        c=frozenset(mul(g,h) for h in H)
        if c not in seen: seen.add(c); reps.append(c)
    return reps
def orbit_action(H):
    C=cosets(H); idx={c:i for i,c in enumerate(C)}
    return [[idx[frozenset(mul(g,x) for x in c)] for c in C] for g in S4]
def build_W(types):
    acts=[orbit_action(SUB[t]) for t in types]; off=0; perms=[[] for _ in S4]; labels=[]
    for t,a in zip(types,acts):
        n=len(a[0])
        for gi in range(24): perms[gi]+= [off+x for x in a[gi]]
        labels += [f'{t}#{i}' for i in range(n)]; off+=n
    return perms, labels
def simplex_orbits(perms, n):
    seen=set(); orbs=[]
    for k in range(1,n+1):
        for s in itertools.combinations(range(n),k):
            fs=frozenset(s)
            if fs in seen: continue
            o={frozenset(p[i] for i in s) for p in perms}; seen|=o; orbs.append(sorted(o, key=sorted))
    return orbs
def homology_chars(simplices, perms):
    """dims des composantes isotypiques de H~_k (reduite) ; simplices : liste de tuples tries"""
    byd={-1:[()]}
    for s in simplices: byd.setdefault(len(s)-1,[]).append(s)
    for k in byd: byd[k].sort()
    idx={k:{s:i for i,s in enumerate(v)} for k,v in byd.items()}
    maxd=max(byd)
    # matrices d'action et de bord
    def act(k,g):
        M=np.zeros((len(byd[k]),len(byd[k])))
        for j,s in enumerate(byd[k]):
            im=[g[v] for v in s]; order=sorted(range(len(im)),key=lambda a:im[a]); 
            # signe de la permutation qui trie
            sg=1; arr=list(order); 
            for a in range(len(arr)):
                while arr[a]!=a: b=arr[a]; arr[a],arr[b]=arr[b],arr[a]; sg=-sg
            M[idx[k][tuple(sorted(im))],j]=sg
        return M
    def bd(k):
        M=np.zeros((len(byd[k-1]),len(byd[k])))
        for j,s in enumerate(byd[k]):
            for i in range(len(s)): M[idx[k-1][s[:i]+s[i+1:]],j]=(-1)**i
        return M
    D={k:bd(k) for k in range(0,maxd+1)}
    A={k:[act(k,g) for g in perms] for k in range(-1,maxd+1)}
    res={}
    for r,chi in CHI.items():
        E={k:(chi[0]/24)*sum(chi[CL[cls(S4[gi])]]*A[k][gi] for gi in range(24)) for k in range(-1,maxd+1)}
        dimE={k:np.linalg.matrix_rank(E[k]) for k in E}
        rk={}
        for k in range(0,maxd+1):
            rk[k]=np.linalg.matrix_rank(D[k]@E[k]) if E[k].size else 0
        for k in range(-1,maxd+1):
            h=dimE[k]-rk.get(k,0)-rk.get(k+1,0)
            if h: res[(k,r)]=int(round(h/ (CHI[r][0])))
    return res
