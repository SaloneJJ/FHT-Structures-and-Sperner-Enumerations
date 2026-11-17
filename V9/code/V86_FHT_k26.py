import json, itertools, time
from collections import defaultdict
def maximal(facets):
    fs = sorted(set(frozenset(f) for f in facets), key=len, reverse=True)
    out = []
    for f in fs:
        if not any(f < g for g in out): out.append(f)
    return out
def build_K():  # modele en joins V25 sans y ; JJ SALONE, FHT 2026
    d = json.load(open('V25_K3_pieces.json'))
    A = d["A"]; B = [b for b in d["B"] if not any(v.startswith('y') for v in b)]
    G4, E4 = d["generation"]["feuilles"], d["generation"]["K4"]
    I = [frozenset(['L_'+v for v in a]+['R_'+v for v in b]+[g]) for a in A for b in B for g in G4]
    II = [frozenset(['L_'+v for v in b]+['R_'+v for v in a]+e) for b in B for a in A for e in E4]
    return maximal(I+II), maximal(I), maximal(II)
def faces(F):
    s=set()
    for f in F:
        f=tuple(sorted(f))
        for k in range(1,len(f)+1):
            s.update(itertools.combinations(f,k))
    return s
def betti_from_faces(fs):
    byd=defaultdict(list)
    for s in fs: byd[len(s)-1].append(s)
    byd[-1]=[()]
    maxd=max(byd)
    ind={d:{s:i for i,s in enumerate(sorted(byd[d]))} for d in byd}
    def rank(d):
        piv={}; r=0
        for s in ind[d]:
            col=0
            for t in range(len(s)): col ^= 1 << ind[d-1][s[:t]+s[t+1:]]
            while col:
                h=col.bit_length()-1
                if h in piv: col^=piv[h]
                else: piv[h]=col; r+=1; break
        return r
    rk={d:rank(d) for d in range(0,maxd+1)}; rk[maxd+1]=0; rk[-1]=0
    return {d:len(ind[d])-rk[d]-rk[d+1] for d in range(-1,maxd+1) if len(ind[d])-rk[d]-rk[d+1]}
if __name__=='__main__':
    K,I,II=build_K()
    V=sorted(set().union(*K)); print(len(V),'feuilles',len(K),'facettes', max(map(len,K)))
    t=time.time(); fs=faces(K); print(len(fs),'faces', time.time()-t)
    t=time.time(); print(betti_from_faces(fs), time.time()-t)
