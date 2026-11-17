"""Assemblage par secteurs leptons | quarks : X = U_S (A_S^L * B_S^R * G0) U (B_S^L * A_S^R * K4), S in {lep, quark}.
Recherche constructive sur un cote W donne ; filtre Betti Z/2 = (12,36,36,12) en degres 3..6 ; sortie des candidats pour AF9."""
import sys, time, itertools, json, os
sys.path.insert(0, os.path.expanduser('~/mnt/FHT_26/03_code/CalculCodes/articleV_FHT21'))
from fht_tools import betti_z2, maximal
from s4enum import *
from enum_pieces import lefschetz_orbit
from gen_pairs import lam_target
def enum_target(perms, Vsub, want, tl=60):
    Vsub=sorted(Vsub); n=len(Vsub); orbs=[]; seen=set()
    for k in range(1,n+1):
        for s in itertools.combinations(Vsub,k):
            fs=frozenset(s)
            if fs in seen: continue
            o={frozenset(p[i] for i in s) for p in perms}
            if not all(set(x)<=set(Vsub) for x in o): continue
            seen|=o; orbs.append(sorted(tuple(sorted(x)) for x in o))
    key={frozenset(s):i for i,o in enumerate(orbs) for s in o}
    faces=[{key[frozenset(o[0][:i]+o[0][i+1:])] for i in range(len(o[0]))} if len(o[0])>1 else set() for o in orbs]
    L=[lefschetz_orbit(o,perms) for o in orbs]; tg=lam_target(want); out=[]; t0=time.time()
    def dfs(i,inc,lam):
        if time.time()-t0>tl: raise TimeoutError
        if i==len(orbs):
            if lam==tg:
                S=[x for j in inc for x in orbs[j]]
                if homology_chars(S,perms)==want: out.append(frozenset(S))
            return
        if len(orbs[i][0])==1: dfs(i+1,inc+[i],[a+b for a,b in zip(lam,L[i])]); return
        dfs(i+1,inc,lam)
        if faces[i]<=set(inc): dfs(i+1,inc+[i],[a+b for a,b in zip(lam,L[i])])
    try: dfs(0,[],[0]*5); st='complet'
    except TimeoutError: st='incomplet'
    return out, st
def facets(S): return [s for s in S if not any(set(s)<set(t) for t in S)]
def assemble(names, sectors):
    G0=['g0','g1','g2','g3']; K4=[list(e) for e in itertools.combinations(G0,2)]
    F=[]
    for A,B in sectors:
        fa=[[ 'L_'+names[v] for v in s] for s in facets(A)]; fbR=[[ 'R_'+names[v] for v in s] for s in facets(B)]
        fbL=[[ 'L_'+names[v] for v in s] for s in facets(B)]; faR=[[ 'R_'+names[v] for v in s] for s in facets(A)]
        F+=[a+b+[g] for a in fa for b in fbR for g in G0]
        F+=[b+a+e for b in fbL for a in faR for e in K4]
    return maximal([frozenset(f) for f in F])
if __name__=='__main__':
    types=['S3','D4']+['S4']*int(sys.argv[1]); perms,lab=build_W(types)
    names=['x0','x1','x2','x3','p0','p1','p2']+[f'f{i}' for i in range(int(sys.argv[1]))]
    fixed=list(range(7,len(names)))
    t=time.time()
    As=[]; 
    for sub in itertools.combinations(fixed, 3):
        a,st=enum_target(perms,[4,5,6]+list(sub),{(0,'1'):1,(1,'2'):1,(2,'s'):1}); As+=a
    print('A (sans x) :',len(As), round(time.time()-t), flush=True)
    Bq=[]
    for k in range(0,len(fixed)+1):
        for sub in itertools.combinations(fixed,k):
            b,st=enum_target(perms,[0,1,2,3]+list(sub),{(1,'3'):1}, tl=30); Bq+=b
    print('B couleur :',len(Bq), flush=True)
    Bl=[]
    for k in range(2,len(fixed)+1):
        for sub in itertools.combinations(fixed,k):
            b,st=enum_target(perms,list(sub),{(1,'1'):1}); Bl+=b
    print('B lepton :',len(Bl), round(time.time()-t), flush=True)
    target={3:12,4:36,5:36,6:12}; hits=[]; n=0; t0=time.time()
    for Al,Aq in itertools.product(As,As):
        if Al==Aq: continue
        for bl in Bl:
            for bq in Bq:
                n+=1
                if n<=int(sys.argv[3]): continue
                X=assemble(names,[(Al,bl),(Aq,bq)])
                verts=set().union(*X)
                if betti_z2([list(f) for f in X])==target: hits.append((Al,Aq,bl,bq,len(verts)))
                if time.time()-t0>float(sys.argv[2]): break
            if time.time()-t0>float(sys.argv[2]): break
        if time.time()-t0>float(sys.argv[2]): print('arret temps'); break
    print('combinaisons', n, 'depuis', sys.argv[3], 'Betti conformes', len(hits))
    json.dump([[sorted(map(list,h[0])),sorted(map(list,h[1])),sorted(map(list,h[2])),sorted(map(list,h[3])),h[4]] for h in hits[:20]], open('secteurs2_hits.json','w'))
    for h in hits[:5]: print('sommets', h[4])
