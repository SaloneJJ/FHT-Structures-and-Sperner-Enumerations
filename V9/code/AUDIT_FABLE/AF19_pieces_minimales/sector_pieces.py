import sys, time, itertools
from s4enum import *
from gen_pairs import all_complexes, lam_target, multisets, SZ
ROWS=[('1',0),('2',1),('s',2)]          # H~(A0) : singulet deg 0, doublet deg 1, sgn deg 2
COLS=[('1',1),('3',1)]                   # H~(B0) : lepton, couleur en deg 1
targets={}
for k in range(1,4):
    for sub in itertools.combinations(range(3),k):
        for sh in (-1,0,1):
            w={}
            ok=True
            for i in sub:
                r,d=ROWS[i]; d+=sh
                if d<0: ok=False
                w[(d,r)]=w.get((d,r),0)+1
            if ok: targets[('A',tuple(ROWS[i][0] for i in sub),sh)]=w
for k in range(1,3):
    for sub in itertools.combinations(range(2),k):
        for sh in (-1,0,1):
            w={}
            for i in sub:
                r,d=COLS[i]; w[(d+sh,r)]=1
            targets[('B',tuple(COLS[i][0] for i in sub),sh)]=w
best={}; t0=time.time()
for W in multisets(int(sys.argv[1])):
    n=sum(SZ[t] for t in W)
    perms,orbs,res,st=all_complexes(W, tl=8)
    lamset={}
    for lam,inc in res: lamset.setdefault(lam,[]).append(inc)
    for key,w in targets.items():
        if key in best: continue
        tl=tuple(lam_target(w))
        if tl not in lamset: continue
        for inc in lamset[tl]:
            if homology_chars([x for j in inc for x in orbs[j]],perms)==w:
                best[key]=(n,'+'.join(W),st); break
    if time.time()-t0>140: print('arret temps a n =',n); break
for key in sorted(targets, key=lambda k:(k[0],len(k[1]),k[1],k[2])):
    print(key, '->', best.get(key,'non trouve (<= taille max)'))
