from s4enum import *
import json
print({k:len(v) for k,v in SUB.items()})
PAIR=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
def pimg(g,k):
    img=frozenset(frozenset(g[x] for x in b) for b in PAIR[k])
    return [frozenset(frozenset(b) for b in Q) for Q in PAIR].index(img)
names=['x0','x1','x2','x3','p0','p1','p2','z','f','q','w']; ix={v:i for i,v in enumerate(names)}
perms=[]
for g in S4:
    perms.append([g[i] for i in range(4)]+[4+pimg(g,k) for k in range(3)]+[7,8,9,10])
d=json.load(open('/sessions/'+__import__('os').environ.get('X','')+'x' if False else __import__('os').path.expanduser('~/mnt/FHT_26/03_code/CalculCodes/articleV_FHT21/V25_K3_pieces.json')))
def cx(F):
    S=set()
    for f in F:
        f=tuple(sorted(ix[v] for v in f))
        for k in range(1,len(f)+1): S.update(itertools.combinations(f,k))
    return sorted(S)
A=cx(d['A']); B=cx([b for b in d['B'] if not any(v.startswith('y') for v in b)])
print('A0', homology_chars(A,perms)); print('B0', homology_chars(B,perms))
Z=sorted(set(A)&set(B)); print('Z', homology_chars(Z,perms), 'sommets Z', sorted({v for s in Z for v in s}))
