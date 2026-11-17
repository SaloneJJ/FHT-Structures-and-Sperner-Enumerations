import numpy as np, itertools, AF23_dirac_natif_S3 as A
facets=[(i,s) for i in range(4) for s in (0,1)]; tang={f:[a for a in range(4) if a!=f[0]] for f in facets}
flip={}
for f in facets:
    i,s=f; tA=tang[f]
    for d in range(3):
        for e in (1,-1):
            j=tA[d]; g=(j,1 if e>0 else 0); tB=tang[g]
            sj=e*(1 if s==0 else -1); cmap=[tB.index(i) if a==j else tB.index(a) for a in tA]; sig=[sj if a==j else 1 for a in tA]
            V=A.transport(None,cmap,sig); flip[(f,g)]=np.allclose(V@A.FG@V.conj().T,-A.FG)
# 2-coloration : signe eps_f tel que eps_g = -eps_f si flip, +eps_f sinon
eps={facets[0]:1}; stack=[facets[0]]; ok=True
while stack:
    f=stack.pop()
    for (a,b),fl in flip.items():
        if a!=f: continue
        want=-eps[f] if fl else eps[f]
        if b in eps:
            if eps[b]!=want: ok=False
        else: eps[b]=want; stack.append(b)
print('signes de facette coherents :', ok, eps)
