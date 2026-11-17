import numpy as np, itertools, AF23_dirac_natif_S3 as A
facets=[(i,s) for i in range(4) for s in (0,1)]; tang={f:[a for a in range(4) if a!=f[0]] for f in facets}
res={'FG conserve':0,'FG change de signe':0,'autre':0}; G5c={'G5 conserve':0,'autre':0}
for f in facets:
    i,s=f; tA=tang[f]
    for d in range(3):
        for e in (1,-1):
            j=tA[d]; g=(j,1 if e>0 else 0); tB=tang[g]
            sj=e*(1 if s==0 else -1); cmap=[tB.index(i) if a==j else tB.index(a) for a in tA]; sig=[sj if a==j else 1 for a in tA]
            V=A.transport(None,cmap,sig); F2=V@A.FG@V.conj().T
            if np.allclose(F2,A.FG): res['FG conserve']+=1
            elif np.allclose(F2,-A.FG): res['FG change de signe']+=1
            else: res['autre']+=1
            G2=V@A.G5@V.conj().T; G5c['G5 conserve' if np.allclose(G2,A.G5) else 'autre']+=1
print('passages de facette (8 facettes x 6) :', res, G5c)
