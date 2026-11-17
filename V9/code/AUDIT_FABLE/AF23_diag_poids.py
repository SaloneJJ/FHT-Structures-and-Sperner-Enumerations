import numpy as np, AF23_dirac_natif_S3 as A
for N in (2,):
    H, cells = A.build(N,'S3'); E,U = np.linalg.eigh(H.toarray())
    o=np.argsort(np.abs(E))[:16]; lev=np.array([c[0] for c in cells]); fac=[c[1] for c in cells]
    # distance de la cellule aux aretes de la racine : coordonnee au bord de la facette
    for k in o:
        w=np.abs(U[:,k])**2; wc=w.reshape(-1,8).sum(1)
        wl=[wc[lev==n].sum() for n in range(N+1)]
        fin=(lev==N); m=2**N
        edge=np.array([any(x in (0,m-1) for x in c[2]) for c in cells])
        print(f'E={E[k]:+.5f}  poids par niveau {np.round(wl,3)}  poids niveau fin pres des aretes de la racine {wc[fin&edge].sum()/max(wc[fin].sum(),1e-12):.3f}')
    print('fraction des cellules fines pres des aretes :', np.mean([any(x in (0,2**N-1) for x in c[2]) for c in cells if c[0]==N]))
