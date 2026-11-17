import sys, os; _H = os.path.dirname(os.path.abspath(__file__)); [sys.path.insert(0, os.path.join(_H, "..", p)) for p in ("AF23_dirac_S3", "V352", "AF24_holonomies")]; sys.path.insert(0, _H)
import numpy as np, itertools
exec(open(os.path.join(_H,'cov.py')).read().split("N=1\n")[0])
N=1; H,cells=Hor(N,0.0); Hd=H.toarray(); fac=[(i,s) for i in range(4) for s in (0,1)]
lev=np.repeat([c[0] for c in cells],8); cf=np.repeat([fac.index(c[1]) for c in cells],8)
P,chi=sym(N,cells,[0,2,3,1],set())
D=np.abs(P@Hd@P.conj().T-Hd)
same=cf[:,None]==cf[None,:]; samelev=lev[:,None]==lev[None,:]; diag=np.eye(len(Hd),dtype=bool)
blk=np.kron(np.eye(len(cells)),np.ones((8,8))).astype(bool)
print('sur site       :',D[blk].max())
print('meme facette, meme niveau, hors site :',D[same&samelev&~blk].max())
print('meme facette, niveaux differents (vertical) :',D[same&~samelev].max())
print('facettes differentes :',D[~same].max())
# test d'une seule facette : la rotation 1->2->3->1 envoie (0,0) sur elle-meme
