import sys, os; _H = os.path.dirname(os.path.abspath(__file__)); [sys.path.insert(0, os.path.join(_H, "..", p)) for p in ("AF23_dirac_S3", "V352", "AF24_holonomies")]; sys.path.insert(0, _H)
import numpy as np, itertools, sys
exec(open(os.path.join(_H,'cov.py')).read().split("N=1\n")[0])
import AF24_theoreme_holonomies as T
N=1; H,cells=Hor(N,0.0); Hd=H.toarray(); fac=[(i,s) for i in range(4) for s in (0,1)]; F=len(fac)
cfi=np.array([fac.index(c[1]) for c in cells]); cf=np.repeat(cfi,8)
pairs=sorted({(min(a,b),max(a,b)) for a in range(F) for b in range(F) if a!=b and fac[a][0]!=fac[b][0]})
pid={p:k for k,p in enumerate(pairs)}; print('paires adjacentes :',len(pairs))
def img(f,sigma,refl):
    i,s=fac[f]; i2=sigma[i]; return fac.index((i2,(1-s) if i2 in refl else s))
gens=[('cycle',[0,2,3,1],set()),('12r0',[0,2,1,3],{0}),('01r2',[1,0,2,3],{2}),('r0r1',[0,1,2,3],{0,1})]
rows=[]; rhs=[]; nT=len(pairs)
for gi,(name,sigma,refl) in enumerate(gens):
    P,chi=sym(N,cells,sigma,refl); R=P@Hd@P.conj().T
    for (a,b) in pairs:
        ga,gb=img(a,sigma,refl),img(b,sigma,refl)
        m=np.outer(cf==ga,cf==gb); h=Hd[m]; r=R[m]; nz=np.abs(h)>1e-9
        e=r[nz]/h[nz]; assert np.allclose(e,e[0]) and np.isclose(abs(e[0]),1), (name,a,b,e[:4])
        E=0 if e[0].real>0 else 1
        row=np.zeros(nT+F*len(gens),np.uint8); row[pid[(min(ga,gb),max(ga,gb))]]^=1; row[pid[(a,b)]]^=1
        row[nT+gi*F+ga]^=1; row[nT+gi*F+gb]^=1; rows.append(row); rhs.append(E)
M=np.array(rows); b=np.array(rhs,np.uint8)
x=T.solve2(M,b); print('choix de signes symetrique (4 generateurs) : existe =',x is not None)
if x is not None:
    tsol=x[:nT]; print('faces a changer de signe :',[ (fac[p[0]],fac[p[1]]) for p,k in pid.items() if tsol[k]])
    np.save('tsol.npy',tsol); import pickle; pickle.dump((pairs,fac),open('pairs.pkl','wb'))
print('rang',T.rank2(M),'/ inconnues',M.shape[1])
