import sys, os; _H = os.path.dirname(os.path.abspath(__file__)); [sys.path.insert(0, os.path.join(_H, "..", p)) for p in ("AF23_dirac_S3", "V352", "AF24_holonomies")]; sys.path.insert(0, _H)
import numpy as np, itertools, sys
src=open(os.path.join(_H,'cov3.py')).read()
pre=src.split("rows=[]")[0]; exec(pre)
def system(gs):
    rows=[]; rhs=[]
    for gi,(name,sigma,refl) in enumerate(gs):
        P,chi=sym(N,cells,sigma,refl); R=P@Hd@P.conj().T
        for (a,b) in pairs:
            ga,gb=img(a,sigma,refl),img(b,sigma,refl)
            m=np.outer(cf==ga,cf==gb); h=Hd[m]; r=R[m]; nz=np.abs(h)>1e-9; e=r[nz]/h[nz]
            if not (np.allclose(e,e[0]) and np.isclose(abs(e[0]),1)): return 'non scalaire'
            row=np.zeros(nT+F*len(gs),np.uint8); row[pid[(min(ga,gb),max(ga,gb))]]^=1; row[pid[(a,b)]]^=1
            row[nT+gi*F+ga]^=1; row[nT+gi*F+gb]^=1; rows.append(row); rhs.append(0 if e[0].real>0 else 1)
    return T.solve2(np.array(rows),np.array(rhs,np.uint8)) is not None
nT=len(pairs)
allg=[('cycle 123',[0,2,3,1],set()),('12 & r0',[0,2,1,3],{0}),('01 & r2',[1,0,2,3],{2}),('r0 r1',[0,1,2,3],{0,1}),('r1 r2',[0,1,2,3],{1,2}),('cycle 012',[1,2,0,3],set()),('0123 cycle4 & r?',[1,2,3,0],{0})]
for g in allg: print(f'{g[0]:16s} seul :',system([g]))
for a,b in itertools.combinations(allg[:6],2): print(f'{a[0]} + {b[0]} :',system([a,b]))
