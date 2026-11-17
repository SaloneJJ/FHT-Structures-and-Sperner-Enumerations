"""V148 (e), complements apres lecture de e1/e2 : (1) Schur raffine (commutation avec les 6 generateurs, support) ;
(2) SLAC : N pair/impair, option Nyquist, comptage des noeuds ; (3) commutateurs avec algebre tronquee (Toeplitz)."""
import itertools, numpy as np, sys
sys.argv=[sys.argv[0],'none']
exec(open('V148_FHT_e_raccord.py').read().split("for p, fn in (('Ea'")[0])   # reutilise le modele F (meme graine)
print('=== (1) Schur raffine ===')
gens=[]
for _ in range(4):
    a=alea_alg(); Ja=Japp(UF,a)
    gens += [a+a.conj().T, 1j*(a-a.conj().T), Ja+Ja.conj().T, 1j*(Ja-Ja.conj().T)]
w1,V1=np.linalg.eigh(gens[0]*1.0+0.37*gens[2]+0.11*gens[4])
groups=[];cur=[0]
for i in range(1,dF):
    if abs(w1[i]-w1[cur[0]])<1e-8: cur.append(i)
    else: groups.append(cur);cur=[i]
groups.append(cur)
basis=[]
for g in groups:
    for a_ in g:
        for b_ in g:
            E=np.zeros((dF,dF),complex);E[a_,b_]=1;basis.append(V1@E@V1.conj().T)
def kern(odd):
    ncol=2*len(basis); G=np.zeros((ncol,ncol))
    cl=[lambda B:B-B.conj().T, lambda B:Japp(UF,B)-B]+([lambda B:B@GF+GF@B] if odd else [])+[(lambda B,X=X:B@X-X@B) for X in gens[1:]]
    for f in cl:
        cols=[]
        for Bm in basis:
            for c in (1.0,1.0j):
                v=f(c*Bm).ravel(); cols.append(np.concatenate([v.real,v.imag]))
        M=np.array(cols).T; G+=M.T@M
    ev,evec=np.linalg.eigh(G); return evec[:,ev<1e-9*ev.max()].T
ns=kern(True); ns0=kern(False)
print('  noyau avec les 16 generateurs : impair ->',len(ns),'; temoin sans la condition impair ->',len(ns0))
def toB(v): return sum(((v[2*i]+1j*v[2*i+1])*basis[i]) for i in range(len(basis)))
supp=np.zeros(dF,bool); Bs=[toB(v) for v in ns]
for B in Bs: supp|=(np.abs(B).sum(axis=0)>1e-8)|(np.abs(B).sum(axis=1)>1e-8)
sec=sorted(set(((i>=dI),(i%dI)//(nB*nG),((i%dI)//nG)%nB) for i in np.where(supp)[0]))
print('  etats touches :',int(supp.sum()),'sur',dF,'; secteurs (antiparticule?, a, b) :',sec)
print('\n=== (2) SLAC : noeuds de masse nulle selon N et option Nyquist ===')
def nodes(N,opt):
    m=np.arange(-(N//2),N-(N//2)); p=2*np.pi*m/N
    def sd(x):
        if N%2==0 and opt=='zero': x=np.where(np.isclose(np.abs(x),np.pi),0.0,x)
        if N%2==0 and opt=='pi':   x=np.where(np.isclose(x,-np.pi),np.pi,x)    # branche asymetrique -pi -> +pi
        return x
    s=sd(p); # noeuds : valeurs de m ou s=0 (en 1D), puis produit en 3D
    n1=int(np.sum(np.abs(s)<1e-12)); 
    # KO-1 : s(-k) = -s(k) ?
    idx=lambda mm:(mm+N//2)%N
    mm=m; ko=max(abs(sd(p)[idx(-mm)] + sd(p)[idx(mm)]) for _ in [0])
    return n1**3, float(np.max(np.abs(sd(2*np.pi*(-m)/N)+sd(p))))
for N in (6,8,9,16,17):
    for opt in (('zero','pi') if N%2==0 else ('-',)):
        n3,kerr=nodes(N,opt)
        print(f'  N = {N:2d} Nyquist->{opt:4s} : noeuds s=0 en 3D = {n3} ; ecart KO-1 (s(-k)+s(k)) = {kerr:.2e}')
print('  (comparaison : naif N pair -> 8 noeuds, N impair -> 8 aussi : sin p = 0 en p = 0 et p = pi, pi absent si N impair => 1 noeud)')
for N in (8,9):
    m=np.arange(-(N//2),N-(N//2)); p=2*np.pi*m/N; print(f'  naif N = {N} : noeuds 3D = {int(np.sum(np.abs(np.sin(p))<1e-12))**3}')

print('\n=== (3) commutateur [D, P f P] avec algebre tronquee (Toeplitz), SLAC, 1D ===')
def comm_trunc(N,f_hat):
    m=np.arange(-(N//2),N-(N//2)); p=2*np.pi*m/N
    d=N*np.where(np.isclose(np.abs(p),np.pi),0.0,p)    # D_phys diagonal en Fourier
    T=np.zeros((N,N),complex)
    for i,mi in enumerate(m):
        for j,mj in enumerate(m):
            T[i,j]=f_hat.get(mi-mj,0.0)
    return np.linalg.norm(np.diag(d)@T-T@np.diag(d),2)
cos_hat={1:0.5,-1:0.5}
tri_hat={0:0.5}; 
for k in range(1,400,2): tri_hat[k]=-2/(np.pi**2*k**2)/1.0*1; tri_hat[-k]=tri_hat[k]    # |2x-1| = 1/2 + sum_{k impair}  (4/(pi^2k^2)) cos(2 pi k x)/... normalisation ci-dessous
tri_hat={0:0.5}
for k in range(1,400,2): tri_hat[k]=2/(np.pi**2*k**2); tri_hat[-k]=tri_hat[k]
print('  N      : 16 32 64 128 256 ;  (derivee physique : D = N p, p = 2 pi m/N)')
for nom,fh in (('cos',cos_hat),('triangle ~|2x-1| (serie tronquee)',tri_hat)):
    print('  ',nom,' : ',' '.join(f'{comm_trunc(N,fh):8.3f}' for N in (16,32,64,128,256)), '  (naif tronque aussi :', ' '.join(f'{(lambda N:np.linalg.norm(np.diag(N*np.sin(2*np.pi*np.arange(-(N//2),N-(N//2))/N))@np.array([[fh.get(a-b,0.0) for b in range(-(N//2),N-(N//2))] for a in range(-(N//2),N-(N//2))])-np.array([[fh.get(a-b,0.0) for b in range(-(N//2),N-(N//2))] for a in range(-(N//2),N-(N//2))])@np.diag(N*np.sin(2*np.pi*np.arange(-(N//2),N-(N//2))/N)),2))(N):8.3f}' for N in (16,64,256)),')')
