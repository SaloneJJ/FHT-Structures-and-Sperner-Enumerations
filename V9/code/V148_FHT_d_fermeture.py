"""V148 (d) : D1 spectre de D_Sigma pour les 8 structures spin eps et action de l'inflation ; D2 coloriage L/R (bipartisme) du tore 2^n ;
(e, complement) commutateur SLAC N impair, algebre tronquee.   JJ SALONE, FHT 2026"""
import itertools, numpy as np
print('=== D1 : spectre du Dirac continu -i sigma.grad sur T^3 (periode 1) pour eps in (Z/2)^3 : valeurs propres +-2 pi |m + eps/2| ===')
for eps in itertools.product((0,1),repeat=3):
    best=min(np.linalg.norm(np.array(m)+np.array(eps)/2) for m in itertools.product(range(-3,4),repeat=3))
    mult=sum(1 for m in itertools.product(range(-3,4),repeat=3) if abs(np.linalg.norm(np.array(m)+np.array(eps)/2)-best)<1e-12)
    print(f'  eps = {eps} : plus petite |valeur propre| = 2 pi x {best:.4f} = {2*np.pi*best:.4f} ; multiplicite des vecteurs d\'onde = {mult} (x2 signes x2 spin)')
print('  inflation E : x -> 2x mod 1. Mode exp(2 i pi k.x) o E = exp(2 i pi (2k).x) ; 2k = 2m + eps entier => periodique : E* envoie toute structure eps sur 000.')
print('  verification : 2(m + eps/2) entier pour tous m, eps :', all(float(2*(m+e/2))==int(2*(m+e/2)) for m in range(-5,6) for e in (0,1)))
print('\n=== D2 : coloriage 2 couleurs L/R des cellules du tore (Z/2^n)^3, adjacence par faces, echange L<->R a chaque pas ===')
def bip(n):
    s=2**n; col={}; ok=True
    start=(0,0,0); col[start]=0; stack=[start]
    while stack:
        c=stack.pop()
        for ax in range(3):
            for d in (1,-1):
                v=list(c); v[ax]=(v[ax]+d)%s; v=tuple(v)
                if v==c: ok=False; continue          # cellule collee a elle-meme (n = 0)
                if v not in col: col[v]=1-col[c]; stack.append(v)
                elif col[v]==col[c]: ok=False
    return ok
for n in range(0,5): print(f'  n = {n} (cote {2**n}) : coloriage L/R coherent (holonomie triviale) : {bip(n)}')
print('\n=== (e, complement) SLAC N impair : [D, P f P] (Toeplitz) et [D, f] avec enroulement, 1D, D = N p, p = 2 pi m / N ===')
def comm(N,fh,trunc):
    m=np.arange(-(N//2),N-(N//2)); p=2*np.pi*m/N; d=N*p
    T=np.zeros((N,N),complex)
    for i,mi in enumerate(m):
        for j,mj in enumerate(m):
            if trunc: T[i,j]=fh.get(mi-mj,0.0)
            else: T[i,j]=sum(fh.get(mi-mj+k*N,0.0) for k in (-1,0,1))
    return np.linalg.norm(np.diag(d)@T-T@np.diag(d),2)
cos={1:0.5,-1:0.5}
for trunc in (True,False):
    print('  cos 2 pi x,',('algebre tronquee' if trunc else 'multiplication avec enroulement'),': N = 17,33,65,129,257 ->',' '.join(f'{comm(N,cos,trunc):8.3f}' for N in (17,33,65,129,257)),'(2 pi = 6.283)')
