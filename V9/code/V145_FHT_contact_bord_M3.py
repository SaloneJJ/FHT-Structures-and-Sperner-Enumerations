"""V145 (complement M3) : dimension de marche d_w et dimension spectrale d_s du graphe de contact du tore T3
construit combinatoirement (meme construction que V145_FHT_contact_bord.py).
d_w par la loi d'echelle de la premiere valeur propre non nulle (lambda_1 ~ diam^-d_w), d_s = 2 d_H / d_w,
et controle par la probabilite de retour P(t) dans une fenetre intermediaire. JJ SALONE, FHT 2026"""
import numpy as np, sys, io, contextlib
from collections import defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    import V145_FHT_contact_bord as B
lam1 = {}; Pret = {}
for n in (2, 3, 4):
    cells, uf, glued = B.construire(n, True)[n]
    pos = {w: i for i, w in enumerate(cells)}
    byc = defaultdict(list)
    for w in cells:
        for x in B.C: byc[uf.f((w, x))].append(pos[w])
    A = np.zeros((len(cells),) * 2)
    for lst in byc.values():
        for a in lst:
            for b in lst:
                if a != b: A[a, b] = 1
    L = np.diag(A.sum(1)) - A
    ev = np.sort(np.linalg.eigvalsh(L))
    lam1[n] = ev[ev > 1e-9][0]
    Pret[n] = ev
    print(f'n={n} : N={len(cells)}, diam={2**(n-1)}, lambda_1={lam1[n]:.6f}', flush=True)
for n in (3, 4):
    dw = np.log(lam1[n-1] / lam1[n]) / np.log(2)
    print(f'd_w (n={n-1}->{n}) = {dw:.4f} ; d_s = 2*3/d_w = {6/dw:.4f}')
ev = Pret[4]; N = len(ev)
ts = np.logspace(np.log10(0.12), np.log10(0.6), 10)  # fenetre 1/D < t < L^2/(4 pi^2 D), D ~ 8.9 (lambda_1 = D (2 pi/16)^2)
P = np.array([np.exp(-ev * t).sum() / N for t in ts])
# soustraire la contribution du mode zero (1/N) pour la fenetre intermediaire
Pc = P - 1.0 / N
sl = np.polyfit(np.log(ts), np.log(Pc), 1)[0]
print(f'retour P(t) (n=4, t de {ts[0]:.2f} a {ts[-1]:.2f}, fenetre etroite, indicatif) : pente {sl:.3f}, d_s = {-2*sl:.3f}')
