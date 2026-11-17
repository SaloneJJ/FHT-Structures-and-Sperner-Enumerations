# M3 (V145 §2) : graphe du roi sur le tore (L = 2^(n-1)) ; valeurs propres exactes lambda(k) = 27 - prod_i (1 + 2 cos k_i).
# Meme procedure que V145 (pente de log N(lambda) sur le quantile 8 %, et P(t) sur t in [0.12, 0.6]), puis tailles croissantes.
import numpy as np
for L in (16, 32, 64, 128):
    k = 2 * np.pi * np.arange(L) / L; c = 1 + 2 * np.cos(k)
    ev = np.sort((27 - np.einsum('i,j,k->ijk', c, c, c)).ravel())
    lam = ev[1:]; N = np.arange(1, len(lam) + 1)
    sel = (lam > 0) & (lam < np.quantile(lam, 0.08))
    p = np.polyfit(np.log(lam[sel]), np.log(N[sel]), 1)[0]
    # fenetre "petits lambda" a fraction fixe du spectre mais loin de lambda_1 : quantile 0.2 % a 2 %
    sel2 = (lam > np.quantile(lam, 0.002)) & (lam < np.quantile(lam, 0.02))
    p2 = np.polyfit(np.log(lam[sel2]), np.log(N[sel2]), 1)[0] if sel2.sum() > 10 else np.nan
    ts = np.logspace(np.log10(0.12), np.log10(0.6), 10)
    P = np.array([np.exp(-ev * t).mean() - 1 / ev.size for t in ts])
    s = np.polyfit(np.log(ts), np.log(P), 1)[0]
    # fenetre adaptee a la taille : 1 << D t << L^2/(4 pi^2), D = 9
    ts2 = np.logspace(np.log10(1.0), np.log10(L**2 / (4 * np.pi**2 * 9) / 4), 10)
    P2 = np.array([np.exp(-ev * t).mean() - 1 / ev.size for t in ts2])
    s2 = np.polyfit(np.log(ts2), np.log(P2), 1)[0]
    lam1 = ev[ev > 1e-12][0]
    print(f"L={L:4d} : lambda_1 = {lam1:.6f} ; d_s comptage (8 %) = {2*p:.3f} ; comptage (0,2-2 %) = {2*p2:.3f} ; "
          f"retour t in [0.12,0.6] : {-2*s:.3f} ; retour fenetre adaptee [{ts2[0]:.1f},{ts2[-1]:.1f}] : {-2*s2:.3f}")
