"""V357 (relecture de V341) : chaîne réduite de profondeur (N = 60), μ = r = t_c = 1, m_K = √3(ε − 1), η = +1 ; secteur s = +1.
Modes nuls par SVD du bloc chiral (seuil 1e-6) : goût f, chiralité, niveau du maximum, décroissance mesurée vers l'intérieur."""
import numpy as np
from V357_chaine_profondeur import chain
N = 60
print('ε      m_K     λ_surv  λ_jum   modes de bord (niveau ≥ 50) : goût, chiralité, λ mesuré | |E| min par goût')
for eps in (1.0, 0.9, 0.8, 0.6, 0.5, 0.45, 0.43, 0.42, 0.40, 0.35, 0.3, 0.2, 0.1, 0.0, -0.3, -0.6, -0.75):
    mK = np.sqrt(3) * (eps - 1); out = []; emin = {}
    for f in (+1, -1):
        H = chain(N, mK, 1.0, np.ones(N + 1), +1, f)
        sz = np.tile([1., -1.], N + 1); ip, im = np.where(sz > 0)[0], np.where(sz < 0)[0]
        Q = H[np.ix_(ip, im)]; U, sv, Vh = np.linalg.svd(Q)
        emin[f] = sv.min()
        for j in np.where(sv < 1e-6)[0]:
            for c, p in ((+1, np.abs(U[:, j]) ** 2), (-1, np.abs(Vh[j]) ** 2)):
                n0 = int(np.argmax(p))
                if n0 >= 50:
                    lam = np.sqrt(p[N - 3] / p[N - 2]) if p[N - 2] > 1e-200 else 0.0
                    out.append(f'f={f:+d} c={c:+d} λ={lam:.3f}')
    print(f'{eps:+.2f}  {mK:+.3f}  {abs(mK):.3f}   {abs(2 + mK):.3f}   ' + ('; '.join(out) if out else '— aucun —') + f'  | {emin[1]:.1e} {emin[-1]:.1e}')
