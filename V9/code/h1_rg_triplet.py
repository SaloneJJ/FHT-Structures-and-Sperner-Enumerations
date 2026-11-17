"""V72 : RG a une boucle du MS (V = lambda |H|^4, m_H = v sqrt(2 lambda)), condition au bord spectrale
lambda(Lambda) = 4 g^2 (n^2+3)/(n+3)^2 (CCM 2012, g = g2(Lambda)), seuil du triplet en M_Delta : lambda -> lambda - dl.
dl pour m_H = 125 GeV ; r = M_R/M_Delta et m_II/m_I (V71) avec c3 omega k = xi = 2 (hypothese O(1)).
Simplifications : pas de y_nu ni du triplet dans la RG entre Lambda et M_Delta ; y_t mesure. JJ SALONE, FHT 2026"""
import numpy as np
v, mt = 246.22, 173.0
def beta(y):
    g1, g2, g3, yt, l = y; k = 1 / (16 * np.pi**2)
    return k * np.array([41/6 * g1**3, -19/6 * g2**3, -7 * g3**3,
        yt * (4.5 * yt**2 - 17/12 * g1**2 - 9/4 * g2**2 - 8 * g3**2),
        24 * l**2 + 12 * l * yt**2 - 6 * yt**4 - 3 * l * (3 * g2**2 + g1**2) + 3/8 * (2 * g2**4 + (g2**2 + g1**2)**2)])
def run(y, t0, t1, n=1500):
    h = (t1 - t0) / n; y = np.array(y, float)
    for _ in range(n):
        k1 = beta(y); k2 = beta(y + h/2*k1); k3 = beta(y + h/2*k2); k4 = beta(y + h*k3)
        y = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    return y
y0 = [0.3587, 0.6483, 1.1666, 0.9369, 0.126]
print('controle : lambda(mt) = 0.126 -> m_H =', round(v*np.sqrt(2*0.126), 1))
for Lam in (1e16, 1e17):
    tL = np.log(Lam / mt); yL = run(y0, 0, tL); g2L, ytL = yL[1], yL[3]
    print(f'Lambda = {Lam:.0e} : g1 {yL[0]:.3f} g2 {g2L:.3f} g3 {yL[2]:.3f} y_t {ytL:.3f} lambda_MS(Lambda) {yL[4]:.4f}')
    for n in (0.1, 0.25, 0.5, 1.0, 2.0):
        lL = 4 * g2L**2 * (n**2 + 3) / (n + 3)**2
        for MD in (1e13, 1e14):
            def mH(dl):
                y = yL.copy(); y[4] = lL
                y = run(y, tL, np.log(MD / mt)); y[4] -= dl
                y = run(y, np.log(MD / mt), 0); return v * np.sqrt(2 * max(y[4], 1e-9))
            m0 = mH(0); lo, hi = 0.0, lL
            for _ in range(40):
                mid = (lo + hi) / 2
                (lo, hi) = (mid, hi) if mH(mid) > 125 else (lo, mid)
            dl = (lo + hi) / 2; xi = 2.0
            r = np.sqrt(dl * (n + 3)**2 / (4 * g2L**2 * xi * n**2)); ratio = dl * (n + 3) / (xi * n**2 * ytL**2)
            print(f'  n={n:4.2f} lam(L)={lL:.3f} M_D={MD:.0e} : m_H sans seuil {m0:5.1f} ; dl(125)={dl:.3f} ; r~{r:.2f} ; m_II/m_I~{ratio:.1f}')
