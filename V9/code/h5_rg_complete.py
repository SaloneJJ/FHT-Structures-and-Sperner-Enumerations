"""V75 (2) : RG a une boucle au-dessus de M_Delta (modele a triplet, equations de arXiv:1510.00799 App. A, g_B-L = 0, sans y_N),
puis MS sous M_Delta avec seuil dl = mu^2/M_Delta^2. Conditions au bord spectrales a Lambda (V73-V74) :
lam1 = 4g^2(n^2+3)/(n+3)^2, lam2 = 0, lam3 = 8g^2 Tr(mhat^4), lam4 = 0, lam5 = 16 g^2 n/(n+3), mu = 2 sqrt2 g n/(n+3) M_R,
y_Delta = kap * 2 sqrt2 g mhat (normalisation de y_Delta dans la reference non verifiee : kap in {1/sqrt2, 1}).
Simplifications : y_nu, y_e negliges dans la RG ; g = g2(Lambda). JJ SALONE, FHT 2026"""
import numpy as np
v, mt, P16 = 246.22, 173.0, 16 * np.pi**2
def beta_sm(y):
    g1, g2, g3, yt, l = y
    return np.array([41/6*g1**3, -19/6*g2**3, -7*g3**3, yt*(4.5*yt**2 - 17/12*g1**2 - 9/4*g2**2 - 8*g3**2),
        24*l**2 + 12*l*yt**2 - 6*yt**4 - 3*l*(3*g2**2 + g1**2) + 3/8*(2*g2**4 + (g2**2 + g1**2)**2)]) / P16
def beta_tr(y, ys):
    """y = g1 g2 g3 yt l1 l2 l3 l4 l5 lmu(=ln mu) ; ys = valeurs propres de y_Delta (diag, reelles)"""
    g1, g2, g3, yt, l1, l2, l3, l4, l5, lmu = y
    T2, T4 = np.sum(ys**2), np.sum(ys**4)
    b = np.zeros(10)
    b[0] = 47/6*g1**3; b[1] = -5/2*g2**3; b[2] = -7*g3**3
    b[3] = yt*(4.5*yt**2 - 8*g3**2 - 9/4*g2**2 - 17/12*g1**2)
    b[4] = 24*l1**2 + 3*l4**2 + 3*l4*l5 + 5/4*l5**2 + 3/8*(2*g2**4 + (g2**2 + g1**2)**2) - 3*l1*(3*g2**2 + g1**2) - 6*yt**4 + 12*l1*yt**2
    b[5] = 28*l2**2 + 24*l2*l3 + 6*l3**2 + 2*l4**2 + 2*l4*l5 + 15*g2**4 + 12*g1**4 - 12*g1**2*g2**2 + 4*l2*(-6*g2**2 - 3*g1**2 + 2*T2)
    b[6] = 18*l3**2 + 24*l2*l3 + l5**2 - 6*g2**4 - 12*g1**4 + 24*g1**2*g2**2 - 4*T4 + 4*l3*(-6*g2**2 - 3*g1**2 + 2*T2)
    b[7] = (4*l4**2 + 12*l1*l4 + 16*l2*l4 + 12*l3*l4 + 4*l1*l5 + 6*l2*l5 + 2*l3*l5 + l5**2 + 6*g2**4 + 3*g1**4 - 6*g1**2*g2**2
            + l4*(-33/2*g2**2 - 15/2*g1**2 + 6*yt**2 + 4*T2))
    b[8] = l5*(4*l1 + 4*l2 + 8*l3 + 8*l4 + 4*l5 - 33/2*g2**2 - 15/2*g1**2 + 6*yt**2 + 4*T2) + 12*g1**2*g2**2
    b[9] = 4*l1 + 4*l4 + 6*l5 - 21/2*g2**2 - 9/2*g1**2 + 2*T2 + 6*yt**2          # d ln mu / dt (de lambda0)
    return b / P16
def beta_ys(y, ys):
    g1, g2 = y[0], y[1]; T2 = np.sum(ys**2)
    return ys * (6*ys**2 - 3/2*(g1**2 + 3*g2**2) + 2*T2) / P16      # (1/2+1/2... ) termes y_l, y_D negliges ; 3 y^3 x 2 (y^T y + y y^T)
def rk(f, y, t0, t1, n):
    h = (t1 - t0) / n
    for _ in range(n):
        k1 = f(y); k2 = f(y + h/2*k1); k3 = f(y + h/2*k2); k4 = f(y + h*k3); y = y + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return y
y0 = np.array([0.3587, 0.6483, 1.1666, 0.9369, 0.126])
def solve(Lam, n, MD, mhat, kap):
    tL, tD = np.log(Lam/mt), np.log(MD/mt)
    yD_ = rk(beta_sm, y0, 0, tD, 800)                                       # jauge/yt jusqu'a M_D (MS)
    # au-dessus de M_D : on remonte jauge+yt avec le contenu triplet pour obtenir g(Lambda)
    up = lambda Y: np.concatenate([beta_tr(np.concatenate([Y[:4], [0, 0, 0, 0, 0, 0]]), np.zeros(3))[:4]])
    yL4 = rk(up, yD_[:4], tD, tL, 400); g = yL4[1]
    l1 = 4*g**2*(n**2 + 3)/(n + 3)**2; l3 = 8*g**2*np.sum(mhat**4); l5 = 16*g**2*n/(n + 3)
    ys0 = kap*2*np.sqrt(2)*g*mhat
    def mH_of_r(r):
        MR = r*MD; mu = 2*np.sqrt(2)*g*n/(n + 3)*MR
        Y = np.concatenate([yL4, [l1, 0, l3, 0, l5, np.log(mu)]]); ys = ys0.copy()
        h = (tD - tL)/600
        for _ in range(600):
            f = lambda Z: beta_tr(Z, ys)
            k1 = f(Y); k2 = f(Y + h/2*k1); k3 = f(Y + h/2*k2); k4 = f(Y + h*k3); Y = Y + h/6*(k1 + 2*k2 + 2*k3 + k4)
            ys = ys + h*beta_ys(Y, ys)
        mu_D = np.exp(Y[9]); lam = Y[4] - mu_D**2/MD**2
        stab = Y[4] > 0 and Y[5] + Y[6] > 0 and Y[5] + Y[6]/2 > 0 and Y[7] + 2*np.sqrt(max(Y[4]*(Y[5]+Y[6]), 0)) > 0
        yy = rk(beta_sm, np.concatenate([Y[:4], [lam]]), tD, 0, 800)
        return v*np.sqrt(2*max(yy[4], 1e-9)), Y, stab
    lo, hi = 0.2, 30.0
    for _ in range(36):
        mid = np.sqrt(lo*hi); m, *_ = mH_of_r(mid)
        (lo, hi) = (mid, hi) if m > 125 else (lo, mid)
    r = np.sqrt(lo*hi); m, Y, stab = mH_of_r(r)
    return r, Y, stab, g
for kap in (1/np.sqrt(2), 1.0):
    for mname, mhat in [('hier.', np.array([0, 0, 1.0])), ('degen.', np.ones(3)/np.sqrt(3))]:
        for n in (0.25, 0.5, 1.0):
            r, Y, stab, g = solve(1e16, n, 1e14, mhat, kap)
            r0 = np.sqrt(0.26*(n + 3)**2/(8*g**2*n**2))
            print(f'kap={kap:.3f} m^ {mname:6s} n={n:4.2f} : r = M_R/M_D = {r:5.2f} (sans RG au-dessus : ~{r0:4.2f}) ; m_II/m_I = 2r^2 = {2*r**2:5.1f} ; a M_D : l1 {Y[4]:.3f} l3 {Y[6]:.3f} l5 {Y[8]:.3f} stable {stab}', flush=True)
