"""V75 (1) : sans H3. Y_nu invariant par C3 (regle H1), reel : alpha P + beta (1-P) + gamma J (P projecteur sur l'axe de C3).
M_R et Y_Delta : textures symetriques Z2-invariantes (Z2 = (0 1), drapeau). m_nu = m_I + m_II, ||m_II||/||m_I|| = R.
Deviation de la colonne TM1 et de sin^2 th12 en fonction de R. JJ SALONE, FHT 2026"""
import numpy as np
rng = np.random.default_rng(11)
Bz = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float); Bz = (Bz.T / np.linalg.norm(Bz, axis=1)).T
def rho(g):
    P = np.zeros((4, 4))
    for i in range(4): P[g[i], i] = 1
    return Bz @ P @ Bz.T
C = rho((0, 2, 3, 1)); Z = rho((1, 0, 2, 3))
w, Ue = np.linalg.eig(C); Ue = Ue[:, np.argsort(np.angle(w))]
ax = np.real(Ue[:, np.argmin(np.abs(np.angle(w)))]); Pax = np.outer(ax, ax) / (ax @ ax)
Jr = (C - C.T) / np.sqrt(3)
def z2tex(spec):
    wz, Uz = np.linalg.eigh(Z); v1 = Uz[:, np.argmin(wz)] if np.sum(wz < 0) == 1 else Uz[:, np.argmax(wz)]
    Q = np.linalg.qr(np.column_stack([v1, rng.normal(size=(3, 2))]))[0]
    th = rng.uniform(0, np.pi); c, s = np.cos(th), np.sin(th); Q = Q @ np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return Q @ np.diag(np.array(spec)[rng.permutation(3)]) @ Q.T
TM1 = np.array([2/3, 1/6, 1/6])
def col_dev(m):
    wv, V = np.linalg.eigh(m)
    best = min(np.abs(np.sort(np.abs(Ue.conj().T @ V[:, k])**2) - np.sort(TM1)).max() for k in range(3))
    return best
s13 = 0.022
for spec_name, spec in [('M_R degenere (1,1,1)', (1, 1, 1)), ('M_R (0.5,1,1)', (0.5, 1, 1)), ('M_R hierarchique (0.1,0.5,1)', (0.1, 0.5, 1))]:
    print(spec_name)
    for R in (0.3, 1, 3, 10, 30, 100):
        devs = []
        for _ in range(300):
            al, be, ga = rng.uniform(0.5, 1.5, 2).tolist() + [rng.uniform(-0.5, 0.5)]
            Y = al * Pax + be * (np.eye(3) - Pax) + ga * Jr
            MR = z2tex(spec); YD = z2tex(rng.permutation([0.02, 0.3, 1.0]))
            mI = Y.T @ np.linalg.inv(MR) @ Y; mII = YD / np.linalg.norm(YD) * np.linalg.norm(mI) * R
            devs.append(col_dev(mI + mII))
        d = np.array(devs)
        # deviation de |U_e1|^2 = 2/3 : d(sin2 th12) ~ d / (cos2 th13) (ordre de grandeur)
        print(f'   R = {R:5.1f} : deviation colonne TM1 mediane {np.median(d):.3f}, q90 {np.quantile(d, .9):.3f}  ->  |d sin2 th12| ~ {np.median(d)/(1-s13):.3f}')
