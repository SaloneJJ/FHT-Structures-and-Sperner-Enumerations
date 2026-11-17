"""V77 (1) : Majorana REELLE dans la base geometrique (rep 3 de permutation), Z2-invariante (drapeau, TM1) ;
leptons charges : base propre complexe de C3. Predit-on theta23 = 45 deg et delta = +-90 deg ? JJ SALONE, FHT 2026"""
import numpy as np
rng = np.random.default_rng(3)
Bz = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float); Bz = (Bz.T / np.linalg.norm(Bz, axis=1)).T
def rho(g):
    P = np.zeros((4, 4))
    for i in range(4): P[g[i], i] = 1
    return Bz @ P @ Bz.T
C = rho((0, 2, 3, 1)); Z = rho((1, 0, 2, 3))
w, Ue = np.linalg.eig(C); Ue = Ue[:, np.argsort(np.angle(w))]
def obs(M):
    s13 = abs(M[0, 2]); c13 = np.sqrt(1 - s13**2); s12 = abs(M[0, 1]) / c13; s23 = abs(M[1, 2]) / c13
    c12, c23 = np.sqrt(1 - s12**2), np.sqrt(1 - s23**2)
    J = np.imag(M[0, 0] * M[1, 1] * np.conj(M[0, 1]) * np.conj(M[1, 0]))
    sd = J / (s12 * c12 * s23 * c23 * s13 * c13**2)
    cd = (abs(M[1, 0])**2 - s12**2 * c23**2 - c12**2 * s23**2 * s13**2) / (2 * s12 * c12 * s23 * c23 * s13)
    return s13**2, s12**2, s23**2, np.degrees(np.arctan2(sd, cd)) % 360
def arrange(U):
    P2 = np.abs(U)**2
    k1 = int(np.argmin([abs(np.sort(P2[:, k])[-1] - 2/3) for k in range(3)]))   # colonne TM1
    e = int(np.argmax(P2[:, k1]))
    rest = [k for k in range(3) if k != k1]; k3 = min(rest, key=lambda k: P2[e, k]); k2 = [k for k in rest if k != k3][0]
    rows = [e] + [i for i in range(3) if i != e]
    return U[np.ix_(rows, [k1, k2, k3])]
res = []
for _ in range(2000):
    A = rng.normal(size=(3, 3)); m = A + A.T; m = (m + Z.T @ m @ Z) / 2          # reelle symetrique, Z2-invariante
    ev, O = np.linalg.eigh(m)
    U = arrange(Ue.conj().T @ O)
    s13, s12, s23, d = obs(U)
    if 0.015 < s13 < 0.03: res.append((s13, s12, s23, d))
r = np.array(res)
print('tirages avec sin2th13 dans [0.015, 0.03] :', len(r))
print('sin2th23 : min %.4f max %.4f' % (r[:, 2].min(), r[:, 2].max()))
print('delta (deg) : valeurs distinctes ~', sorted(set(np.round(r[:, 3]).astype(int))))
print('sin2th12 : median %.4f' % np.median(r[:, 1]))
# controle : Majorana complexe generique Z2-invariante -> theta23, delta libres
res2 = []
for _ in range(2000):
    A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); m = A + A.T; m = (m + Z.T @ m @ Z) / 2
    H = m.conj().T @ m; ev, V = np.linalg.eigh(H); U = arrange(Ue.conj().T @ V)
    s13, s12, s23, d = obs(U)
    if 0.015 < s13 < 0.03: res2.append((s23, d))
r2 = np.array(res2); print('controle complexe : sin2th23 dans [%.3f, %.3f], delta dans [%.0f, %.0f]' % (r2[:, 0].min(), r2[:, 0].max(), r2[:, 1].min(), r2[:, 1].max()))
