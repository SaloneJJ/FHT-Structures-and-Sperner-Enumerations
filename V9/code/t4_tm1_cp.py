"""V69 : (a) relation de somme TM1 (U = U_TBM . R23(theta, phi)) -> delta_CP en fonction de theta23, theta13 ;
(b) quarks : melange entre bases propres de C3 fixant deux sommets differents. JJ SALONE, FHT 2026"""
import numpy as np, itertools
TBM = np.array([[np.sqrt(2/3), 1/np.sqrt(3), 0], [-1/np.sqrt(6), 1/np.sqrt(3), -1/np.sqrt(2)], [-1/np.sqrt(6), 1/np.sqrt(3), 1/np.sqrt(2)]])
def angles(U):
    s13 = abs(U[0, 2]); c13 = np.sqrt(1 - s13**2)
    s12 = abs(U[0, 1]) / c13; s23 = abs(U[1, 2]) / c13
    c12, c23 = np.sqrt(1 - s12**2), np.sqrt(1 - s23**2)
    J = np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))
    sd = J / (s12 * c12 * s23 * c23 * s13 * c13**2)
    # cos delta par |U_mu1|^2 = s12^2 c23^2 + c12^2 s23^2 s13^2 + 2 s12 c12 s23 c23 s13 cos d
    cd = (abs(U[1, 0])**2 - s12**2 * c23**2 - c12**2 * s23**2 * s13**2) / (2 * s12 * c12 * s23 * c23 * s13)
    return s12**2, s13**2, s23**2, cd, sd
pts = []
for th in np.linspace(0.05, 0.4, 400):
    for ph in np.linspace(-np.pi, np.pi, 721):
        R = np.array([[1, 0, 0], [0, np.cos(th), np.sin(th) * np.exp(-1j * ph)], [0, -np.sin(th) * np.exp(1j * ph), np.cos(th)]])
        pts.append(angles(TBM @ R))
pts = np.array(pts)
for s23t in (0.47, 0.50, 0.561):
    m = (np.abs(pts[:, 1] - 0.0220) < 4e-4) & (np.abs(pts[:, 2] - s23t) < 4e-3)
    if m.any():
        cd = pts[m, 3]; s12 = pts[m, 0]
        print(f'TM1, sin2th13 = 0.022, sin2th23 = {s23t} : sin2th12 = {s12.mean():.4f}, cos delta in [{cd.min():.3f}, {cd.max():.3f}] -> |delta| ~ {np.degrees(np.arccos(np.clip(cd.mean(), -1, 1))):.0f} deg')
# formule analytique candidate
for s23t in (0.47, 0.50, 0.561):
    s13 = np.sqrt(0.022); t23 = np.arcsin(np.sqrt(s23t))
    cd = -(1 - 5 * s13**2) / (2 * np.sqrt(2) * s13 * np.sqrt(1 - 3 * s13**2)) / np.tan(2 * t23)
    print(f'  formule : cos delta = {cd:.3f}')
# (b) quarks
S4 = list(itertools.permutations(range(4)))
Bz = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float); Bz = (Bz.T / np.linalg.norm(Bz, axis=1)).T
def rho(g):
    P = np.zeros((4, 4))
    for i in range(4): P[g[i], i] = 1
    return Bz @ P @ Bz.T
def basis(c):  # vecteurs propres du 3-cycle c
    w, U = np.linalg.eig(rho(c)); return U[:, np.argsort(np.angle(w))]
Uu = basis((0, 2, 3, 1)); Ud = basis((2, 1, 3, 0))  # C3 fixant 0, puis C3 fixant 1
print('quarks, sommets differents : |V|^2 =\n', np.round(np.abs(Uu.conj().T @ Ud)**2, 3))
