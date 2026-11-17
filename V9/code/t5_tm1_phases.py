"""V76 : TM1 avec phase phi algebrique (racine de l'unite) : prediction de (sin2 th23, delta_CP) a sin2 th13 fixe.
U = U_TBM . R23(theta, phi). JJ SALONE, FHT 2026"""
import numpy as np
TBM = np.array([[np.sqrt(2/3), 1/np.sqrt(3), 0], [-1/np.sqrt(6), 1/np.sqrt(3), -1/np.sqrt(2)], [-1/np.sqrt(6), 1/np.sqrt(3), 1/np.sqrt(2)]])
def U(th, ph):
    R = np.array([[1, 0, 0], [0, np.cos(th), np.sin(th) * np.exp(-1j * ph)], [0, -np.sin(th) * np.exp(1j * ph), np.cos(th)]])
    return TBM @ R
def obs(M):
    s13 = abs(M[0, 2]); c13 = np.sqrt(1 - s13**2); s12 = abs(M[0, 1]) / c13; s23 = abs(M[1, 2]) / c13
    c12, c23 = np.sqrt(1 - s12**2), np.sqrt(1 - s23**2)
    J = np.imag(M[0, 0] * M[1, 1] * np.conj(M[0, 1]) * np.conj(M[1, 0]))
    sd = J / (s12 * c12 * s23 * c23 * s13 * c13**2)
    cd = (abs(M[1, 0])**2 - s12**2 * c23**2 - c12**2 * s23**2 * s13**2) / (2 * s12 * c12 * s23 * c23 * s13)
    return s13**2, s12**2, s23**2, np.degrees(np.arctan2(sd, cd)) % 360
for name, ph in [('0', 0), ('pi/6', np.pi/6), ('pi/4', np.pi/4), ('pi/3', np.pi/3), ('pi/2', np.pi/2), ('2pi/3', 2*np.pi/3), ('3pi/4', 3*np.pi/4), ('5pi/6', 5*np.pi/6), ('pi', np.pi)]:
    ths = np.linspace(0.01, 0.6, 20000); best = min(ths, key=lambda t: abs(obs(U(t, ph))[0] - 0.0220))
    s13, s12, s23, d = obs(U(best, ph))
    s13b, s12b, s23b, db = obs(U(best, -ph))
    print(f'phi = {name:5s} : sin2th23 = {s23:.3f} ; delta = {d:6.1f} deg  (phi -> -phi : delta = {db:6.1f}, sin2th23 = {s23b:.3f})')
