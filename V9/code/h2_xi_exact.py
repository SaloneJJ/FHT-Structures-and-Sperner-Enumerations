"""V73 : facteurs exacts c3, k du trilineaire mu H^T i s2 Delta^+ H issu de Tr(Phi^4) (une generation de leptons, J-double).
D = [[S, T^+],[T, conj(S)]] ; S : Yukawa (L <- nu_R : y_nu Ht, L <- e_R : y_e H) ; T (symetrique) : M_R (nu_R -> nu_R^c), y_D i s2 Delta (L -> L^c).
Calibration CCM : Tr Phi^2 = 4 a |H|^2, Tr Phi^4 = 4 b |H|^4 <-> cinetique (f0/2pi^2) a, quartique (f0/2pi^2) b. JJ SALONE, FHT 2026"""
import numpy as np
s2 = np.array([[0, 1], [-1, 0]], complex)
def Phi(h0, hp, d0, dp, dpp, ynu=0.7, ye=0.3, yD=1e-3, MR=1.0):
    H = np.array([hp, h0]); Ht = s2 @ H.conj()          # Ht = i s2 H*
    Dl = np.array([[dp / np.sqrt(2), dpp], [d0, -dp / np.sqrt(2)]])
    S = np.zeros((4, 4), complex)                        # base particules : nuL, eL, nuR, eR
    S[0:2, 2] = ynu * Ht; S[0:2, 3] = ye * H
    S = S + S.conj().T
    T = np.zeros((4, 4), complex)
    T[0:2, 0:2] = yD * (s2 @ Dl); T[2, 2] = MR
    assert np.allclose(T, T.T)
    return np.block([[S, T.conj().T], [T, S.conj()]])
tr = lambda M, p: np.real(np.trace(np.linalg.matrix_power(M, p)))
ynu, ye, yD, MR = 0.7, 0.3, 1e-3, 1.0
a = ynu**2 + ye**2
# calibration
h = 0.37
print('Tr Phi^2 / (4 a |H|^2) =', tr(Phi(h, 0, 0, 0, 0), 2) / (4 * a * h**2))
Tq = tr(Phi(h, 0, 0, 0, 0, MR=0), 4) / (4 * h**4); print('Tr Phi^4 (H seul) / 4|H|^4 =', Tq, ' b = ynu^4 + ye^4 =', ynu**4 + ye**4)
# normalisation du triplet : Tr Phi^2 = c_D |yD|^2 Tr(Delta^+ Delta)
for f in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
    d = 0.5; P = Phi(0, 0, *(d * np.array(f)), MR=0)
    print('Delta comp', f, ': Tr Phi^2 / (|yD|^2 |d|^2) =', tr(P, 2) / (yD**2 * d**2))
# trilineaire : terme de Tr Phi^4 lineaire en d0 et quadratique en h0, a M_R fixe
def coef(h0, d0):
    return tr(Phi(h0, 0, d0, 0, 0), 4) - tr(Phi(h0, 0, 0, 0, 0), 4) - tr(Phi(0, 0, d0, 0, 0), 4) + tr(Phi(0, 0, 0, 0, 0), 4)
for (h0, d0) in [(0.1, 0.01), (0.2, 0.01), (0.1, 0.02), (0.1, 0.01j)]:
    c = coef(h0, d0)
    print(f'h0={h0} d0={d0} : terme croise / (h0^2 Re d0 ou Im) =', c / (h0**2 * (d0.real if d0.real else d0.imag)), ' attendu ~ ynu^2 yD MR x facteur')

# ---- controles et consequences
print('calibration (M_R = 0) : Tr Phi^2 / (4 a |H|^2) =', tr(Phi(h, 0, 0, 0, 0, MR=0), 2) / (4 * a * h**2))
c = coef(0.1, 0.01) / (0.1**2 * 0.01)
print('Tr Phi^4 croise = %.4f x ynu^2 yD MR' % (c / (ynu**2 * yD * MR)))
# mu_brut = (f0/2pi^2) * (1/4) * 8 ynu^2 yD MR / 2 = (f0/2pi^2) ynu^2 yD MR ; Z_H = (f0/2pi^2) a ; Z_D = (f0/2pi^2) yD^2/2
# => mu_c = 2 sqrt2 g (ynu^2/a) M_R  (c3 = sqrt2, k = 1) ; m_II/m_I = 2 r^2 ; dl = 8 g^2 (n/(n+3))^2 r^2
# verification directe de m_II/m_I par minimisation (une generation), en unites ou f0/2pi^2 = 1/(4 g^2)
g = 0.52; F = 1 / (4 * g**2); n_ = ynu**2  # ici a = ynu^2 + ye^2
ZH, ZD = F * a, F * yD**2 / 2
for r in (1.0, 2.0, 3.0):
    MD = MR / r; muc = np.sqrt(2) * ynu**2 * MR / (a * np.sqrt(F))
    hc = 1e-3; dc = muc * hc**2 / MD**2
    mII = yD * dc / np.sqrt(ZD); mI = ynu**2 * (hc**2 / ZH) / MR
    print(f'r = {r} : m_II/m_I = {mII / mI:.4f}  (2 r^2 = {2 * r**2})')
