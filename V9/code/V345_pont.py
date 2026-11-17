"""V345 K1 : chaîne de profondeur (modèle réduit de V330/V331) au point K, fibre jouet γ_F = σz, D_F = σx, insertion ε_D O ⊗ σx à la seule couche h."""
import numpy as np, sys
sys.path.insert(0, '.')
from rvA_common import clifford, H as Hb, KPT
C = clifford(); G5, al, A, B, G4 = C['G5'], C['al'], C['A'], C['B'], C['G4']
Y123 = al[0]@al[1]@al[2]; O = A@Y123; FG = 1j*A@B@G5
assert np.allclose(O@O, np.eye(8)) and np.allclose(O@G4, -G4@O) and np.allclose(O@FG, -FG@O) and np.allclose(O@G5, G5@O)
r = 1.0; mu = 1.0; Nz = 40
def Hsec(s, mK):
    T = np.diag(np.ones(Nz-1), -1)
    Hl = s*mK*G5 + mu*FG                                     # au point K (q = 0) : partie paire = masse ; terme de goût nu
    Hv = np.kron(-0.5j*(T - T.T), G4) + s*np.kron(0.5*r*(2*np.eye(Nz) - T - T.T), G5)
    return np.kron(np.eye(Nz), Hl) + Hv
def mass(eps, h, eD):
    mK = np.sqrt(3)*(eps - 1)
    P = np.zeros((Nz, Nz)); P[h, h] = 1
    X = eD*np.kron(P, O)
    H = np.block([[Hsec(+1, mK), X], [X, Hsec(-1, mK)]])
    E, V = np.linalg.eigh(H)
    W = (np.abs(V)**2).reshape(2, Nz, 8, -1).sum(axis=(0, 2))
    loc = W[:Nz//2].sum(0) > 0.5                              # localisés au bord n = 0
    return np.min(np.abs(E[loc]))
for eps in (0.8, 0.7):
    lam = np.sqrt(3)*(1 - eps); eD = 0.01
    ms = [mass(eps, h, eD) for h in range(8)]
    sl = np.polyfit(range(1, 7), np.log(ms[1:7]), 1)[0]
    print(f"ε = {eps} : λ′ = {lam:.3f} ; m(h), h = 0..7 : {np.array2string(np.array(ms), precision=6)}")
    print(f"   pente ln m / h (h = 1..6) = {sl:.4f} ; prédite 2 ln λ′ = {2*np.log(lam):.4f} ; écart {abs(sl/(2*np.log(lam)) - 1)*100:.1f} %")
    print(f"   m(0) = {ms[0]:.6f} ; prédit ε_D(1 − λ′²) = {eD*(1 - lam**2):.6f} ; écart {abs(ms[0]/(eD*(1-lam**2)) - 1)*100:.1f} %")
