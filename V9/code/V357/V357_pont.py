"""V357 (relecture de V342 M3 et V345 K1) : chaîne réduite de profondeur, deux secteurs γ_F = s = ±1, espace interne spineur ⊗ goût :
Γ4 = σ_y⊗1, Γ5 = σ_x⊗1, F = 1⊗τ_z, O = A·Y1Y2Y3 réduit = σ_x⊗τ_x (commute avec Γ5, anticommute avec Γ4 et F, O² = 1), contrôle O = Γ5.
H_s = (S−S†)/2i ⊗ Γ4 + s[(r/2)(2−S−S†) + m_K] ⊗ Γ5 + μ FΓ5 (nu) ; H = H_+ ⊕ H_− + ε_D (P_h ⊗ O) ⊗ σ_x (insertion au seul niveau N − h).
Masse du survivant = plus petite |E| ; prédiction (V345) : ε_D (1 − λ²) λ^(2h), λ = |1 + m_eff/r| de la copie qui porte le Weyl."""
import numpy as np
s1 = np.eye(2); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.])
G4 = np.kron(sy, s1); G5 = np.kron(sx, s1); F = np.kron(s1, sz); OAY = np.kron(sx, sx)
assert np.allclose(OAY @ G5, G5 @ OAY) and np.allclose(OAY @ G4, -G4 @ OAY) and np.allclose(OAY @ F, -F @ OAY)
def Hs(N, mK, s, mu=1.0, r=1.0):
    L = N + 1; S = np.diag(np.ones(N), -1)
    return np.kron((S - S.T) / 2j, G4) + s * np.kron((r / 2) * (2 * np.eye(L) - S - S.T) + mK * np.eye(L), G5) + mu * np.kron(np.eye(L), F @ G5)
def masse(N, eps, epsD, h, O=OAY):
    mK = np.sqrt(3) * (eps - 1); L = N + 1
    P = np.zeros((L, L)); P[N - h, N - h] = 1.0
    X = epsD * np.kron(P, O)
    H = np.block([[Hs(N, mK, +1), X], [X, Hs(N, mK, -1)]])
    E, V = np.linalg.eigh(H)
    # états de bord : poids > 0,5 sur les 10 niveaux les plus fins
    w = (np.abs(V) ** 2).reshape(2, L, 4, -1).sum(axis=(0, 2))
    sel = w[N - 9:].sum(0) > 0.5
    return np.abs(E[sel]).min() if sel.any() else np.abs(E).min()
N = 40; epsD = 1e-3
for eps in (0.9, 0.8, 0.7, 0.6, 0.2, 0.1):
    mK = np.sqrt(3) * (eps - 1); lam = abs(mK) if eps > 1 - 1 / np.sqrt(3) else abs(2 + mK)
    ms = [masse(N, eps, epsD, h) for h in range(6)]
    pred = [epsD * (1 - lam ** 2) * lam ** (2 * h) for h in range(6)]
    g5 = masse(N, eps, epsD, 0, O=G5)
    print(f'ε = {eps:.2f} (λ = {lam:.3f}, {"survivant" if eps > 0.4226 else "jumeau"}) : m(h)/préd = {np.round(np.array(ms) / np.array(pred), 4).tolist()} ; '
          f'm(0) = {ms[0]:.3e} ; contrôle O = Γ5 : {g5:.1e}')
