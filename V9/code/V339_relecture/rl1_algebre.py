"""V339 relecture, partie 1 : contrôles algébriques (point 1 et point 4).
(a) orientation (signe de chiralité) du Cl(5) natif de V335 contre le modèle emprunté V162 ;
(b) convention « cos » : quelles matrices sont habillées ; (c) F = -Γ4 S et FΓ5 = i U_I ; (d) t_f = t_c : flux π par plaquette, 8 zéros en K, H² = q² ;
(e) structure réelle J = U K sur l'octet, compatible avec les termes verticaux pondérés (ω réel) et le terme nu."""
import numpy as np, sys
sys.path.insert(0, '.')
from rl_common import *
np.set_printoptions(precision=4, suppress=True)
c = clifford_own(); G5, al, A, B = c['G5'], c['al'], c['A'], c['B']
F = 1j * A @ B
print('anticommutation des 6 générateurs {G5, α1..3, A, B} :', max(np.abs(X @ Y + Y @ X).max() for i, X in enumerate([G5, *al, A, B]) for j, Y in enumerate([G5, *al, A, B]) if i < j))
# (a) Γ4 de V335 : G4 = G5 Y1 Y2 Y3 avec Yg = parties sin ; α_d (cinétique effectif au cône) = -t_c Y_d
sxm = np.array([[0, 1], [1, 0]], complex); sym = np.array([[0, -1j], [1j, 0]]); szm = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, cc: np.kron(np.kron(a, b), cc)
Xg = [k3(sxm, I2, I2), k3(szm, sxm, I2), k3(szm, szm, sxm)]; Yg = [k3(sym, I2, I2), k3(szm, sym, I2), k3(szm, szm, sym)]
G4_v335 = (sum(Xg) / np.sqrt(3)) @ Yg[0] @ Yg[1] @ Yg[2]
H0, M = octet_mats(1, 1)
print('mes X_d = M_d + M_d^T = Xg de V335 :', all(np.allclose(M[d] + M[d].T, Xg[d]) for d in range(3)), '; mes Y_d = i(M_d - M_d^T) = Yg :', all(np.allclose(1j * (M[d] - M[d].T), Yg[d]) for d in range(3)))
print('mon G5 = ΣXg/√3 :', np.allclose(G5, sum(Xg) / np.sqrt(3)))
for gauge in (False, True):
    cc_ = clifford_own(alt=gauge); a_ = cc_['al']
    print(f"jauge {'alternée (cône en Γ)' if gauge else 'K'} : α_d au cône = {np.round([np.trace(a_[d] @ Yg[d]).real / 8 for d in range(3)], 6)} × Y_d ; masse au cône (t_f=0,4) = {np.trace(Hk(np.zeros(3) if gauge else np.full(3, np.pi), 0.4, 1, gauge) @ sum(Xg)/np.sqrt(3)).real/8:+.4f} Γ5(V335)")
Gb = borrowed_gammas()
print('modèle emprunté : ε = tr(Γ1Γ2Γ3Γ4Γ5)/4 =', np.round(np.trace(Gb[0] @ Gb[1] @ Gb[2] @ Gb[3] @ Gb[4]).real / 4, 6))
# (b) parties paires / impaires : F commute avec G5, α, Γ4 ; anticommute avec A, B
for nm, Z in (('G5', G5), ('α1', al[0]), ('Γ4', G4_v335), ('A', A), ('B', B)):
    print(f'  [F, {nm}] = {np.abs(F @ Z - Z @ F).max():.1e}   {{F, {nm}}} = {np.abs(F @ Z + Z @ F).max():.1e}')
FG = F @ G5
# conventions exactes de V335 pour les identités (c) : Γ5 = ΣX/√3, A, B de V335, F = iAB, Γ4 = Γ5 Y1Y2Y3
G5v = sum(Xg) / np.sqrt(3); Av = (Xg[0] - Xg[1]) / np.sqrt(2); Bv = (Xg[0] + Xg[1] - 2 * Xg[2]) / np.sqrt(6); Fv = 1j * Av @ Bv; FGv = Fv @ G5v
print('mon Γ5 = −ΣXg/√3 (signe opposé à V335, sans effet) :', np.allclose(G5, -G5v), '; mon F = F de V335 :', np.allclose(F, Fv))
av = [-Yg[d] for d in range(3)]
print('orientation (V335) : ε_natif = tr(α1α2α3 Γ4 Γ5)/8 avec α_d = −Y_d (cinétique effectif), Γ4 = G4, Γ5 = G5 de V335 :',
      np.round(np.trace(av[0] @ av[1] @ av[2] @ G4_v335 @ G5v).real / 8, 6))
# (c) F = -Γ4 S ; FΓ5 = i U_I
S = np.diag(SUB)
print('conventions V335 : |F + Γ4 S| =', np.abs(Fv + G4_v335 @ S).max(), '; |F − Γ4 S| =', np.abs(Fv - G4_v335 @ S).max())
P = c['P']; P /= np.sqrt(np.real(np.trace(P @ P.conj().T)) / 8)
anti = np.zeros((8, 8))
for s in SITES: anti[IDX[tuple(1 - np.array(s))], IDX[s]] = 1
D = np.diag(np.diag(FGv @ anti.T) / 1j)
print('FΓ5 (V335) = i D P_antipode, D diagonale de signes :', np.allclose(FGv, 1j * D @ anti), np.round(np.diag(D).real, 3), '; D réel, D² = 1 :', np.allclose(D @ D, np.eye(8)))
print('conventions V335 : |FΓ5 + Γ4 S Γ5| =', np.abs(FGv + G4_v335 @ S @ G5v).max())
print('S = (Cl(6) chiralité) ? S ∝ X1X2X3Y1Y2Y3 :', np.round(np.trace(S @ Xg[0] @ Xg[1] @ Xg[2] @ Yg[0] @ Yg[1] @ Yg[2]) / 8, 6))
# (d) t_f = t_c : flux par plaquette du réseau uniforme
def hop(R1, R2):   # amplitude du saut entre sites globaux R1 -> R2 (voisins), réseau uniforme, jauge K
    x1, s1 = np.divmod(np.array(R1), 2); x2, s2 = np.divmod(np.array(R2), 2)
    Hm = H0 if np.all(x1 == x2) else None
    if Hm is None:
        d = int(np.nonzero(x2 - x1)[0][0]); Hm = M[d] if (x2 - x1)[d] == 1 else M[d].T
    return Hm[IDX[tuple(s1)], IDX[tuple(s2)]]
fl = []
for R in [np.array(v) for v in it.product(range(4), repeat=3)]:
    for a, b in ((0, 1), (0, 2), (1, 2)):
        ea, eb = np.eye(3, dtype=int)[a], np.eye(3, dtype=int)[b]
        p = hop(R, R + ea) * hop(R + ea, R + ea + eb) * hop(R + ea + eb, R + eb) * hop(R + eb, R)
        fl.append(np.sign(p))
print('t_f = t_c : produit des sauts sur les plaquettes élémentaires :', set(fl), '(−1 = flux π)')
print('spectre de H(K) à t_f = t_c :', np.round(np.linalg.eigvalsh(Hk(np.full(3, np.pi))), 8))
q = np.array([0.013, -0.02, 0.007]); H2 = Hk(np.full(3, np.pi) + q) @ Hk(np.full(3, np.pi) + q)
print('H(K+q)² - |q|² 1 (ordre 2) :', np.abs(H2 - (q @ q) * np.eye(8)).max())
# (e) structure réelle : J = U K, U = G5 A B
U = G5 @ A @ B
odd = list(al) + [G4_v335]; even = [G5, A, B] + [M[d] + M[d].T for d in range(3)]
# NB : en impulsion, une matrice « sin » Z (imaginaire) vient d'un terme d'espace réel (T - T^†) ⊗ (réel antisymétrique) ;
# J H_impair J⁻¹ = −H_impair en espace réel  <=>  U Z̄ U⁻¹ = +Z en impulsion. Pour le terme sur site FΓ5 : U conj(FΓ5) U⁻¹ = −FΓ5.
print('J² = U Ū =', np.round((U @ U.conj())[0, 0].real, 6), '; U Ū ∝ 1 :', np.allclose(U @ U.conj(), (U @ U.conj())[0, 0] * np.eye(8)))
print('impairs (α_d, Γ4) : U Z̄ U⁻¹ = +Z :', [np.allclose(U @ Z.conj() @ np.linalg.inv(U), Z) for Z in odd])
print('pairs (Γ5, A, B, X_d) : U Z̄ U⁻¹ = +Z :', [np.allclose(U @ Z.conj() @ np.linalg.inv(U), Z) for Z in even])
print('J FΓ5 J⁻¹ = −FΓ5 (terme nu) :', np.allclose(U @ FG.conj() @ np.linalg.inv(U), -FG))
print('les sauts réels de l\'arbre (T, S pondéré par ω = ±1, sans champ de fond) sont invariants par K : la même J vaut sur l\'arbre (vérifié numériquement dans rl2, θ = 0, q = 0).')
