"""Relecture R4 (V153 §3, V154 §1) avec F reconstruit (rel_F.py) et Gamma propres (rel_tree.py).
J = (J_F gamma_F) (x) J_X, J_X = U_X o conj, U_X construit pour la representation de la relecture."""
import sys, time, numpy as np, scipy.sparse as sp
from rel_F import *
from rel_tree import G, CHI, build
t0 = time.time(); rng = np.random.default_rng(2026)
DF = build_DF(rng)
# ---------- axiomes KO-6 de F (V146)
print('F : J_F^2 = 1 :', np.allclose(UF @ UF.conj(), np.eye(dF)), '; J_F gamma_F = -gamma_F J_F :', np.allclose(Jconj(UF, GF), -GF),
      '; J_F D_F = D_F J_F :', np.allclose(Jconj(UF, DF), DF), '; D_F impair :', np.allclose(DF @ GF + GF @ DF, 0), '; D_F hermitien :', np.allclose(DF, DF.conj().T))
basis = real_basis(); Jb = [Jconj(UF, b) for b in basis]
o0 = max(np.abs(a @ b - b @ a).max() for a in basis for b in Jb)
o1 = max(np.abs((DF @ a - a @ DF) @ b - b @ (DF @ a - a @ DF)).max() for a in basis for b in Jb)
gcomm = max(np.abs(GF @ a - a @ GF).max() for a in basis)
print(f'F : ordre zero max|[a, JbJ^-1]| = {o0:.1e} ; ordre un max|[[D_F,a],JbJ^-1]| = {o1:.1e} (base reelle complete 24 x 24) ; [gamma_F, a] = {gcomm:.1e}')
print('   (comparaison) max |[D_F, a]| =', f'{max(np.abs(DF @ a - a @ DF).max() for a in basis):.2f}')
# ---------- J_X dans la representation de la relecture : Gamma_1, Gamma_3, Gamma_4 reels ; Gamma_2, Gamma_5 imaginaires
real = [np.allclose(g.imag, 0) for g in G]; print('Gamma reels :', real)
UX = G[1] @ G[4]
print('U_X conj(Gamma_a) U_X^-1 = Gamma_a :', [bool(np.allclose(UX @ g.conj() @ np.linalg.inv(UX), g)) for g in G], '; J_X^2 =', np.round((UX @ UX.conj())[0, 0].real, 6),
      '; J_X chi J_X^-1 = -chi :', np.allclose(UX @ CHI.conj() @ np.linalg.inv(UX), -CHI))
UFG = UF @ GF; print('(J_F gamma_F)^2 =', np.round((UFG @ UFG.conj())[0, 0].real, 6))
# unicite de U_X (a une phase pres) : solutions de U conj(G_a) = G_a U
M = np.concatenate([np.kron(np.eye(4), g) - np.kron(g.conj().T, np.eye(4)) for g in G])   # vec en colonnes : (I (x) G) - (conj(G)^T (x) I)
print('dimension des U avec U conj(Gamma_a) = Gamma_a U :', 16 - np.linalg.matrix_rank(M))
# ---------- systeme complet, N = 3
N = 3; cells, index, K, W, depth = build(N, -1.0); Nc = len(cells)
eps = 0.2 / np.abs(np.linalg.eigvalsh(DF)).max()
G5 = sp.kron(sp.identity(Nc), sp.csr_matrix(G[4]))
H = sp.kron(K, sp.identity(dF)) + sp.kron(W, sp.csr_matrix(GF)) + eps * sp.kron(G5, sp.csr_matrix(DF)); H = H.tocsr()
U = sp.kron(sp.kron(sp.identity(Nc), sp.csr_matrix(UX)), sp.csr_matrix(UFG)).tocsr()
print(f'\nN = {N}, dim {H.shape[0]}, eps = {eps:.4f} : ||H - H^+|| = {abs(H - H.getH()).max():.1e} ; ||U conj(U) - 1|| = {abs(U @ U.conj() - sp.identity(H.shape[0])).max():.1e} (J^2 = +1) ;'
      f' ||U conj(H) U^+ + H|| = {abs(U @ H.conj() @ U.getH() + H).max():.1e} (JH = -HJ) ; ||U conj(H) U^+ - H|| = {abs(U @ H.conj() @ U.getH() - H).max():.2f}')
# ordre un pour A_F, base reelle complete
worst = 0.0
for a in basis[::3]:
    A1 = sp.kron(sp.identity(4 * Nc), sp.csr_matrix(a)); C = (H @ A1 - A1 @ H)
    for b in Jb[1::4]:
        B1 = sp.kron(sp.identity(4 * Nc), sp.csr_matrix(b)); worst = max(worst, abs(C @ B1 - B1 @ C).max())
print(f'ordre un pour A_F (8 x 6 elements de base) : {worst:.1e}')
# ---------- identite exacte H^2 = H0^2 + eps^2 (1 (x) D_F^2), H0 = K (x) 1 + W (x) gamma_F
H0 = (sp.kron(K, sp.identity(dF)) + sp.kron(W, sp.csr_matrix(GF))).tocsr()
R = H @ H - H0 @ H0 - eps ** 2 * sp.kron(sp.identity(4 * Nc), sp.csr_matrix(DF @ DF))
print(f'identite operatorielle ||H^2 - H0^2 - eps^2 D_F^2|| = {abs(R).max():.1e}  => E^2 = E0^2 + (eps mu)^2 pour TOUS les etats (bord, racine, partenaires, bulk)')
# et H0(gamma=-1) = - Gamma5 H0(gamma=+1) Gamma5 : spectres identiques au signe pres
Hp = (K + W); Hm = (K - W); print(f'Gamma5 H+ Gamma5 = -H- : {abs(G5 @ Hp @ G5 + Hm).max():.1e}')
# ---------- R4-P2 : J envoie les vecteurs propres de H+ sur ceux de H- avec E -> -E ; J_d^2 sur les modes nuls de bord
UXf = sp.kron(sp.identity(Nc), sp.csr_matrix(UX)).tocsr()
E, V = np.linalg.eigh(Hp.toarray()); sel = np.argsort(np.abs(E))[:40]
err = max(np.abs(Hm @ (UXf @ V[:, a].conj()) + E[a] * (UXf @ V[:, a].conj())).max() for a in sel)
chi_op = sp.kron(sp.identity(Nc), sp.csr_matrix(CHI))
flip = [(np.real(V[:, a].conj() @ (chi_op @ V[:, a])), np.real((UXf @ V[:, a].conj()).conj() @ (chi_op @ (UXf @ V[:, a].conj())))) for a in sel]
print(f'R4-P2 (N = 3, 40 plus petites |E| de H+) : max ||H-(J psi) + E J psi|| = {err:.1e} ; max |chi(J psi) + chi(psi)| = {max(abs(x + y) for x, y in flip):.1e}')
Jd = G5 @ UXf; print('J_d = Gamma5 J_X (identifie H- a H+) : J_d^2 =', np.round((Jd @ Jd.conj()).diagonal()[:4].real, 6))
print(f'temps {time.time()-t0:.0f} s')
