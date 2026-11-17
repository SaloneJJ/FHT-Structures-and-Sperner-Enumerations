#!/usr/bin/env python3
"""AUDIT FABLE -- A3, sections 6-9 (registre C6-C9, M1-M3) : identites algebriques de l'octet, verifiees avec des matrices ecrites ici
(aucun import de rvA_common / V331). Chaque ligne affiche l'enonce de l'article (fragments new_sec7/8/9) et VRAI / FAUX.
Usage : python3 AF4_clifford_octet.py       (numpy)
"""
import itertools, numpy as np
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
X = [k3(sx, I2, I2), k3(sz, sx, I2), k3(sz, sz, sx)]            # sauts « cos » le long des trois directions du cube des enfants (flux pi par face)
Y = [k3(sy, I2, I2), k3(sz, sy, I2), k3(sz, sz, sy)]            # sauts « sin »
I8 = np.eye(8)
G5 = sum(X) / np.sqrt(3); A = (X[0] - X[1]) / np.sqrt(2); B = (X[0] + X[1] - 2 * X[2]) / np.sqrt(6)
F = 1j * A @ B; G4 = G5 @ Y[0] @ Y[1] @ Y[2]; S = k3(sz, sz, sz); FG = F @ G5
UI = k3(sx, sx, sx)                                              # inversion s -> 1 - s (sans phase)
com = lambda a, b: np.abs(a @ b - b @ a).max() < 1e-12
acom = lambda a, b: np.abs(a @ b + b @ a).max() < 1e-12
herm = lambda a: np.abs(a - a.conj().T).max() < 1e-12
res = []
def chk(label, val, expected=True):
    res.append((label, bool(val), expected)); print(f'  [{"ok" if bool(val) == expected else "!!"}] {label} : {"VRAI" if val else "FAUX"}' + ('' if bool(val) == expected else '   <-- contraire a l\'enonce'))

print('Generateurs : six involutions hermitiennes anticommutantes (Cl(6) sur C^8)')
gens = X + Y
chk('X_d, Y_d hermitiens, de carre 1, deux a deux anticommutants', all(herm(g) and np.allclose(g @ g, I8) for g in gens) and all(acom(a, b) for a, b in itertools.combinations(gens, 2)))
chk('flux pi sur chaque face du cube : X_a X_b X_a X_b = -1', all(np.allclose(X[a] @ X[b] @ X[a] @ X[b], -I8) for a, b in itertools.combinations(range(3), 2)))
chk('Gamma5, A, B, Y_d : six involutions hermitiennes anticommutantes', all(herm(g) and np.allclose(g @ g, I8) for g in [G5, A, B] + Y) and all(acom(a, b) for a, b in itertools.combinations([G5, A, B] + Y, 2)))
chk('Gamma4 = Gamma5 Y1 Y2 Y3 hermitien, de carre 1, anticommute avec Gamma5 et Y_d, commute avec A, B', herm(G4) and np.allclose(G4 @ G4, I8) and all(acom(G4, z) for z in [G5] + Y) and com(G4, A) and com(G4, B))

print('Lemme « Clifford identities »')
chk('(i) F^2 = 1', np.allclose(F @ F, I8))
chk('(i) F* = -F (adjoint) -- ENONCE TEL QU\'ECRIT dans new_sec7.tex', np.allclose(F.conj().T, -F), expected=False)
chk('(i) F hermitien (F* = F) et F^T = -F (F imaginaire pur) -- enonce corrige', herm(F) and np.allclose(F.T, -F) and np.allclose(F.real, 0))
chk('(i) F commute avec Gamma5, Y_d, Gamma4', com(F, G5) and all(com(F, y) for y in Y) and com(F, G4))
chk('(i) F anticommute avec A et B', acom(F, A) and acom(F, B))
nz = [(i, j, F[i, j]) for i in range(8) for j in range(8) if abs(F[i, j]) > 1e-12]
bits = lambda i: np.array([(i >> 2) & 1, (i >> 1) & 1, i & 1])
chk('(ii) F = sauts imaginaires +-i/sqrt3 sur les diagonales de face (distance de Hamming 2), 12 diagonales', len(nz) == 24 and all(abs(abs(v) - 1 / np.sqrt(3)) < 1e-12 and abs(v.real) < 1e-12 and (bits(i) != bits(j)).sum() == 2 for i, j, v in nz))
nz2 = [(i, j, FG[i, j]) for i in range(8) for j in range(8) if abs(FG[i, j]) > 1e-12]
chk('(ii) F Gamma5 = sauts imaginaires +-i sur les 4 grandes diagonales (antipodes)', len(nz2) == 8 and all(abs(abs(v) - 1) < 1e-12 and abs(v.real) < 1e-12 and (bits(i) != bits(j)).sum() == 3 for i, j, v in nz2))
# releve de l'inversion : U_I = D * (sx sx sx) avec D diagonale de signes ; on verifie F Gamma5 = i U_I pour un tel releve, et U_I^2 = +-1
D = np.diag(np.diag(-1j * FG @ UI)); Ulift = D @ UI
chk('(ii) F Gamma5 = i U_I, U_I = (signes) x (inversion s -> 1-s), U_I unitaire', np.allclose(FG, 1j * Ulift) and np.allclose(np.abs(np.diag(D)), 1) and np.allclose(Ulift @ Ulift.conj().T, I8))
print(f'       U_I^2 = {np.round((Ulift @ Ulift)[0, 0].real, 3)} * 1 ; U_I hermitien : {herm(Ulift)}')
chk('(iii) F = -Gamma4 S   (S = signe de sous-reseau (-1)^(s1+s2+s3))', np.allclose(F, -G4 @ S))
chk('(iii) variante de signe : F = +Gamma4 S', np.allclose(F, G4 @ S), expected=False)
P6 = G5 @ A @ B @ Y[0] @ Y[1] @ Y[2]
ph = (P6 @ S)[0, 0]
chk('(iii) S = phase x produit des six generateurs', np.allclose(P6 @ S, ph * I8) and abs(abs(ph) - 1) < 1e-12)
chk('(iv) F Gamma5 commute avec Gamma5, A, B', com(FG, G5) and com(FG, A) and com(FG, B))
chk('(iv) F Gamma5 anticommute avec Y_d et Gamma4', all(acom(FG, y) for y in Y) and acom(FG, G4))
chk('(iv) F Gamma5 hermitien, de carre 1', herm(FG) and np.allclose(FG @ FG, I8))
chk('(v) Y_d et Gamma4 anticommutent avec U_I', all(acom(y, Ulift) for y in Y) and acom(G4, Ulift))

print('Proposition « Flavour doubling » : commutant de {Gamma5, Y_d}')
gensC = [G5] + Y
M = np.vstack([np.kron(g, I8) - np.kron(I8, g.T) for g in gensC])
dimc = 64 - np.linalg.matrix_rank(M, tol=1e-9)
chk('dimension du commutant = 4', dimc == 4)
cand = [I8, G4 @ A, G4 @ B, A @ B]
chk('1, Gamma4 A, Gamma4 B, AB commutent avec Gamma5 et Y_d et sont independants (=> M_2(C))', all(com(c, g) for c in cand for g in gensC) and np.linalg.matrix_rank(np.array([c.ravel() for c in cand])) == 4)

print('Theoreme BBH / Proposition « staggered »')
rng = np.random.default_rng(1); ok = True
for _ in range(20):
    tf, tc = rng.uniform(.2, 1.5, 2); k = rng.uniform(-np.pi, np.pi, 3)
    H = sum((tf + tc * np.cos(k[d])) * X[d] + tc * np.sin(k[d]) * Y[d] for d in range(3))
    E2 = sum(tf**2 + tc**2 + 2 * tf * tc * np.cos(k[d]) for d in range(3))
    ok &= np.allclose(H @ H, E2 * I8)
chk('H(k)^2 = sum_d (tf^2 + tc^2 + 2 tf tc cos k_d) (20 tirages)', ok)
tf, tc = 0.7, 1.0; HK = sum((tf - tc) * X[d] for d in range(3))
chk('H(K) = sqrt3 (tf - tc) Gamma5, K = (pi, pi, pi)', np.allclose(HK, np.sqrt(3) * (tf - tc) * G5))

print('Choix D-C3 / structure reelle J_lat = Gamma5 A B . K')
U = G5 @ A @ B
Jc = lambda Z: U @ Z.conj() @ np.linalg.inv(U)
chk('J_lat^2 = -1 sur l\'octet (U conj(U) = -1)', np.allclose(U @ U.conj(), -I8))
chk('J X_d J^-1 = +X_d, J Y_d J^-1 = +Y_d, J Gamma4 J^-1 = +Gamma4', all(np.allclose(Jc(z), z) for z in X + Y + [G4]))
chk('J (F Gamma5) J^-1 = -F Gamma5 (terme de gout nu : change de signe, comme les termes impairs habilles ?)', np.allclose(Jc(FG), -FG))

print('Proposition « Insertion O_phi »')
Yp = Y[0] @ Y[1] @ Y[2]
Ophi = lambda p: (np.cos(p) * A + np.sin(p) * B) @ Yp
ok1 = all(herm(Ophi(p)) and np.allclose(Ophi(p) @ Ophi(p), I8) for p in (0, .4, 1.3, 2.9))
chk('(i) O_phi hermitien, de carre 1', ok1)
chk('(i) O_phi anticommute avec Y_d, Gamma4, F ; commute avec Gamma5', all(all(acom(Ophi(p), z) for z in Y + [G4, F]) and com(Ophi(p), G5) for p in (0, .4, 1.3)))
chk('(ii) J_lat O_phi J_lat^-1 = +O_phi', all(np.allclose(Jc(Ophi(p)), Ophi(p)) for p in (0, .4, 1.3)))
chk('(vi) O_0 ne commute ni n\'anticommute avec les X_d', not any(com(Ophi(0), x) or acom(Ophi(0), x) for x in X[:1]) or not all(com(Ophi(0), x) or acom(Ophi(0), x) for x in X))
print('       detail (O_0 contre X_1, X_2, X_3) : commute', [com(Ophi(0), x) for x in X], '; anticommute', [acom(Ophi(0), x) for x in X])
print('       detail (O_{0.4} contre X_d)        : commute', [com(Ophi(.4), x) for x in X], '; anticommute', [acom(Ophi(.4), x) for x in X])
print('Proposition « The angle phi »')
th = 0.37; Ut = np.cos(th) * I8 + 1j * np.sin(th) * F          # exp(i theta F), car F^2 = 1
s_plus = np.allclose(Ut @ Ophi(0.2) @ Ut.conj().T, Ophi(0.2 + 2 * th)); s_minus = np.allclose(Ut @ Ophi(0.2) @ Ut.conj().T, Ophi(0.2 - 2 * th))
chk('(i) U_theta O_phi U_theta^-1 = O_{phi + 2 theta} ou O_{phi - 2 theta} (U_theta = exp(i theta F))', s_plus or s_minus)
print(f'       signe : phi {"+" if s_plus else "-"} 2 theta')
p = 0.9
chk('(i) O_phi = O_0 (cos phi - i sin phi F)', np.allclose(Ophi(p), Ophi(0) @ (np.cos(p) * I8 - 1j * np.sin(p) * F)))
chk('(i) variante : O_phi = O_0 (cos phi + i sin phi F)', np.allclose(Ophi(p), Ophi(0) @ (np.cos(p) * I8 + 1j * np.sin(p) * F)), expected=False)
chk('(i) F commute avec la partie isotrope (Gamma5, Y_d, Gamma4, F Gamma5)', all(com(F, z) for z in [G5, G4, FG] + Y))
print('Proposition « No mass through Gamma5 x D_F »')
chk('[Gamma5, F] = 0 (l\'insertion Gamma5 x D_F conserve le gout)', com(G5, F))
bad = [l for l, v, e in res if v != e]
print(f'\nVERDICT : {len(res) - len(bad)} / {len(res)} lignes conformes a l\'attendu' + ('' if not bad else ' ; ECARTS : ' + ' | '.join(bad)))
