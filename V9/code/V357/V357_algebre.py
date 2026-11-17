"""V357 (relecture de V342 M1, M2, M4 et de V352 §3.5) : algèbre de l'octet, construite indépendamment (Jordan–Wigner à ordre inversé).
X_d réels, Y_d imaginaires purs (K : X -> X, Y -> −Y). Γ5 = ΣX/√3, A = (X1 − X2)/√2, B = (X1 + X2 − 2X3)/√6, F = iAB, Γ4 = Γ5 Y1Y2Y3, O_φ = (cos φ A + sin φ B) Y1Y2Y3."""
import numpy as np
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
X = [k3(I2, I2, sx) @ np.eye(8), k3(I2, sx, sz), k3(sx, sz, sz)]; X = [k3(sz, sz, sx), k3(sz, sx, I2), k3(sx, I2, I2)]   # ordre inversé
Y = [k3(sz, sz, sy), k3(sz, sy, I2), k3(sy, I2, I2)]
gens = X + Y
for i in range(6):
    for j in range(6):
        ac = gens[i] @ gens[j] + gens[j] @ gens[i]
        assert np.allclose(ac, 2 * np.eye(8) * (i == j)), (i, j)
G5 = sum(X) / np.sqrt(3); A = (X[0] - X[1]) / np.sqrt(2); B = (X[0] + X[1] - 2 * X[2]) / np.sqrt(6)
Y123 = Y[0] @ Y[1] @ Y[2]; G4 = G5 @ Y123; F = 1j * A @ B
O = lambda ph: (np.cos(ph) * A + np.sin(ph) * B) @ Y123
com = lambda a, b: np.linalg.norm(a @ b - b @ a); acom = lambda a, b: np.linalg.norm(a @ b + b @ a)
print('== M1')
for nm, Z in (('Γ5', G5), ('A·Y123', O(0)), ('B·Y123', O(np.pi / 2)), ('O_φ (φ = 0,7)', O(0.7))):
    print(f'  {nm:14s} herm {np.allclose(Z, Z.conj().T)}  Z²=1 {np.allclose(Z @ Z, np.eye(8))}  '
          f'{{Z,Y_d}}=0 {all(acom(Z, y) < 1e-12 for y in Y)}  {{Z,Γ4}}=0 {acom(Z, G4) < 1e-12}  [Z,Γ5]=0 {com(Z, G5) < 1e-12}  '
          f'{{Z,F}}=0 {acom(Z, F) < 1e-12}  [Z,F]=0 {com(Z, F) < 1e-12}  ; relation à X_d : comm {[round(com(Z, x), 3) for x in X]} anticomm {[round(acom(Z, x), 3) for x in X]}')
print('  F hermitien', np.allclose(F, F.conj().T), ' F²=1', np.allclose(F @ F, np.eye(8)), ' [F, X_d] = 0', all(com(F, x) < 1e-12 for x in X), ' [F, Γ4] = 0', com(F, G4) < 1e-12)
print('== M2 : J = U K, U = Γ5 A B')
U = G5 @ A @ B; print('  J² = U U* =', np.round((U @ U.conj())[0, 0], 6), ' (∝ 1 :', np.allclose(U @ U.conj(), (U @ U.conj())[0, 0] * np.eye(8)), ')')
for nm, Z in (('Γ5', G5), ('A·Y123', O(0)), ('B·Y123', O(np.pi / 2)), ('O_φ', O(0.7)), ('Y_d', Y[0]), ('X_d', X[1]), ('Γ4', G4), ('FΓ5', F @ G5)):
    JZ = U @ Z.conj() @ np.linalg.inv(U); sgn = '+' if np.allclose(JZ, Z) else ('−' if np.allclose(JZ, -Z) else '?')
    print(f'  J {nm} J⁻¹ = {sgn}{nm}')
print('== M4 : automorphismes de Cl(6) relevant les rotations (au niveau de Clifford)')
# C2 autour de l'axe 1 : Y2 -> −Y2, Y3 -> −Y3, X inchangés (k2, k3 -> −k2, −k3) ; relevé intérieur U = Y2 Y3
for ax in range(3):
    o = [d for d in range(3) if d != ax]; Uc = Y[o[0]] @ Y[o[1]]
    ok = all(np.allclose(Uc @ x @ np.linalg.inv(Uc), x) for x in X) and np.allclose(Uc @ Y[ax] @ np.linalg.inv(Uc), Y[ax]) and all(np.allclose(Uc @ Y[d] @ np.linalg.inv(Uc), -Y[d]) for d in o)
    print(f'  C2 (axe {ax + 1}) relevé correct {ok} ; O_φ invariant (φ = 0,7) : {np.allclose(Uc @ O(0.7) @ np.linalg.inv(Uc), O(0.7))}')
# C3 : permutation cyclique des axes 1 -> 2 -> 3 : on cherche U avec U X_d U⁻¹ = X_{d+1}, U Y_d U⁻¹ = Y_{d+1} (intertwiner, par moyenne de groupe)
src = X + Y; dst = [X[1], X[2], X[0], Y[1], Y[2], Y[0]]
import itertools
Ms = []
for r in range(1 << 6):
    P = np.eye(8, dtype=complex); Q = np.eye(8, dtype=complex)
    for i in range(6):
        if r >> i & 1: P = P @ src[i]; Q = Q @ dst[i]
    Ms.append((P, Q))
for seed in range(5):
    R = np.random.default_rng(seed).normal(size=(8, 8)) + 0j
    Uc = sum(Q @ R @ np.linalg.inv(P) for P, Q in Ms)
    if np.linalg.norm(Uc) > 1e-6: break
Uc = Uc / np.sqrt(abs(np.linalg.det(Uc)) ** (1 / 4))
print('  C3 : U X_d U⁻¹ = X_{d+1}', all(np.allclose(Uc @ s @ np.linalg.inv(Uc), d_) for s, d_ in zip(src, dst)))
for ph in (0.0, 0.7):
    Ot = Uc @ O(ph) @ np.linalg.inv(Uc)
    # projeter sur O_ψ : Ot = (cos ψ A + sin ψ B) Y123 ?
    c = np.real(np.trace(Ot @ np.linalg.inv(O(0)))) / 8; s_ = np.real(np.trace(Ot @ np.linalg.inv(O(np.pi / 2)))) / 8
    print(f'  C3 : O_φ (φ = {np.degrees(ph):.1f}°) -> O_ψ avec ψ = {np.degrees(np.arctan2(s_, c)) % 360:.1f}° (|composantes| = {np.hypot(c, s_):.3f}) ; Γ5 invariant {np.allclose(Uc @ G5 @ np.linalg.inv(Uc), G5)} ; F -> {np.round(np.real(np.trace(Uc @ F @ np.linalg.inv(Uc) @ F)) / 8, 3)}·F')
