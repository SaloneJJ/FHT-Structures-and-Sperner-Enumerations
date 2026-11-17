"""V353 piste (e) : règles e1, e2 de 01_synthese/V353 §1 (écrites avant). numpy seul.   JJ SALONE, FHT 2026"""
import numpy as np
rng = np.random.default_rng(353)
def ckm(Mu, Md):
    _, Uu = np.linalg.eigh(Mu@Mu.conj().T); _, Ud = np.linalg.eigh(Md@Md.conj().T); V = Uu.conj().T@Ud
    J = np.imag(V[0, 1]*V[1, 2]*np.conj(V[0, 2])*np.conj(V[1, 1])); return np.abs(V), J
print("== e1 : phase de secteur commune aux générations : |V| et J indépendants de (φ_u, φ_d)")
for t in range(3):
    Yu, Yd = rng.normal(size=(3, 3)), rng.normal(size=(3, 3)); Yu = Yu + Yu.T; Yd = Yd + Yd.T
    V0, J0 = ckm(Yu, Yd)
    for pu, pd in ((0, np.pi/2), (0.3, 1.9), (np.pi/3, 0)):
        V1, J1 = ckm(np.exp(1j*pu)*Yu, np.exp(1j*pd)*Yd)
        print(f"  tirage {t}, (φ_u, φ_d) = ({np.rad2deg(pu):.0f}°, {np.rad2deg(pd):.0f}°) : max| |V|−|V0| | = {np.abs(V1 - V0).max():.1e} ; J = {J1:.1e} (J0 = {J0:.1e})")
print("== e2 : insertion A·Y123 ⊗ Y^A + B·Y123 ⊗ Y^B ⇒ bloc ⟨f=−1|·|f=+1⟩ = (bloc de O_0) ⊗ (Y^A − i Y^B)")
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
X = [k3(sx, I2, I2), k3(sz, sx, I2), k3(sz, sz, sx)]; Y = [k3(sy, I2, I2), k3(sz, sy, I2), k3(sz, sz, sy)]
A = (X[0]-X[1])/np.sqrt(2); B = (X[0]+X[1]-2*X[2])/np.sqrt(6); Y123 = Y[0]@Y[1]@Y[2]; F = 1j*A@B
w, V = np.linalg.eigh(F); Pp, Pm = V[:, w > 0], V[:, w < 0]
YA, YB = rng.normal(size=(3, 3)), rng.normal(size=(3, 3)); YA += YA.T; YB += YB.T
Ins = np.kron(A@Y123, YA) + np.kron(B@Y123, YB)
blk = np.kron(Pm, np.eye(3)).conj().T@Ins@np.kron(Pp, np.eye(3)); ref = Pm.conj().T@(A@Y123)@Pp
print("  écart au produit (bloc O_0) ⊗ (Y^A − iY^B) :", "%.1e" % np.abs(blk - np.kron(ref, YA - 1j*YB)).max())
Md = YA - 1j*YB; Mu = rng.normal(size=(3, 3)); Mu += Mu.T
print("  exemple : Y^A, Y^B réelles symétriques non proportionnelles, u réel : J =", "%.2e" % ckm(Mu, Md)[1], "(≠ 0 : l'octet fournit la phase)")
print("  contrôle : Y^B ∝ Y^A (direction E unique dans le secteur) : J =", "%.1e" % ckm(Mu, YA - 1j*0.7*YA)[1])
