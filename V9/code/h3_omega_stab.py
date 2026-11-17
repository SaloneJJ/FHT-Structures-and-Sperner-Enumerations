"""V74 : (a) omega et rapport type II / type I a trois generations, textures invariantes par le Z2 de TM1 ;
(b) quartiques canoniques (lambda1..5 du type II) issues de Tr Phi^4 (une generation de leptons), stabilite a Lambda.
JJ SALONE, FHT 2026"""
import numpy as np
rng = np.random.default_rng(7)
# (a) Z2 = transposition (0 1) dans la rep 3 (base standard), textures symetriques invariantes : rho^T M rho = M
Bz = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float); Bz = (Bz.T / np.linalg.norm(Bz, axis=1)).T
P = np.eye(4)[[1, 0, 2, 3]]; R = Bz @ P @ Bz.T
def z2_sym(Mlist):
    M = sum(Mlist) ; return (M + R.T @ M @ R) / 2
def rand_tex(spec):
    """matrice symetrique reelle Z2-invariante de spectre donne (vecteurs propres : v fixe par R + rotation libre)"""
    w, U = np.linalg.eigh(R); v1 = U[:, np.argmin(w)] if np.sum(w < 0) == 1 else U[:, np.argmax(w)]
    Q = np.linalg.qr(np.column_stack([v1, rng.normal(size=(3, 2))]))[0]
    th = rng.uniform(0, np.pi); c, s = np.cos(th), np.sin(th)
    Q = Q @ np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    perm = rng.permutation(3); return Q @ np.diag(np.array(spec)[perm]) @ Q.T
print('(a) Y_nu = y.1 ; r3 = M3/M_Delta = 2 ; rapport ||m_II||/||m_I|| et omega = <Yhat, M_R/||M_R||>')
for spec in [(1, 1, 1), (0.5, 1, 1), (0.1, 0.5, 1), (0.01, 0.1, 1)]:
    out = []
    for _ in range(200):
        MR = rand_tex(spec); Y = rand_tex(rng.permutation([0.01, 0.3, 1.0])); Y /= np.linalg.norm(Y)
        assert np.allclose(R.T @ MR @ R, MR) and np.allclose(R.T @ Y @ R, Y)
        MD = 1 / 2.0  # M3 = 1
        mI = np.linalg.inv(MR); mII = 2 / MD**2 * Y * np.sum(Y * MR)
        out.append((np.linalg.norm(mII) / np.linalg.norm(mI), np.sum(Y * MR) / np.linalg.norm(MR)))
    o = np.array(out)
    print(f'  spectre M_R {spec} : rapport median {np.median(np.abs(o[:, 0])):.2f} [q10 {np.quantile(np.abs(o[:, 0]), .1):.2f}, q90 {np.quantile(np.abs(o[:, 0]), .9):.2f}] ; |omega| median {np.median(np.abs(o[:, 1])):.2f}')
# (b) quartiques
s2 = np.array([[0, 1], [-1, 0]], complex)
def Phi(H, Dl, ynu, ye, yD):
    Ht = s2 @ H.conj(); S = np.zeros((4, 4), complex); S[0:2, 2] = ynu * Ht; S[0:2, 3] = ye * H; S = S + S.conj().T
    T = np.zeros((4, 4), complex); T[0:2, 0:2] = yD * (s2 @ Dl)
    return np.block([[S, T.conj().T], [T, S.conj()]])
ynu, ye, yD = 0.7, 0.3, 0.05
rows, ys = [], []
for _ in range(60):
    H = rng.normal(size=2) + 1j * rng.normal(size=2); d = rng.normal(size=3) + 1j * rng.normal(size=3)
    Dl = np.array([[d[1] / np.sqrt(2), d[2]], [d[0], -d[1] / np.sqrt(2)]])
    hh = np.real(H.conj() @ H); tdd = np.real(np.trace(Dl.conj().T @ Dl))
    rows.append([hh**2, hh * tdd, np.real(H.conj() @ Dl @ Dl.conj().T @ H), tdd**2, np.real(np.trace(Dl.conj().T @ Dl @ Dl.conj().T @ Dl))])
    M = Phi(H, Dl, ynu, ye, yD); ys.append(np.real(np.trace(np.linalg.matrix_power(M, 4))))
c, res, *_ = np.linalg.lstsq(np.array(rows), np.array(ys), rcond=None)
print('(b) Tr Phi^4 = c1 (H+H)^2 + c4 (H+H)Tr(D+D) + c5 H+DD+H + c2 (TrD+D)^2 + c3 Tr(D+DD+D) ; residu', res)
print('   c =', np.round(c, 6))
print('   / (ynu^4+ye^4, ., ., yD^4, yD^4) :', np.round(c[0] / (ynu**4 + ye**4), 4), np.round(c[3] / yD**4, 4), np.round(c[4] / yD**4, 4))
print('   c4/(yD^2), c5/(yD^2) en fonction de ynu^2, ye^2 :', np.round(c[1] / yD**2, 4), np.round(c[2] / yD**2, 4), ' (ynu^2 =', ynu**2, ', ye^2 =', ye**2, ')')
