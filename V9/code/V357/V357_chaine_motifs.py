"""V357 (relecture de V340) : tous les motifs η ∈ {±1}^(N+1) (η_0 = +1) sur la chaîne réduite de profondeur (R-unif, μ = r = 1, m_K = 0).
Base chirale (σ_z) : H = [[0, Q], [Q†, 0]] ; modes nuls de chiralité + : noyau de Q† ; de chiralité − : noyau de Q (SVD, seuil 1e-8).
Pour chaque motif : nombre de modes nuls (deux goûts), chiralité nette pondérée sur les deux niveaux les plus fins (N−1, N) et sur le niveau N."""
import numpy as np, itertools, sys
from V357_chaine_profondeur import chain
def chiral_zero(N, eta, f, s=+1, tol=1e-8):
    H = chain(N, 0.0, 1.0, eta, s, f)
    P = np.kron(np.eye(N + 1), np.array([[1, 1], [1, -1]]) / np.sqrt(2))       # σ_x -> σ_z ... base où σ_z est diagonale : on diagonalise directement
    sz = np.kron(np.eye(N + 1), np.diag([1., -1.]))
    ip = np.where(np.diag(sz) > 0)[0]; im = np.where(np.diag(sz) < 0)[0]
    Q = H[np.ix_(ip, im)]
    assert np.allclose(H[np.ix_(ip, ip)], 0) and np.allclose(H[np.ix_(im, im)], 0)
    U, sv, Vh = np.linalg.svd(Q)
    modes = []
    # noyau de Q (chiralité −) : lignes de Vh à valeur singulière nulle ; noyau de Q† (chiralité +) : colonnes de U
    for j in range(N + 1):
        if sv[j] < tol:
            modes.append((-1, np.abs(Vh[j].conj()) ** 2)); modes.append((+1, np.abs(U[:, j]) ** 2))
    return modes
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
rows = []
for m in itertools.product((1, -1), repeat=N):
    eta = (1,) + m; nz = 0; netN = 0.0; net2 = 0.0; deep = 0
    for f in (+1, -1):
        for c, p in chiral_zero(N, eta, f):
            nz += 1; netN += c * p[N]; net2 += c * p[N - 1:].sum()
            if p[1:N - 1].sum() > 0.5: deep += 1
    rows.append((eta, nz, round(netN, 3), round(net2, 3), deep))
ok = [r for r in rows if r[0][N - 1] == r[0][N]]; ko = [r for r in rows if r[0][N - 1] != r[0][N]]
def summary(rs):
    from collections import Counter
    return Counter((r[2], r[3]) for r in rs), Counter(r[1] for r in rs), Counter(r[4] for r in rs)
print(f'N = {N} : {len(rows)} motifs')
for nm, rs in (('η_(N−1) = η_N', ok), ('η_(N−1) ≠ η_N', ko)):
    a, b, c = summary(rs)
    print(f'  {nm} : (chiralité nette sur N, sur N−1..N) -> {dict(a)} ; nombre de modes nuls -> {dict(b)} ; modes intérieurs (1..N−2) -> {dict(c)}')
for r in rows[:: max(1, len(rows) // 12)]:
    print('   ', ''.join('+' if e > 0 else '-' for e in r[0]), 'modes nuls', r[1], 'net(N)', r[2], 'net(N−1..N)', r[3], 'intérieurs', r[4])
