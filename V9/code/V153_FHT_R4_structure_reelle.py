"""V153 R4 : structure reelle J du systeme complet (F de V146) x (arbre, fermions de paroi) ; plongement du Yukawa ; masses au bord.
Criteres et predictions : V153 par. 3 (ecrits avant calcul).   Usage : python3 V153_FHT_R4_structure_reelle.py
JJ SALONE, FHT 2026"""
import sys, time, itertools
import numpy as np, scipy.sparse as sp
t0 = time.time()
np.set_printoptions(precision=5, suppress=True, linewidth=160)
# --- modele F de V146 (meme graine), repris de V148_FHT_e_raccord.py
sys.argv = [sys.argv[0], 'none']
exec(open('V148_FHT_e_raccord.py').read().split("for p, fn in (('Ea'")[0])
# --- arbre et hamiltonien de paroi (memes definitions que V151_FHT_B1_paroi.py)
sxp = np.array([[0, 1], [1, 0]], complex); syp = np.array([[0, -1j], [1j, 0]]); szp = np.array([[1, 0], [0, -1]], complex); I2p = np.eye(2, dtype=complex)
GAM = [np.kron(sxp, sxp), np.kron(sxp, syp), np.kron(sxp, szp), np.kron(syp, I2p), np.kron(szp, I2p)]
CHI = np.kron(sxp, I2p)
r = 1.0; m = -1.0

def grid(s):
    I, J, K_ = np.meshgrid(np.arange(s), np.arange(s), np.arange(s), indexing='ij'); return I.ravel(), J.ravel(), K_.ravel()

def tree_ops(N):
    sides = [2 ** n for n in range(N + 1)]; offs = np.concatenate([[0], np.cumsum([s ** 3 for s in sides])]).astype(np.int64); Nc = int(offs[-1])
    rT = [[], [], []]; cT = [[], [], []]; rS, cS, vS = [], [], []
    for n, s in enumerate(sides):
        I, J, K_ = grid(s); base = offs[n] + (I * s + J) * s + K_
        for d, (a, b, c) in enumerate([((I + 1) % s, J, K_), (I, (J + 1) % s, K_), (I, J, (K_ + 1) % s)]):
            rT[d].append(base); cT[d].append(offs[n] + (a * s + b) * s + c)
        if n >= 1:
            q = sides[n - 1]; rS.append(base); cS.append(offs[n - 1] + ((I // 2) * q + J // 2) * q + K_ // 2); vS.append(np.full(len(base), 1 / np.sqrt(8)))
    T = [sp.coo_matrix((np.ones(Nc), (np.concatenate(rT[d]), np.concatenate(cT[d]))), shape=(Nc, Nc)).tocsr() for d in range(3)]
    S = sp.coo_matrix((np.concatenate(vS), (np.concatenate(rS), np.concatenate(cS))), shape=(Nc, Nc)).tocsr()
    depth = np.concatenate([np.full(s ** 3, n) for n, s in enumerate(sides)])
    return sides, offs, Nc, T, S, depth

def KW(T, S, Nc):
    kr = lambda A, g: sp.kron(A, sp.csr_matrix(g), format='csr')
    K = sum(kr((T[d] - T[d].T) / 2j, GAM[d]) for d in range(3)) + kr((S - S.T) / 2j, GAM[3])
    W = -(r / 2) * sum(kr(T[d] + T[d].T, GAM[4]) for d in range(3)) - (r / 2) * kr(S + S.T, GAM[4]) + (m + 4 * r) * kr(sp.identity(Nc), GAM[4])
    return K.tocsr(), W.tocsr()

def sector_basis(sides, offs, Nc, eps):
    rows, cols, vals = [], [], []; col = 0
    for n, s in enumerate(sides):
        if s == 1:
            if eps == (0, 0, 0): rows.append([offs[n]]); cols.append([col]); vals.append([1.0]); col += 1
            continue
        h = s // 2; I, J, K_ = grid(h); nrep = len(I)
        for dl in itertools.product((0, 1), repeat=3):
            sgn = (-1) ** (eps[0] * dl[0] + eps[1] * dl[1] + eps[2] * dl[2]) / np.sqrt(8)
            rows.append(offs[n] + ((I + dl[0] * h) * s + (J + dl[1] * h)) * s + (K_ + dl[2] * h)); cols.append(col + np.arange(nrep)); vals.append(np.full(nrep, sgn))
        col += nrep
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(Nc, col)).tocsr()

UX = GAM[1] @ GAM[3]                       # Gamma_2 Gamma_4
print('=== R4 : structure reelle ===')
print('  U_X Gamma_a^* U_X^-1 = Gamma_a (a = 1..5) :', [bool(np.allclose(UX @ g.conj() @ np.linalg.inv(UX), g)) for g in GAM])
print('  J_X^2 = U_X conj(U_X) =', np.round(UX @ UX.conj(), 6)[0, 0], '(attendu -1) ; J_X chi J_X^-1 = -chi :', bool(np.allclose(UX @ CHI.conj() @ np.linalg.inv(UX), -CHI)))
UFG = UF @ GF                               # partie unitaire de J_F gamma_F
print('  (J_F gamma_F)^2 =', np.round((UFG @ UFG.conj())[0, 0], 6), '(attendu -1) ; J_F D_F J_F^-1 = D_F :', bool(np.allclose(Japp(UF, DF), DF)), '; J_F gamma_F = -gamma_F J_F :', bool(np.allclose(Japp(UF, GF), -GF)))

# ---------- R4-P1 : identite J H J^-1 = -H avec le F complet, arbre N = 3
N = 3; sides, offs, Nc, T, S, depth = tree_ops(N); K, W = KW(T, S, Nc)
muF = np.abs(np.linalg.eigvalsh(DF)).max(); eps_y = 0.2 / muF
G5 = sp.kron(sp.identity(Nc), sp.csr_matrix(GAM[4]), format='csr')
H = sp.kron(K, sp.identity(dF), format='csr') + sp.kron(W, sp.csr_matrix(GF), format='csr') + eps_y * sp.kron(G5, sp.csr_matrix(DF), format='csr')
U = sp.kron(sp.kron(sp.identity(Nc), sp.csr_matrix(UX)), sp.csr_matrix(UFG), format='csr')
print(f'\n  R4-P1 (N = {N}, dim {H.shape[0]}, epsilon = {eps_y:.4f}) : hermiticite {abs(H - H.getH()).max():.1e} ; '
      f'max |U conj(H) U^+ + H| = {abs(U @ H.conj() @ U.getH() + H).max():.1e} ; max |U conj(U) - 1| = {abs(U @ U.conj() - sp.identity(U.shape[0])).max():.1e}')
# premier ordre par rapport a A_F : [[H, a], J b J^-1]
A1 = sp.kron(sp.identity(4 * Nc), sp.csr_matrix(alea_alg()), format='csr'); B1 = sp.kron(sp.identity(4 * Nc), sp.csr_matrix(Japp(UF, alea_alg())), format='csr')
C = H @ A1 - A1 @ H; print(f'  ordre un pour A_F : max |[[H, a], J b J^-1]| = {abs(C @ B1 - B1 @ C).max():.1e}')

# ---------- R4-P2 : J envoie les modes de bord (gamma_F = +1) sur ceux de gamma_F = -1, E -> -E (sans Yukawa)
N = 4; sides, offs, Nc, T, S, depth = tree_ops(N); K, W = KW(T, S, Nc)
Hp = (K + W).tocsr(); Hm = (K - W).tocsr(); UXf = sp.kron(sp.identity(Nc), sp.csr_matrix(UX), format='csr')
B = sp.kron(sector_basis(sides, offs, Nc, (1, 0, 0)), sp.identity(4), format='csr')
E, V = np.linalg.eigh((B.getH() @ Hp @ B).toarray()); sel = np.argsort(np.abs(E))[:16]
worst = 0.0; chis = []
for a in sel:
    psi = B @ V[:, a]; phi = UXf @ psi.conj()
    worst = max(worst, np.abs(Hm @ phi + E[a] * phi).max())
    p = phi.reshape(-1, 4); chis.append(np.real(np.einsum('ia,ab,ib->', p.conj(), CHI, p)))
print(f'\n  R4-P2 (N = 4, secteur 100, 16 plus petites |E| de H+) : max |H- (J psi) + E (J psi)| = {worst:.1e} ; helicite des images : [{min(chis):+.3f}, {max(chis):+.3f}]')

# ---------- R4-P3 : masses de bord.  SVD de la partie impaire de D_F -> blocs 2 x 2 exacts
Pp = np.where(gamma_F > 0)[0]; Pm = np.where(gamma_F < 0)[0]
M = DF[np.ix_(Pm, Pp)]; mus = np.linalg.svd(M, compute_uv=False)
print(f'\n  R4-P3 : valeurs singulieres de la partie impaire de D_F ({len(mus)}) : {np.round(mus, 3)}')
print(f'  epsilon = {eps_y:.4f} ; masses attendues epsilon*mu dans [{eps_y*mus.min():.4f}, {eps_y*mus.max():.4f}]')
sz2 = np.array([[1, 0], [0, -1]]); sx2 = np.array([[0, 1], [1, 0]])
B0 = sp.kron(sector_basis(sides, offs, Nc, (0, 0, 0)), sp.identity(4), format='csr'); K0 = (B0.getH() @ K @ B0).toarray(); W0 = (B0.getH() @ W @ B0).toarray(); G50 = (B0.getH() @ sp.kron(sp.identity(Nc), sp.csr_matrix(GAM[4])) @ B0).toarray()
B1s = B; K1 = (B1s.getH() @ K @ B1s).toarray(); W1 = (B1s.getH() @ W @ B1s).toarray(); G51 = (B1s.getH() @ sp.kron(sp.identity(Nc), sp.csr_matrix(GAM[4])) @ B1s).toarray()
E1_shell = np.sort(np.abs(np.linalg.eigvalsh(K1 + W1)))[:8]
print(f'  energies |E| les plus basses du secteur 100 sans Yukawa : {np.round(E1_shell, 5)}')
E1 = E1_shell[0]
for mu in [mus.min(), np.median(mus), mus.max()]:
    Hk0 = np.kron(K0, np.eye(2)) + np.kron(W0, sz2) + eps_y * mu * np.kron(G50, sx2)
    e0 = np.sort(np.abs(np.linalg.eigvalsh(Hk0)))[:8]
    Hk1 = np.kron(K1, np.eye(2)) + np.kron(W1, sz2) + eps_y * mu * np.kron(G51, sx2)
    e1 = np.sort(np.abs(np.linalg.eigvalsh(Hk1)))[:8]
    print(f'  mu = {mu:8.4f}, epsilon*mu = {eps_y*mu:.5f} : k = 0, 8 plus petites |E| = {np.round(e0, 5)}')
    print(f'      couche |n|^2 = 1 : 8 plus petites |E| = {np.round(e1, 5)} ; sqrt(E1^2 + (eps mu)^2) = {np.sqrt(E1**2 + (eps_y*mu)**2):.5f}')
print(f'\ntemps {time.time()-t0:.0f} s')
