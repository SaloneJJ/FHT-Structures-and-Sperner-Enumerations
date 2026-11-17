"""Relecture independante de B1 (V151/V152) : construction propre de l'hamiltonien de paroi dans l'arbre.
Representation des Gamma differente de celle des auteurs (Gamma_j = sigma_j (x) tau_z, Gamma_4 = 1 (x) tau_x, Gamma_5 = 1 (x) tau_y),
helicite chi = -i Gamma_4 Gamma_5 = 1 (x) tau_z. Diagonalisation dans l'espace COMPLET (pas de reduction par secteurs) :
dense pour N = 3, shift-invert (eigsh) pour N = 4.
Usage : python3 rel_B1.py N m [tree|flat] [k]"""
import sys, time, itertools
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
t0 = time.time()
N = int(sys.argv[1]); m = float(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'tree'
kk = int(sys.argv[4]) if len(sys.argv) > 4 else 200
r = 1.0
s0 = np.eye(2); s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1., -1.]).astype(complex)
G = [np.kron(s1, s3), np.kron(s2, s3), np.kron(s3, s3), np.kron(s0, s1), np.kron(s0, s2)]
for a in range(5):
    for b in range(5):
        assert np.allclose(G[a] @ G[b] + G[b] @ G[a], 2 * (a == b) * np.eye(4))
CHI = -1j * G[3] @ G[4]
print('Gamma1..5 =', np.round((G[0] @ G[1] @ G[2] @ G[3] @ G[4])[0, 0], 6), '(auteurs : -1 aussi => meme classe de representation de Cl(5))')
assert np.allclose(CHI, np.kron(s0, s3))

# ---- cellules
side = (lambda n: 2 ** N) if mode == 'flat' else (lambda n: 2 ** n)
cells = []; index = {}
for n in range(N + 1):
    s = side(n)
    for c in itertools.product(range(s), repeat=3):
        index[(n, c)] = len(cells); cells.append((n, c))
Nc = len(cells); depth = np.array([n for n, c in cells])
blocks = {}   # (i, j) -> 4x4
def add(i, j, M):
    blocks[(i, j)] = blocks.get((i, j), 0) + M
hop_lat = [G[d] / 2j - (r / 2) * G[4] for d in range(3)]   # coefficient de psi(x + e_d) dans (H psi)(x)
hop_ver = G[3] / 2j - (r / 2) * G[4]                       # coefficient de psi(parent) dans (H psi)(enfant), x poids
for i, (n, c) in enumerate(cells):
    s = side(n)
    add(i, i, (m + 4 * r) * G[4])
    for d in range(3):
        c2 = list(c); c2[d] = (c2[d] + 1) % s; j = index[(n, tuple(c2))]
        add(i, j, hop_lat[d]); add(j, i, hop_lat[d].conj().T)
    if n >= 1:
        if mode == 'flat': p = index[(n - 1, c)]; wgt = 1.0
        else: p = index[(n - 1, tuple(x // 2 for x in c))]; wgt = 1 / np.sqrt(8)
        add(i, p, wgt * hop_ver); add(p, i, wgt * hop_ver.conj().T)
rows, cols, vals = [], [], []
for (i, j), M in blocks.items():
    for a in range(4):
        for b in range(4):
            if M[a, b] != 0: rows.append(4 * i + a); cols.append(4 * j + b); vals.append(M[a, b])
H = sp.csr_matrix((vals, (rows, cols)), shape=(4 * Nc, 4 * Nc))
print(f'N = {N}, m = {m}, {mode} : {Nc} cellules, dim {4*Nc}, ||H - H^+|| = {abs(H - H.getH()).max():.1e}')

# ---- symetries : demi-translations (construction propre)
def half_perm(d):
    P = np.zeros(Nc, int)
    for i, (n, c) in enumerate(cells):
        s = side(n); c2 = list(c); c2[d] = (c2[d] + (s // 2 if s > 1 else 0)) % s; P[i] = index[(n, tuple(c2))]
    M = sp.csr_matrix((np.ones(Nc), (np.arange(Nc), P)), shape=(Nc, Nc))
    return sp.kron(M, sp.identity(4), format='csr')
for d in range(3):
    U = half_perm(d); print(f'  ||[tau_{d}, H]|| = {abs(U @ H - H @ U).max():.1e}', end='')
print()

# ---- diagonalisation dans l'espace complet
if 4 * Nc <= 4000:
    E, V = np.linalg.eigh(H.toarray())
else:
    E, V = sla.eigsh(H.tocsc(), k=kk, sigma=1.0e-4, which='LM', tol=1e-10)
o = np.argsort(E); E, V = E[o], V[:, o]
V, _R = np.linalg.qr(V)          # eigsh ne garantit pas l'orthonormalite dans les sous-espaces degeneres
V = V * np.sign(np.real(np.diag(_R)))[None, :]
print(f'  residu max ||H v - E v|| = {np.abs(H @ V - V * E[None, :]).max():.1e}')
# rotation dans chaque sous-espace degenere : diagonaliser chi + 1e-3 * profondeur (leve le melange arbitraire)
OBS = sp.kron(sp.diags(depth.astype(float)), sp.identity(4)) * 1e-3 + sp.kron(sp.identity(Nc), sp.csr_matrix(CHI))
i0 = 0
while i0 < len(E):
    i1 = i0 + 1
    while i1 < len(E) and abs(E[i1] - E[i0]) < 1e-7: i1 += 1
    if i1 - i0 > 1:
        Vb = V[:, i0:i1]; C = Vb.conj().T @ (OBS @ Vb); w_, R_ = np.linalg.eigh((C + C.conj().T) / 2); V[:, i0:i1] = Vb @ R_
    i0 = i1
o = np.argsort(np.abs(E)); E, V = E[o], V[:, o]
print(f'  {len(E)} valeurs propres ; max |E| calcule = {np.abs(E).max():.4f}   ({time.time()-t0:.0f} s)')
L = 2 ** N; win = 1.1 * np.sin(2 * np.pi * 2 / L)
rec = []
for a in range(len(E)):
    p = V[:, a].reshape(Nc, 4); w = np.sum(np.abs(p) ** 2, 1)
    wn = np.bincount(depth, weights=w, minlength=N + 1); nb = float(np.arange(N + 1) @ wn)
    chi = float(np.real(np.einsum('ia,ab,ib->', p.conj(), CHI, p)))
    rec.append((E[a], nb, chi, wn))
def classe(nb): return 'bord' if nb >= N - 1.5 else ('racine' if nb <= 1.0 else 'interieur')
def count(cut, lab):
    cnt = {}
    for (e, nb, chi, wn) in rec:
        if abs(e) < cut:
            k = (classe(nb), '+' if chi > 0.3 else ('-' if chi < -0.3 else '0')); cnt[k] = cnt.get(k, 0) + 1
    print(f'  {lab} |E| < {cut:.4f} : ' + ' ; '.join(f'{a} {b}: {v}' for (a, b), v in sorted(cnt.items())))
count(win, 'fenetre des auteurs')
if np.abs(E).max() > 0.9: count(0.9, 'fenetre elargie')
# niveaux distincts (multiplicites) pour E > 0, helicite +
lev = []
for (e, nb, chi, wn) in sorted(rec, key=lambda x: x[0]):
    if e > 1e-6:
        if lev and abs(e - lev[-1][0]) < 1e-6: lev[-1][1] += 1; lev[-1][2].append(chi); lev[-1][3].append(nb)
        else: lev.append([e, 1, [chi], [nb]])
print('  niveaux E > 0 (E, multiplicite, chi moyen, n_moyen moyen) :')
for e, mult, chis, nbs in lev[:14]:
    print(f'    {e:.5f}  x{mult:2d}  chi {np.mean(chis):+.3f}  n_moy {np.mean(nbs):.2f}')
zero = [x for x in rec if abs(x[0]) < 1e-5]
print(f'  etats |E| < 1e-5 : {len(zero)} : ' + ', '.join(f'(chi {c:+.2f}, n_moy {nb:.2f})' for e, nb, c, wn in zero))
neg = [x for x in rec if x[2] < -0.3 and x[1] > 1.0]
print(f'  etats d helicite - hors racine (n_moy > 1) calcules : {len(neg)} ; les plus bas : ' + ', '.join(f'{abs(e):.5f}(n {nb:.2f})' for e, nb, c, wn in sorted(neg, key=lambda x: abs(x[0]))[:14]))
mixed = [x for x in rec if abs(x[2]) <= 0.3 and abs(x[0]) > 1e-5]
print(f'  etats d helicite melangee |chi| <= 0.3 (hors nuls) : {len(mixed)} ; plus bas : ' + ', '.join(f'{abs(e):.4f}' for e, nb, c, wn in sorted(mixed, key=lambda x: abs(x[0]))[:6]))
if neg:
    e, nb, c, wn = min(neg, key=lambda x: abs(x[0])); print(f'  profil du premier partenaire : w(n) = {np.round(wn, 3)}')
print(f'temps {time.time()-t0:.0f} s')
