# Relecture V195 : modes de coin. Séparabilité exacte H = h⊗1⊗1 + Z⊗h⊗1 + Z⊗Z⊗h (jauge de Kawamoto–Smit, amas ouvert).
import numpy as np, itertools as it

def ssh(N, t_of_bond):
    h = np.zeros((N, N))
    for x in range(N - 1): h[x, x + 1] = h[x + 1, x] = t_of_bond(x)
    return h
def ks(d, r): return [1, (-1) ** r[0], (-1) ** (r[0] + r[1])][d]
def cluster3d(N, t_of_bond):
    sites = list(it.product(range(N), repeat=3)); idx = {s: i for i, s in enumerate(sites)}
    H = np.zeros((N ** 3,) * 2)
    for s in sites:
        for d in range(3):
            if s[d] + 1 < N:
                u = list(s); u[d] += 1; H[idx[s], idx[tuple(u)]] = H[idx[tuple(u)], idx[s]] = t_of_bond(s[d]) * ks(d, s)
    return H

# 1) identité opératorielle de séparabilité
N = 6; eps = 0.5; tb = lambda x: eps if x % 2 == 0 else 1.0
h = ssh(N, tb); Z = np.diag([(-1) ** x for x in range(N)]); I = np.eye(N)
Hsep = np.kron(np.kron(h, I), I) + np.kron(np.kron(Z, h), I) + np.kron(np.kron(Z, Z), h)
print("H_3D = h⊗1⊗1 + Z⊗h⊗1 + Z⊗Z⊗h (N=6) :", np.allclose(cluster3d(N, tb), Hsep), "; donc H² = Σ_d h_d² et E² = e_a²+e_b²+e_c²")

# 2) énergie des modes de coin : E = √3 δ_n exactement, δ_n = énergie du mode de bord de la chaîne SSH ouverte de 2n sites
print(" ε    n  | E (3D)      √3·δ_n (1D)  | √3(1−ε²)ε^n  | rapport E/[√3(1−ε²)ε^n]")
for eps in [0.3, 0.5, 0.7]:
    tb = lambda x, e=eps: e if x % 2 == 0 else 1.0
    for n in [2, 3, 4, 5, 8]:
        dn = np.sort(np.abs(np.linalg.eigvalsh(ssh(2 * n, tb))))[0]
        E3 = np.sort(np.abs(np.linalg.eigvalsh(cluster3d(2 * n, tb))))[7] if n <= 4 else np.nan
        f = np.sqrt(3) * (1 - eps ** 2) * eps ** n
        print(" %.1f  %d  | %.6f    %.6f     | %.6f      | %.4f" % (eps, n, E3, np.sqrt(3) * dn, f, np.sqrt(3) * dn / f))

# 3) base localisée des coins et octet effectif (E1-c et « octet tordu effectif »)
n = 3; eps = 0.5; tb = lambda x: eps if x % 2 == 0 else 1.0; N = 2 * n
w, V = np.linalg.eigh(ssh(N, tb)); o = np.argsort(np.abs(w))[:2]; P = V[:, o]
# modes de bord localisés : diagonaliser la position projetée
xop = P.T @ np.diag(np.arange(N)) @ P; _, U = np.linalg.eigh(xop); L = P @ U          # L[:,0] gauche, L[:,1] droite
for j in range(2): L[:, j] *= np.sign(L[np.argmax(np.abs(L[:, j])), j])
print("modes de bord 1D (n=3, ε=0.5) : poids pair du mode gauche = %.4f, du mode droit = %.4f"
      % (np.sum(L[::2, 0] ** 2), np.sum(L[::2, 1] ** 2)))
heff = L.T @ ssh(N, tb) @ L; zeff = L.T @ np.diag([(-1) ** x for x in range(N)]) @ L
print("h effectif 2×2 =", np.round(heff, 6).tolist(), "; Z effectif =", np.round(zeff, 6).tolist())
# octet effectif 3D dans la base produit (coins c ∈ {G,D}^3)
H3 = cluster3d(N, tb); B = np.zeros((N ** 3, 8))
for c, (a, b, cc) in enumerate(it.product(range(2), repeat=3)):
    B[:, c] = np.kron(np.kron(L[:, a], L[:, b]), L[:, cc])
Heff = B.T @ H3 @ B
print("octet effectif : H_eff² ∝ 1 : %s ; couplages non nuls entre coins à distance de Hamming 1 seulement : %s ; |t_f'| = %.6f"
      % (np.allclose(Heff @ Heff, (Heff @ Heff)[0, 0] * np.eye(8)),
         all(abs(Heff[i, j]) < 1e-12 or bin(i ^ j).count("1") == 1 for i in range(8) for j in range(8)), abs(Heff[0, 1])))
# flux par face du cube effectif : produit des signes autour de chaque face
faces = []
for ax in range(3):
    for v in range(2):
        cs = [c for c in range(8) if (c >> (2 - ax)) & 1 == v]
        a0 = cs[0]; nb = [c for c in cs if bin(c ^ a0).count("1") == 1]; opp = [c for c in cs if bin(c ^ a0).count("1") == 2][0]
        faces.append(np.sign(Heff[a0, nb[0]] * Heff[nb[0], opp] * Heff[opp, nb[1]] * Heff[nb[1], a0]))
print("produit des sauts autour des 6 faces de l'octet effectif :", faces, "(−1 = flux π)")
# chiralité de sous-réseau des 8 modes de coin localisés
par = np.array([sum(s) % 2 for s in it.product(range(N), repeat=3)])
print("poids sur le sous-réseau pair de chaque mode de coin localisé :", np.round([np.sum(B[par == 0, c] ** 2) for c in range(8)], 4).tolist())
print("parité du site de coin correspondant :", [((a * (N - 1)) + (b * (N - 1)) + (cc * (N - 1))) % 2 for a, b, cc in it.product(range(2), repeat=3)])
