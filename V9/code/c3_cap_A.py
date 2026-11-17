"""V66 : cap-produit de la piece A (V25) : H^1(A) x H_k(A) -> H_{k-1}(A), candidat Yukawa (singulet -> doublet -> singulet).
Homologie reelle non reduite, harmoniques ; Alexander-Whitney : [v0..vk] cap phi = phi(v0 v1) [v1..vk]. JJ SALONE, FHT 2026"""
import sys, json, itertools, numpy as np
sys.path.insert(0, '.')
from fht_tools import actside
d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'V25_K3_pieces.json'))
A = [tuple(a) for a in d['A']]
V = sorted(set(v for a in A for v in a)); o = {v: i for i, v in enumerate(V)}
faces = set()
for a in A:
    a = tuple(sorted(a, key=o.get))
    for k in range(1, len(a) + 1): faces.update(itertools.combinations(a, k))
byd = {}
for s in faces: byd.setdefault(len(s) - 1, []).append(s)
for k in byd: byd[k].sort(key=lambda s: tuple(map(o.get, s)))
ind = {k: {s: i for i, s in enumerate(byd[k])} for k in byd}
top = max(byd)
def bd(k):
    M = np.zeros((len(byd[k-1]), len(byd[k])))
    for j, s in enumerate(byd[k]):
        for t in range(len(s)): M[ind[k-1][s[:t] + s[t+1:]], j] = (-1) ** t
    return M
B = {k: bd(k) for k in range(1, top + 1)}
H = {}
for k in range(top + 1):
    L = (B[k+1] @ B[k+1].T if k + 1 in B else 0) + (B[k].T @ B[k] if k in B else 0)
    L = np.zeros((len(byd[k]),) * 2) + L
    w, U = np.linalg.eigh(L); H[k] = U[:, w < 1e-9]
print('Betti non reduits', {k: H[k].shape[1] for k in H})
def cap(phi, k):
    """matrice C_k -> C_{k-1} de c -> c cap phi (phi 1-cochaine indexee par byd[1])"""
    M = np.zeros((len(byd[k-1]), len(byd[k])))
    for j, s in enumerate(byd[k]):
        M[ind[k-1][s[1:]], j] = phi[ind[1][s[:2]]]
    return M
red0 = np.ones((len(byd[0]), 1)) / np.sqrt(len(byd[0]))   # direction constante (H_0 non reduit -> reduit)
for i in range(H[1].shape[1]):
    phi = H[1][:, i]
    M21 = H[1].T @ cap(phi, 2) @ H[2]
    M10 = H[0].T @ cap(phi, 1) @ H[1]
    P0 = H[0].T @ red0; M10r = M10 - P0 @ (P0.T @ M10)  # partie reduite
    print(f'phi_{i} : H2->H1 {np.round(M21.ravel(), 4)}   H1->H0(red.) norme {np.linalg.norm(M10r):.4f} rang {np.linalg.matrix_rank(M10r, 1e-9)}')
# cup H^1 x H^1 -> H^2 (Alexander-Whitney : (a u b)(v0 v1 v2) = a(v0 v1) b(v1 v2))
for i, j in itertools.product(range(H[1].shape[1]), repeat=2):
    a, b = H[1][:, i], H[1][:, j]
    c = np.array([a[ind[1][s[:2]]] * b[ind[1][s[1:]]] for s in byd[2]])
    print(f'cup(phi_{i}, phi_{j}) -> H^2 : {np.round(H[2].T @ c, 4)}')
# composantes et localisation des classes
comp = {v: v for v in V}
def f(x):
    while comp[x] != x: x = comp[x]
    return x
for a in A:
    for v in a[1:]: comp[f(v)] = f(a[0])
cs = {}
for v in V: cs.setdefault(f(v), []).append(v)
print('composantes', list(cs.values()))
for k in (1, 2):
    for i in range(H[k].shape[1]):
        supp = set(v for j, s in enumerate(byd[k]) if abs(H[k][j, i]) > 1e-9 for v in s)
        print(f'classe H_{k}[{i}] portee par', sorted(supp, key=o.get))
