"""Constructeur independant de l'arbre augmente (T^3, profondeurs 0..N) et des parties K (cinetique) et W (masse + Wilson)
de l'hamiltonien de paroi, avec la representation des Gamma de la relecture."""
import itertools, numpy as np, scipy.sparse as sp
s0 = np.eye(2); s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1., -1.]).astype(complex)
G = [np.kron(s1, s3), np.kron(s2, s3), np.kron(s3, s3), np.kron(s0, s1), np.kron(s0, s2)]
CHI = -1j * G[3] @ G[4]

def build(N, m, r=1.0, flat=False):
    side = (lambda n: 2 ** N) if flat else (lambda n: 2 ** n)
    cells = []; index = {}
    for n in range(N + 1):
        for c in itertools.product(range(side(n)), repeat=3):
            index[(n, c)] = len(cells); cells.append((n, c))
    Nc = len(cells)
    Kr, Kc, Kv, Wr, Wc, Wv = [], [], [], [], [], []
    def put(R, C, V, i, j, M):
        for a in range(4):
            for b in range(4):
                if M[a, b] != 0: R.append(4 * i + a); C.append(4 * j + b); V.append(M[a, b])
    for i, (n, c) in enumerate(cells):
        s = side(n)
        put(Wr, Wc, Wv, i, i, (m + 4 * r) * G[4])
        for d in range(3):
            c2 = list(c); c2[d] = (c2[d] + 1) % s; j = index[(n, tuple(c2))]
            put(Kr, Kc, Kv, i, j, G[d] / 2j); put(Kr, Kc, Kv, j, i, -G[d] / 2j)
            put(Wr, Wc, Wv, i, j, -r / 2 * G[4]); put(Wr, Wc, Wv, j, i, -r / 2 * G[4])
        if n >= 1:
            p = index[(n - 1, c if flat else tuple(x // 2 for x in c))]; wgt = 1.0 if flat else 1 / np.sqrt(8)
            put(Kr, Kc, Kv, i, p, wgt * G[3] / 2j); put(Kr, Kc, Kv, p, i, -wgt * G[3] / 2j)
            put(Wr, Wc, Wv, i, p, -wgt * r / 2 * G[4]); put(Wr, Wc, Wv, p, i, -wgt * r / 2 * G[4])
    K = sp.csr_matrix((Kv, (Kr, Kc)), shape=(4 * Nc, 4 * Nc)); W = sp.csr_matrix((Wv, (Wr, Wc)), shape=(4 * Nc, 4 * Nc))
    depth = np.array([n for n, c in cells])
    return cells, index, K, W, depth
