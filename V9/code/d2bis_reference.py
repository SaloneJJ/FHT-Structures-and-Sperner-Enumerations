#!/usr/bin/env python3
"""D2-bis -- reference : Wilson-Dirac sans masse (4 composantes, r = 1) au seul niveau fin, sur l'OREILLER de D2_oreiller.py (meme recollement, repere miroir de B,
transports V_f = s_f alpha_a alpha_b, classe de signes fixee par le meme principe : holonomie S(+pi) aux 12 aretes) et sur le TORE (8 cubes periodiques).
Criteres : CRITERES_D2bis.md (ecrits avant calcul). Usage : python3 d2bis_reference.py signes | N [N ...]   (sortie : niveaux |E| * 2^N et multiplicites)"""
import sys, itertools, json, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
AL = [np.kron(sx, s) for s in (sx, sy, sz)]; BE = np.kron(sz, I2); I4 = np.eye(4); R_W = 1.0

def V_face(d, sign=1):
    a, b = [c for c in range(3) if c != d]; return sign * AL[a] @ AL[b]

def neighbour(mode, m, q, co, d, e, signs):
    c2 = list(co); c2[d] += e
    if mode == 'tore':
        g = [q[i] * m + co[i] for i in range(3)]; g[d] = (g[d] + e) % (2 * m)
        return tuple(x // m for x in g), tuple(x % m for x in g), I4
    if 0 <= c2[d] < m: return q, tuple(c2), I4
    q2 = 'B' if q == 'A' else 'A'
    co2 = tuple((0 if e > 0 else m - 1) if i == d else m - 1 - co[i] for i in range(3))
    face = (d, e) if q == 'A' else (d, -e)
    V = V_face(d, signs[face]); W = V.conj().T if q == 'A' else V
    return q2, co2, W

def S_rot(a, b, th): return np.cos(th / 2) * I4 - np.sin(th / 2) * AL[a] @ AL[b]

def edge_holonomies(signs):
    m = 2; out = {}
    for d in range(3):
        a, b = [c for c in range(3) if c != d]
        for ea, eb in itertools.product((1, -1), repeat=2):
            co = [0, 0, 0]; co[a] = m - 1 if ea > 0 else 0; co[b] = m - 1 if eb > 0 else 0
            P = I4.copy(); q, c = 'A', tuple(co)
            q, c, W = neighbour('oreiller', m, q, c, a, ea, signs); P = W.conj().T @ P
            q, c, W = neighbour('oreiller', m, q, c, b, -eb, signs); P = W.conj().T @ P
            assert q == 'A'
            Sp = S_rot(a, b, ea * eb * np.pi)
            out[(d, ea, eb)] = +1 if np.allclose(P, Sp) else (-1 if np.allclose(P, -Sp) else 0)
    return out

def signes():
    faces = [(d, e) for d in range(3) for e in (1, -1)]; res = []
    for bits in itertools.product((1, -1), repeat=5):
        sg = dict(zip(faces, (1,) + bits))
        if all(v == +1 for v in edge_holonomies(sg).values()): res.append(sg)
    return res

def build(N, mode, signs):
    m = 2 ** N; cubes = ['A', 'B'] if mode == 'oreiller' else list(itertools.product((0, 1), repeat=3))
    cells = [(q, co) for q in cubes for co in itertools.product(range(m), repeat=3)]; idx = {c: i for i, c in enumerate(cells)}
    I, J, V = [], [], []
    def add(u, v, B):
        rr, cc = np.nonzero(np.abs(B) > 1e-14); I.extend(4 * u + rr); J.extend(4 * v + cc); V.extend(B[rr, cc])
    for (q, co) in cells:
        u = idx[(q, co)]; add(u, u, 3 * R_W * BE)
        for d in range(3):
            q2, co2, W = neighbour(mode, m, q, co, d, +1, signs); v = idx[(q2, co2)]
            B = (AL[d] / 2j - R_W / 2 * BE) @ W
            add(u, v, B); add(v, u, B.conj().T)
    return sp.csr_matrix((V, (I, J)), shape=(4 * len(cells), 4 * len(cells)))

def low(H, k):
    n = H.shape[0]
    if n <= 4000: E = np.linalg.eigvalsh(H.toarray())
    else: E = sla.eigsh(H.tocsc(), k=min(k, n - 2), sigma=1e-6, which='LM', return_eigenvectors=False)
    return np.sort(np.abs(E))

def packets(E, scale, tol=2e-3):
    out = []
    for x in E * scale:
        if out and abs(x - out[-1][0]) < tol * max(1, x): out[-1][1] += 1
        else: out.append([x, 1])
    return [(round(a, 4), c) for a, c in out]

if __name__ == '__main__':
    cls = signes()
    print('classes de signes satisfaisant le principe (S(+pi) aux 12 aretes) :', [{f'{k[0]}{"+" if k[1] > 0 else "-"}': v for k, v in c.items()} for c in cls], flush=True)
    if sys.argv[1:] == ['signes']: sys.exit()
    sg = cls[0]
    H = build(2, 'oreiller', sg); print('hermitien :', abs(H - H.conj().T).max() < 1e-12)
    for N in map(int, sys.argv[1:]):
        t0 = time.time(); sc = 2 ** N
        Et = low(build(N, 'tore', sg), 120); Eo = low(build(N, 'oreiller', sg), 80)
        print(f'N = {N} ({time.time() - t0:.0f} s)\n  tore     (E 2^N, multiplicite) : {packets(Et, sc)[:6]}\n  oreiller (E 2^N, multiplicite) : {packets(Eo, sc)[:6]}', flush=True)
