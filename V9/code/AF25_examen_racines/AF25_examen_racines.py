#!/usr/bin/env python3
"""AF25 -- examen des deux racines de A3 sur S^3 (tesseracte / oreiller). Criteres et predictions : 02_theoreme_holonomies/EXAMEN_RACINES.md §2
(ecrits avant calcul).  (1) symetries dans SO(4) via 2O x 2O ; (2) quotients S^3/Q8 et RP^3 comme un cube recolle (theoreme H2') ;
(3) Laplacien combinatoire du graphe des cubes (niveau fin) : multiplicites et paquets.
Usage : python3 AF25_examen_racines.py | tee AF25.log      (quelques minutes ; numpy + scipy)
"""
import sys, os, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'AF24_holonomies')); sys.path.insert(0, HERE)
import AF24_theoreme_holonomies as T

# ---------------------------------------------------------------- quaternions, 2O
def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return np.array([a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2, a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0])
conj = lambda a: a * np.array([1, -1, -1, -1])
def group_2O():
    S = set()
    for i in range(4):
        for s in (1, -1): v = [0.] * 4; v[i] = s; S.add(tuple(v))
    for s in itertools.product((.5, -.5), repeat=4): S.add(s)
    r = 1 / np.sqrt(2)
    for i, j in itertools.combinations(range(4), 2):
        for a in (r, -r):
            for b in (r, -r): v = [0.] * 4; v[i] = a; v[j] = b; S.add(tuple(np.round(v, 12)))
    return [np.array(x) for x in S]
G = group_2O(); assert len(G) == 48
Q8 = [g for g in G if np.isclose(np.max(np.abs(g)), 1)]; assert len(Q8) == 8
def mat(l, r):   # x -> l x conj(r)
    return np.array([qmul(qmul(l, e), conj(r)) for e in np.eye(4)]).T
def key(M): return tuple(np.round(M, 8).ravel())
def inQ8(q): return any(np.allclose(q, p) for p in Q8)

def symetries():
    E = np.eye(4); tess, ore = set(), set(); tess_pairs = 0; cond_ok = True
    for l in G:
        for r in G:
            M = mat(l, r)
            # tesseracte : envoie Q8 (centres des cubes) sur Q8
            if all(any(np.allclose(M @ p, q) for q in Q8) for p in Q8):
                tess.add(key(M)); tess_pairs += 1
                cond_ok &= inQ8(qmul(l, conj(r)))
            # oreiller : e0 -> +-e0 (hemispheres) et {+-i, +-j, +-k} -> lui-meme (cube equatorial)
            if abs(abs(M[0, 0]) - 1) < 1e-9 and all(any(np.allclose(M @ p, q) for q in Q8 if abs(q[0]) < .5) for p in Q8 if abs(p[0]) < .5):
                ore.add(key(M))
    signed = {key(np.array(P) * np.array(s)[:, None]) for P in itertools.permutations(np.eye(4)) for s in itertools.product((1, -1), repeat=4)
              if np.isclose(np.linalg.det(np.array(P) * np.array(s)[:, None]), 1)}
    print(f'(1) tesseracte : {tess_pairs} paires (l, r) de 2O x 2O, {len(tess)} rotations ; = toutes les permutations signees de det +1 ({len(signed)}) : {tess == signed} ;'
          f' condition l r^-1 dans Q8 : {cond_ok}')
    print(f'    oreiller : {len(ore)} rotations ; forme (l, +-l) :',
          all(any(np.allclose(mat(l, l) * s, np.array(k).reshape(4, 4)) for l in G for s in (1, -1)) for k in ore))
    # translations a gauche par Q8 : symetries du tesseracte ; a gauche par i : symetrie de l'oreiller ?
    i_ = np.array([0, 1, 0, 0.])
    print('    translation a gauche par i : symetrie du tesseracte', key(mat(i_, np.array([1, 0, 0, 0.]))) in tess,
          '; de l\'oreiller', key(mat(i_, np.array([1, 0, 0, 0.]))) in ore)
    print('    2O/Q8 : ordre', len(G) // len(Q8), '; classes l mod Q8 =', len({tuple(sorted(key(np.array([qmul(l, q)])) for q in Q8)) for l in G}))

# ---------------------------------------------------------------- quotients : un cube recolle
def cube_quotient(m, kind):
    """[0, m]^3 subdivise ; faces opposees recollees par quart de tour (Q8) ou demi-tour (RP^3), autour de l'axe normal, centre du cube"""
    c = m / 2
    def rot(v, ax, k):   # rotation de k quarts de tour autour de l'axe ax passant par le centre
        v = list(v); a, b = (ax + 1) % 3, (ax + 2) % 3
        for _ in range(k % 4): v[a], v[b] = c - (v[b] - c), v[a]   # (a, b) -> (-b, a) recentres
        return tuple(v)
    k = 1 if kind == 'Q8' else 2
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    pts = list(itertools.product(range(m + 1), repeat=3))
    for p in pts:
        for ax in range(3):
            if p[ax] == m:          # face x_ax = m  ~  face x_ax = 0 par rotation
                q = list(p); q[ax] = 0; q = rot(q, ax, k)
                q = tuple(int(round(t)) for t in q); parent[find(p)] = find(q)
    cubes = [{cc: find(tuple(co[i] + cc[i] for i in range(3))) for cc in itertools.product((0, 1), repeat=3)} for co in itertools.product(range(m), repeat=3)]
    return T.Cx(f'{"S^3/Q8 (quarts de tour)" if kind == "Q8" else "RP^3 (demi-tours)"}, un cube subdivise {m}^3', cubes)

# ---------------------------------------------------------------- spectres
def spectre(C, k=40, label=''):
    import scipy.sparse as sp, scipy.sparse.linalg as sla
    n = len(C.cubes); I, J = [], []
    for f, lst in C.fcubes.items():
        (a, _, _), (b, _, _) = lst; I += [a, b]; J += [b, a]
    A = sp.csr_matrix((np.ones(len(I)), (I, J)), shape=(n, n)); L = sp.diags(np.asarray(A.sum(1)).ravel()) - A
    if n <= 3000: E = np.linalg.eigvalsh(L.toarray())[:k]
    else: E = np.sort(sla.eigsh(L.tocsc(), k=k, sigma=-1e-3, which='LM', return_eigenvectors=False))
    E = E[1:]; mult = []; cur = [E[0]]
    for x in E[1:]:
        if abs(x - cur[0]) < 1e-8 * max(1, cur[0]): cur.append(x)
        else: mult.append(cur); cur = [x]
    mult.append(cur)
    ex = [(len(c), round(c[0] / E[0], 4)) for c in mult][:10]
    pk = []; cur = [E[0]]
    for x in E[1:]:
        if (x - cur[0]) / cur[0] < 0.05: cur.append(x)
        else: pk.append(cur); cur = [x]
    pk.append(cur)
    print(f'(3) {label} ({n} cubes) : multiplicites exactes (taille, lambda/lambda_1) : {ex}')
    print(f'    paquets a 5 % (taille, moyenne/lambda_1) : {[(len(p), round(float(np.mean(p)) / E[0], 3)) for p in pk[:6]]}   [S^3 : (4, 1), (9, 2,667), (16, 5)]')

if __name__ == '__main__':
    print(__doc__.strip().splitlines()[0])
    symetries()
    for kind in ('Q8', 'RP3'):
        for m in (3, 5):
            T.analyse(cube_quotient(m, kind))
    print()
    for n in (2, 3): spectre(T.bord_4cube(n), label=f'tesseracte niveau {n}')
    for n in (3, 4): spectre(T.oreiller(n), label=f'oreiller niveau {n}')
