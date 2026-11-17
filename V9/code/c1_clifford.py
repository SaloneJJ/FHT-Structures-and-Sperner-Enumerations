"""V62 / etape 1 : compression P c(x) P des operateurs de Clifford de sommet c(x) = eps_x + iota_x
sur l'homologie harmonique d'un representant de la classe de K26 (defaut : K26 vol. min., V37). JJ SALONE, FHT 2026"""
import sys, json, time, numpy as np
sys.path.insert(0, '.')
from simplicial import harmonic
src = sys.argv[1] if len(sys.argv) > 1 else '../../02_objets/V37_FHT_K26volmin.json'
F = [frozenset(f) for f in json.load(open(src))]
V = sorted(set(v for f in F for v in f)); order = {v: i for i, v in enumerate(V)}.get
t0 = time.time(); byd, ind, H = harmonic(F, order)
degs = sorted(d for d in H if H[d].shape[1])
print(src, ' faces', sum(len(byd[d]) for d in byd), ' Betti', {d: H[d].shape[1] for d in degs}, '%.0fs' % (time.time() - t0)); sys.stdout.flush()
def eps(x, d):
    """matrice de eps_x : C_d -> C_{d+1}"""
    M = np.zeros((len(byd[d+1]), len(byd[d])))
    for j, s in enumerate(byd[d]):
        if x in s: continue
        t = tuple(sorted(s + (x,), key=order))
        i = ind[d+1].get(t)
        if i is not None: M[i, j] = (-1) ** t.index(x)
    return M
off = {}; n = 0
for d in degs: off[d] = n; n += H[d].shape[1]
def compressed(x):
    C = np.zeros((n, n))
    for d in degs:
        if d + 1 in off:
            E = H[d+1].T @ eps(x, d) @ H[d]
            C[off[d+1]:off[d+1]+E.shape[0], off[d]:off[d]+E.shape[1]] = E
            C[off[d]:off[d]+E.shape[1], off[d+1]:off[d+1]+E.shape[0]] = E.T
    return C
res = {}
for x in V:
    C = compressed(x); res[x] = C
    blocks = {f'{d}->{d+1}': round(float(np.linalg.norm(C[off[d+1]:off[d+1]+H[d+1].shape[1], off[d]:off[d]+H[d].shape[1]])), 4)
              for d in degs if d + 1 in off}
    sv = np.linalg.svd(C, compute_uv=False)
    print(f'{x:6s} rang {int((sv > 1e-9).sum()):3d}  blocs {blocks}  sv {np.round(sv[:6], 4)}'); sys.stdout.flush()
np.save('c1_compressed.npy', np.array([res[x] for x in V])); json.dump(V, open('c1_vertices.json', 'w'))
