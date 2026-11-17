"""V62 / etape 1 : second ordre de Feshbach, T_x = P c(x) D^+ c(x) P (impair, non local), D = d + d*, D^+ pseudo-inverse.
JJ SALONE, FHT 2026"""
import sys, json, time, numpy as np
sys.path.insert(0, '.')
from simplicial import harmonic
src = sys.argv[1] if len(sys.argv) > 1 else '../../02_objets/V37_FHT_K26volmin.json'
F = [frozenset(f) for f in json.load(open(src))]
V = sorted(set(v for f in F for v in f)); order = {v: i for i, v in enumerate(V)}.get
byd, ind, H = harmonic(F, order)
ds = sorted(byd); sz = {d: len(byd[d]) for d in ds}
off = {}; n = 0
for d in ds: off[d] = n; n += sz[d]
def sl(d): return slice(off[d], off[d] + sz[d])
Dm = np.zeros((n, n))
for d in ds:
    if d - 1 in off:
        B = np.zeros((sz[d-1], sz[d]))
        for j, s in enumerate(byd[d]):
            for t in range(len(s)): B[ind[d-1][s[:t] + s[t+1:]], j] = (-1) ** t
        Dm[sl(d-1), sl(d)] = B; Dm[sl(d), sl(d-1)] = B.T
# pseudo-inverse par blocs de degre : D^+ = D (D^2)^+, D^2 diagonal par degre
L = Dm @ Dm
Lp = np.zeros((n, n))
for d in ds:
    w, U = np.linalg.eigh(L[sl(d), sl(d)]); m = w > 1e-8
    Lp[sl(d), sl(d)] = (U[:, m] / w[m]) @ U[:, m].T
Dp = Dm @ Lp
Hd = [d for d in ds if H[d].shape[1]]
P = np.zeros((n, sum(H[d].shape[1] for d in Hd))); k = 0; hoff = {}
for d in Hd:
    hoff[d] = k; P[sl(d), k:k + H[d].shape[1]] = H[d]; k += H[d].shape[1]
print('dim chaines', n, ' harmonique', P.shape[1]); sys.stdout.flush()
def cliff(x):
    C = np.zeros((n, n))
    for d in ds:
        if d + 1 not in off: continue
        for j, s in enumerate(byd[d]):
            if x in s: continue
            t = tuple(sorted(s + (x,), key=order)); i = ind[d+1].get(t)
            if i is not None:
                C[off[d+1] + i, off[d] + j] = (-1) ** t.index(x); C[off[d] + j, off[d+1] + i] = (-1) ** t.index(x)
    return C
out = {}
pairs = [tuple(a.split(',')) for a in sys.argv[2:]] if len(sys.argv) > 2 else [(x, x) for x in V]
for (x, y) in pairs:
    T = P.T @ cliff(x) @ Dp @ cliff(y) @ P; x = x + '|' + y
    blocks = {}
    for a in Hd:
        for b in Hd:
            nb = np.linalg.norm(T[hoff[b]:hoff[b]+H[b].shape[1], hoff[a]:hoff[a]+H[a].shape[1]])
            if nb > 1e-9: blocks[f'{a}->{b}'] = round(float(nb), 4)
    sv = np.linalg.svd(T, compute_uv=False)
    print(f'{x:12s} rang {int((sv > 1e-9).sum()):3d} blocs {blocks} sv {np.round(sv[:5], 4)}'); sys.stdout.flush()
    out[x] = T
np.save('c2_T.npy', np.array(list(out.values())))
