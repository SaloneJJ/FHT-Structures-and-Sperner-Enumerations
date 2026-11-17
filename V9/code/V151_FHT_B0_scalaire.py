"""V151 B0 : champ scalaire massif sur l'arbre augmente ferme en T^3 ; fonction a deux points au bord.
Criteres et predictions : V151 par. 1 (ecrits avant calcul).   Usage : python3 V151_FHT_B0_scalaire.py
JJ SALONE, FHT 2026"""
import time, itertools
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
t0 = time.time()

def grid(s):
    I, J, K = np.meshgrid(np.arange(s), np.arange(s), np.arange(s), indexing='ij'); return I.ravel(), J.ravel(), K.ravel()

def laplacian(N, kind):
    sides = [2 ** n for n in range(N + 1)]; offs = np.concatenate([[0], np.cumsum([s ** 3 for s in sides])]); Nc = offs[-1]
    A_, B_ = [], []
    for n, s in enumerate(sides):
        I, J, K = grid(s); base = offs[n] + (I * s + J) * s + K
        if n >= 1:
            sp_ = sides[n - 1]; A_.append(base); B_.append(offs[n - 1] + ((I // 2) * sp_ + J // 2) * sp_ + K // 2)
        if kind == 'nu': continue
        dirs = [d for d in itertools.product((-1, 0, 1), repeat=3) if any(d)]
        if kind == 'faces': dirs = [d for d in dirs if sum(map(abs, d)) == 1]
        for d in dirs:
            A_.append(base); B_.append(offs[n] + (((I + d[0]) % s) * s + (J + d[1]) % s) * s + (K + d[2]) % s)
    a = np.concatenate(A_); b = np.concatenate(B_); keep = a != b
    a, b = np.minimum(a[keep], b[keep]), np.maximum(a[keep], b[keep])
    key = np.unique(a * Nc + b); a, b = key // Nc, key % Nc
    Adj = sp.coo_matrix((np.ones(len(a)), (a, b)), shape=(Nc, Nc)); Adj = (Adj + Adj.T).tocsr()
    L = sp.diags(np.asarray(Adj.sum(1)).ravel()) - Adj
    return L.tocsr(), offs, Nc

def green_top(L, offs, Nc, N, mu2, src):
    s = 2 ** N; M = (L + mu2 * sp.identity(Nc)).tocsr()
    b = np.zeros(Nc); b[offs[N] + (src[0] * s + src[1]) * s + src[2]] = 1.0
    x, info = sla.cg(M, b, rtol=1e-12, maxiter=20000); assert info == 0
    return x[offs[N]:].reshape(s, s, s)

def images(v):
    out = set()
    for p in itertools.permutations(range(3)):
        for sg in itertools.product((1, -1), repeat=3):
            out.add(tuple(sg[i] * v[p[i]] for i in range(3)))
    return sorted(out)

families = {9: [(3, 0, 0), (2, 2, 1)], 25: [(5, 0, 0), (4, 3, 0)], 49: [(7, 0, 0), (6, 3, 2)], 81: [(9, 0, 0), (8, 4, 1), (7, 4, 4), (6, 6, 3)]}
for N in (4, 5, 6):
    s = 2 ** N
    srcs = [(0, 0, 0), (s // 2 - 1, s // 2 - 1, s // 2 - 1), (s // 4 + 1, s // 2 - 1, 3 % s), (5 % s, 3 % s, 9 % s)]
    for kind in ('nu', 'faces', 'roi'):
        L, offs, Nc = laplacian(N, kind)
        for mu2 in (1.0, 0.25):
            Gs = {src: green_top(L, offs, Nc, N, mu2, src) for src in srcs}
            val = lambda src, v: Gs[src][(src[0] + v[0]) % s, (src[1] + v[1]) % s, (src[2] + v[2]) % s]
            # inhomogeneite : meme vecteur, sources differentes (et images cubiques)
            inh = []
            for v in [(1, 0, 0), (2, 0, 0), (3, 0, 0), (1, 1, 0), (2, 2, 1)]:
                if max(map(abs, v)) >= s // 2: continue
                vals = [val(src, w) for src in srcs for w in images(v)]; inh.append(max(vals) / min(vals))
            # anisotropie : vecteurs de meme norme, moyennes sur sources et images
            ani = {}
            for n2, fam in families.items():
                if max(max(map(abs, f)) for f in fam) >= s // 2: continue
                means = [np.mean([val(src, w) for src in srcs for w in images(f)]) for f in fam]
                ani[n2] = max(means) / min(means)
            # decroissance radiale (moyenne sur sources et directions), residus en log2 r
            fr = np.fft.fftfreq(s, d=1.0 / s).astype(int)
            R2 = fr[:, None, None] ** 2 + fr[None, :, None] ** 2 + fr[None, None, :] ** 2
            acc = {}
            for src in srcs:
                Gr = np.roll(np.roll(np.roll(Gs[src], -src[0], 0), -src[1], 1), -src[2], 2)
                for q in np.unique(R2):
                    if 4 <= q <= (s // 4) ** 2: acc.setdefault(q, []).append(Gr[R2 == q].mean())
            q = np.array(sorted(acc)); g = np.array([np.mean(acc[x]) for x in q]); lr = 0.5 * np.log2(q); lg = np.log2(g)
            if len(q) > 3:
                pente, b0 = np.polyfit(lr, lg, 1); res = lg - (pente * lr + b0); amp = float(np.ptp(res))
            else:
                pente, amp = np.nan, np.nan
            print(f'N={N} {kind:5s} mu2={mu2:4.2f} : inhomogeneite max/min = {", ".join(f"{x:.2f}" for x in inh)} ; '
                  f'anisotropie |r|^2->rapport {", ".join(f"{k}:{v:.3f}" for k, v in ani.items())} ; pente log2 G / log2 r = {pente:.3f} ; amplitude des residus (log2) = {amp:.3f}', flush=True)
print(f'temps {time.time()-t0:.0f} s')
