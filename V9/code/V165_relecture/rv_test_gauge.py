"""Point 1 : coherence de la construction du flux.
(i) holonomies de boucles aleatoires fermees (laterales + verticales + enroulements) contre exp(i B aire projetee) ;
(ii) lignes de Wilson non contractiles : jauge G1 contre jauge G2 (doivent coincider) ;
(iii) spectres G1 / G2 et H(2 pi) / H(0) ; (iv) comparaison exacte avec les matrices de V158 (meme jauge G1).
Usage : python3 rv_test_gauge.py N"""
import sys, os, numpy as np, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_core import *
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
rng = np.random.default_rng(7)

def adjacency(g, th):
    T, S = g.ops(th)
    S = S / g.wv
    nb = {i: [] for i in range(g.Nc)}       # (voisin, U(voisin <- i), deplacement dans le revetement)
    for j, Tj in enumerate(T):
        Tc = Tj.tocoo()
        for t, s_, u in zip(Tc.row, Tc.col, Tc.data):
            d = np.zeros(3); d[j] = g.dsz[t]
            nb[s_].append((t, u, -d))            # s -> t : t = s - d e_j
            nb[t].append((s_, np.conj(u), d))
    Sc = S.tocoo()
    for c, p, u in zip(Sc.row, Sc.col, Sc.data):
        d = g.cen[c] - g.cen[p]
        nb[p].append((c, u, d)); nb[c].append((p, np.conj(u), -d))
    return nb

for flat in (False, True):
    for gauge in ('G1', 'G2'):
        g = Geo(N, flat, gauge); nb = adjacency(g, 0.7)
        err_c, ncon, wil = 0.0, 0, {}
        from collections import deque
        for trial in range(400):
            i0 = rng.integers(g.Nc); i = i0; pos = g.cen[i0].copy(); path = [pos.copy()]; hol = 1.0 + 0j
            for step in range(rng.integers(2, 25)):
                k, u, d = nb[i][rng.integers(len(nb[i]))]
                hol *= u; pos = pos + d; path.append(pos.copy()); i = k
            # retour par plus court chemin (BFS) vers i0
            prev = {i: None}; dq = deque([i])
            while dq:
                a = dq.popleft()
                if a == i0: break
                for (k, u, d) in nb[a]:
                    if k not in prev: prev[k] = (a, u, d); dq.append(k)
            seq = []; a = i0
            while prev[a] is not None: seq.append((a,) + prev[a]); a = prev[a][0]
            for (k, a_, u, d) in reversed(seq):
                hol *= u; pos = pos + d; path.append(pos.copy())
            i = i0
            # refermer : retour par le meme chemin est trivial ; on ne garde que les boucles revenues sur i0
            if i != i0: continue
            P = np.array(path); disp = np.round((P[-1] - P[0]) / g.L).astype(int)
            area = 0.5 * np.sum(P[:-1, 0] * P[1:, 1] - P[1:, 0] * P[:-1, 1])
            th_term = 0.7 * (P[-1, 2] - P[0, 2]) / g.L
            if np.all(disp == 0):
                ncon += 1; err_c = max(err_c, abs(np.angle(hol * np.exp(-1j * (g.B * area + th_term)))))
        print(f'{"plat " if flat else "arbre"} {gauge} : {ncon} boucles contractiles fermees, ecart max a exp(i B aire) = {err_c:.1e}')

# lignes de Wilson : chemin droit le long de x (resp. y, z) a chaque niveau, G1 vs G2
def wilson_lines(g, th):
    T, S = g.ops(th); out = []
    for n, sd in enumerate(g.sides):
        for j in range(3):
            for x0 in itertools.product(range(sd), repeat=2):
                x = [0, 0, 0]; others = [k for k in range(3) if k != j]
                x[others[0]], x[others[1]] = x0
                hol = 1.0 + 0j
                for step in range(sd):
                    t = g.index[(n, tuple(x))]; y = list(x); y[j] = (y[j] + 1) % sd; s_ = g.index[(n, tuple(y))]
                    hol *= T[j][t, s_]; x = y
                out.append(hol)
    return np.array(out)
import itertools
for flat in (False, True):
    w1 = wilson_lines(Geo(N, flat, 'G1'), 0.7); w2 = wilson_lines(Geo(N, flat, 'G2'), 0.7)
    print(f'{"plat " if flat else "arbre"} : lignes de Wilson G1 vs G2, ecart max {np.abs(w1 - w2).max():.1e} ({len(w1)} lignes)')

# spectres
for flat in (False,):
    g1, g2 = Geo(N, flat, 'G1'), Geo(N, flat, 'G2')
    for name, Hf in (('Wilson m=-1', lambda g, th: H_wilson(g, th)), ('KS + eps (O)', lambda g, th: H_ks_kin(g, th) + sp.diags(mass_profile(g, 1.0, False) * ks_phases(g)[2]))):
        e1 = np.linalg.eigvalsh(Hf(g1, 0.7).toarray()); e2 = np.linalg.eigvalsh(Hf(g2, 0.7).toarray())
        e0 = np.linalg.eigvalsh(Hf(g1, 0.0).toarray()); e2p = np.linalg.eigvalsh(Hf(g1, 2 * np.pi).toarray())
        print(f'arbre {name} : spectre G1 vs G2 {np.abs(e1 - e2).max():.1e} ; H(2pi) vs H(0) {np.abs(e0 - e2p).max():.1e}')

# comparaison avec V158 (meme jauge) : on charge la partie 'geometrie + modeles' du script original sans l'executer en entier
src = open('/home/claude/fht/05_code/code_v30/V158_FHT_R2_staggered.py').read().split('# ---------------------------------------------------------------- controles')[0]
for geom in ('tree', 'flat'):
    for mod in ('bor', 'ksO'):
        ns = {}; sys.argv = ['x', mod, geom, str(N), '-1.0' if mod == 'bor' else '1.0', '8']
        exec(src, ns)
        g = Geo(N, geom == 'flat', 'G1')
        # V158 ordonne les cellules de la meme facon (niveau, x, y, z) : verification
        assert np.all(ns['coords'] == g.coords) and np.all(ns['depth'] == g.depth)
        th = 1.1
        if mod == 'bor': Hm, Hv = H_wilson(g, th), ns['H_of'](th)
        else: Hm, Hv = H_ks_kin(g, th) + sp.diags(np.full(g.Nc, 1.0) * ks_phases(g)[2]), ns['H_of'](th)
        Dm = (H_wilson(g, th, deriv=True) if mod == 'bor' else H_ks_kin(g, th, deriv=True)); Dv = ns['H_of'](th, True)
        print(f'V158 {mod} {geom} : |H_relecture - H_V158| = {abs(Hm - Hv).max():.1e} ; |dH| = {abs(Dm - Dv).max():.1e}')
