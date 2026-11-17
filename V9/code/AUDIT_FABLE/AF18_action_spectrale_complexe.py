#!/usr/bin/env python3
"""AF18 -- Action spectrale du complexe lui-meme : S_Lambda(X) = Tr f(D_X^2 / Lambda^2), D_X = d + d* (Hodge-de Rham) sur les chaines
reelles REDUITES de X (simplexe vide compris, de sorte que ker D_X = homologie reduite = 96 modes pour la classe K26).
Criteres ecrits AVANT calcul : 05_audit/ACTION_SPECTRALE_COMPLEXE_14-11.md (14/11/2026, 14 h 55).
Methode exacte : le spectre non nul de D_X^2 = Laplacien de Hodge total est la reunion, chaque valeur comptee DEUX fois, des carres des
valeurs singulieres non nulles des bords d_k : C_k -> C_{k-1} ; on diagonalise la plus petite des deux matrices de Gram de chaque d_k.
   Tr f(D^2/L^2) = f(0) * b + 2 * sum_k sum_{s>0} f(s^2 / L^2),   b = dim ker = somme des Betti reduits.
Usage :
   python3 AF18_action_spectrale_complexe.py spectre [noms ...]      noms : K4 V36 V37 V32 V129-011 V129-017 V129-018, ou fichiers .json
                                                                        (defaut : les 7 representants) -> AF18_spectre_<nom>.npz
   python3 AF18_action_spectrale_complexe.py bilan                     lit tous les AF18_spectre_*.npz -> AF18_bilan.txt
Memoire : K^(4) demande ~2 Go (matrice de Gram 9 176 x 9 176 et 11 996 x 11 996 en float64).  scipy utilise s'il est present.
"""
import os, sys, json, time, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
CAT = os.path.join(HERE, '..', '..', 'catalogues')
try:
    import scipy.sparse as sp
except Exception:
    sp = None

def load_named(nm):
    import AF15_paysage_puits as P15
    if nm in ('K4', 'V36', 'V37', 'V32'): return nm, [P15.names(f) for f in P15.load(nm)]
    if nm.startswith('V129-'):
        for l in open(os.path.join(CAT, 'V129_candidats_facettes.jsonl')):
            d = json.loads(l)
            if d['id'] == nm: return nm, d['facets']
    d = json.load(open(nm))
    d = d.get('facettes', d.get('facets')) if isinstance(d, dict) else d
    return os.path.basename(nm).replace('.json', ''), d

def chains(F):
    Vs = sorted(set(v for f in F for v in f)); ix = {v: i for i, v in enumerate(Vs)}
    S = set()
    for f in F:
        f = tuple(sorted(ix[v] for v in f))
        for k in range(1, len(f) + 1): S.update(itertools.combinations(f, k))
    byd = {-1: [()]}
    for s in S: byd.setdefault(len(s) - 1, []).append(s)
    for k in byd: byd[k].sort()
    return byd

def boundary_gram_eigs(byd, k):
    """valeurs propres non nulles de d_k^T d_k (= celles de d_k d_k^T), d_k : C_k -> C_{k-1}"""
    rows = {s: i for i, s in enumerate(byd[k - 1])}; cols = byd[k]
    I, J, X = [], [], []
    for j, s in enumerate(cols):
        for i in range(len(s)):
            I.append(rows[s[:i] + s[i + 1:]]); J.append(j); X.append(-1.0 if i % 2 else 1.0)
    m, nn = len(rows), len(cols)
    if sp is not None:
        Dm = sp.csr_matrix((X, (I, J)), shape=(m, nn))
        G = (Dm.T @ Dm if nn <= m else Dm @ Dm.T).toarray()
    else:
        Dm = np.zeros((m, nn)); Dm[I, J] = X
        G = Dm.T @ Dm if nn <= m else Dm @ Dm.T
        del Dm
    w = np.linalg.eigvalsh(G); del G
    tol = 1e-8 * max(1.0, w.max() if len(w) else 1.0)
    return w[w > tol], (m, nn)

def spectre(nm):
    nm, F = load_named(nm); t0 = time.time()
    byd = chains(F); N = sum(len(v) for v in byd.values())
    out = {}; nonzero = 0
    for k in range(0, max(byd) + 1):
        t = time.time(); w, sh = boundary_gram_eigs(byd, k); out[f'd{k}'] = w; nonzero += len(w)
        print(f'   {nm} : d_{k} {sh} rang {len(w)} ({time.time()-t:.0f} s)', flush=True)
    b = N - 2 * nonzero
    fv = {k: len(v) for k, v in byd.items()}
    np.savez(os.path.join(HERE, f'AF18_spectre_{nm}.npz'), betti_total=b, N=N, fvec=json.dumps(fv), **out)
    print(f'{nm} : N = {N} chaines (vide compris), dim ker D = {b} (attendu 96), {nonzero} valeurs singulieres non nulles ({time.time()-t0:.0f} s)', flush=True)

FUNCS = {'chaleur e^-x': lambda x: np.exp(-x), 'coupure 1_{x<=1}': lambda x: (x <= 1).astype(float), 'douce (1-x)_+^2': lambda x: np.clip(1 - x, 0, None) ** 2}

def bilan():
    sp_ = {}
    for f in sorted(os.listdir(HERE)):
        if f.startswith('AF18_spectre_') and f.endswith('.npz'):
            d = np.load(os.path.join(HERE, f)); nm = f[len('AF18_spectre_'):-4]
            s2 = np.concatenate([d[k] for k in d.files if k.startswith('d')])
            sp_[nm] = (int(d['betti_total']), int(d['N']), np.sort(s2))
    L2 = np.logspace(-2, 3, 51)
    lines = []
    def P(*a): s = ' '.join(str(x) for x in a); print(s); lines.append(s)
    P('complexe : dim ker, N, plus petite valeur propre non nulle (trou), multiplicite (x2), plus grande')
    for nm, (b, N, s2) in sp_.items():
        P(f'   {nm:12s} {b:4d} {N:7d}  trou {s2[0]:.6f} (x{2*int(np.sum(np.isclose(s2, s2[0], rtol=1e-6)))})  max {s2[-1]:.3f}  Tr D^2 = {2*s2.sum():.0f}')
    def excess(fn, s2, l2):
        """log de (S - f(0) b) / 2 pour la chaleur (pas de sous-depassement) ; valeur exacte sinon"""
        if fn.startswith('chaleur'):
            x = -s2 / l2; m = x.max(); return m + np.log(np.exp(x - m).sum())
        return FUNCS[fn](s2 / l2).sum()
    for fn in FUNCS:
        P(f'\n=== f = {fn} : S_Lambda - f(0) b = 2 sum f(s^2/Lambda^2) (log pour la chaleur) ; argmin par Lambda^2 (egalite exacte signalee)')
        prev = None
        for l2 in L2:
            E = {nm: excess(fn, s2, l2) for nm, (b, N, s2) in sp_.items()}
            ordre = sorted(E, key=E.get); m = E[ordre[0]]
            ex = [k for k in ordre if abs(E[k] - m) <= 1e-9 * max(1.0, abs(m))]
            best = ex[0] if len(ex) == 1 else 'EGALITE(' + ','.join(ex) + ')'
            if best != prev or l2 in (L2[0], L2[-1]):
                P(f'   Lambda^2 = {l2:9.3g} : argmin {best:24s} ordre {ordre}  ' + ', '.join(f'{k} {E[k]:.4g}' for k in ordre))
            prev = best
    open(os.path.join(HERE, 'AF18_bilan.txt'), 'w').write('\n'.join(lines) + '\n')

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: print(__doc__); sys.exit()
    if a[0] == 'spectre':
        for nm in (a[1:] or ['V37', 'V36', 'V129-017', 'V129-018', 'V32', 'V129-011', 'K4']): spectre(nm)
    elif a[0] == 'bilan': bilan()
