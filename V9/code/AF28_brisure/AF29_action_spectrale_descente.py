#!/usr/bin/env python3
"""AF29 / D18 -- action spectrale S_Lambda = Tr f(D_X^2/Lambda^2) le long de la descente (D_X = d + d* sur les chaines reelles reduites).
Criteres : 08_descente_Xgrand/CRITERES_D18_action_spectrale_descente.md (ecrits avant calcul).
Mesures : N, Tr D^2, ker, spectre bas de D^2 (par degre, solveur creux), trou et multiplicite ; S_Lambda (chaleur) aux petits Lambda, avec borne de troncature.
Usage : python3 AF29_action_spectrale_descente.py [noms] | tee AF29.log      (noms parmi K26 (= K^(4), definition AF3) K26_V25 X_comp X_chi X_SM "X_LR'" "X_PS'" X_grand ; sans nom : tous)
Sorties : AF29_<nom>.npz (valeurs propres basses par degre) et AF29.log.  Memoire : quelques Go pour les plus gros complexes."""
import sys, os, json, itertools, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
H_ = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(H_, 'D17_Dirac_fini.py')).read()
exec(src.split("def aut_sgs_base")[0].split("import AF9_porteurs_triplet as AF9")[0])
sys.path.insert(0, os.path.join(H_, '..', 'AF20_algebres_tordues')); from fht_tools import maximal
exec("d = json.load" + src.split("d = json.load")[1].split("for name in (sys.argv[1:]")[0])
NLOW = 200
sys.path.insert(0, os.path.join(H_, '..', 'AF9_test_triplet'))
from AF3_K26_homologie_automorphismes import K26_definition
cx['K26_V25'] = cx['K26']; cx['K26'] = lambda: [sorted(f) for f in K26_definition()]
def faces(X):
    S = set()
    for f in X:
        f = tuple(sorted(f))
        for r in range(1, len(f) + 1): S.update(itertools.combinations(f, r))
    byd = {}
    for s in S: byd.setdefault(len(s) - 1, []).append(s)
    for k in byd: byd[k].sort()
    return byd
def boundaries(byd):
    idx = {k: {s: i for i, s in enumerate(L)} for k, L in byd.items()}; B = {}
    B[0] = sp.csr_matrix(np.ones((1, len(byd[0]))))                     # augmentation (complexe reduit)
    for k in range(1, max(byd) + 1):
        r, c, v = [], [], []
        for j, s in enumerate(byd[k]):
            for t in range(len(s)): r.append(idx[k - 1][s[:t] + s[t + 1:]]); c.append(j); v.append((-1) ** t)
        B[k] = sp.csr_matrix((np.array(v, dtype=float), (r, c)), shape=(len(byd[k - 1]), len(byd[k])))
    return B
def analyse(name, X):
    t = time.time(); byd = faces(X); B = boundaries(byd); top = max(byd)
    N = sum(len(L) for L in byd.values()) + 1
    trD2 = 2 * sum((k + 1) * len(L) for k, L in byd.items())
    ker = 0; low = []; per = {}
    for k in range(-1, top + 1):
        if k == -1:
            Lk = (B[0] @ B[0].T).tocsc().astype(float)                                 # degre -1 (simplexe vide)
        else:
            Lk = B[k].T @ B[k] + (B[k + 1] @ B[k + 1].T if k + 1 <= top else 0)
            Lk = sp.csc_matrix(Lk, dtype=float)
        n = Lk.shape[0]
        if n <= 2500:
            w = np.linalg.eigvalsh(Lk.toarray())
        else:
            kk = min(n - 2, NLOW + 200)
            w = sla.eigsh(Lk, k=kk, sigma=-1e-3, which='LM', return_eigenvectors=False)
        w = np.sort(np.abs(w)); z = int(np.sum(w < 1e-8)); ker += z
        nz = w[w >= 1e-8]; per[k] = nz[:NLOW]; low.extend(nz[:NLOW].tolist())
    low = np.sort(np.array(low)); lam1 = low[0]; mult = int(np.sum(np.abs(low - lam1) < 1e-6))
    np.savez(os.path.join(H_, f'AF29_{name.replace(chr(39), "p")}.npz'), **{f'deg{k}': v for k, v in per.items()})
    print(f'== {name} : N = {N}, Tr D^2 = {trD2}, ker = {ker}, trou = {lam1:.6f} (x{mult}), 10 plus basses : {np.round(low[:10], 5).tolist()} ({time.time()-t:.0f} s)', flush=True)
    lmax_calc = min(v[-1] for v in per.values() if len(v) == NLOW) if any(len(v) == NLOW for v in per.values()) else np.inf
    out = {}
    for r in (0.05, 0.1, 0.2, 0.5):
        L2 = r * lam1; S = ker + float(np.sum(np.exp(-low / L2)))
        bound = (N - ker - len(low)) * np.exp(-lmax_calc / L2)
        out[r] = (S, bound)
    print('   S_Lambda (chaleur) - 96 aux Lambda^2 = r * trou :', {r: f'{S - 96:.4e} (borne {b:.1e})' for r, (S, b) in out.items()}, flush=True)
    return dict(N=N, trD2=trD2, ker=ker, lam1=lam1, mult=mult, low=low)
res = {}
for name in (sys.argv[1:] or ['K26', 'K26_V25', 'X_comp', 'X_chi', 'X_SM', "X_LR'", "X_PS'", 'X_grand']):
    res[name] = analyse(name, cx[name]())
if len(res) > 1:
    print('\n=== bilan basse energie (chaleur, S - 96 = sum exp(-lambda/Lambda^2) sur le spectre bas)')
    lmin = min(r['lam1'] for r in res.values())
    for r in (0.05, 0.1, 0.2, 0.5, 1.0):
        L2 = r * lmin; vals = {n: float(np.sum(np.exp(-v['low'] / L2))) for n, v in res.items()}
        print(f'   Lambda^2 = {L2:.4f} : ordre {sorted(vals, key=vals.get)}  ', {n: f'{x:.3e}' for n, x in vals.items()})
    print('   grand Lambda (volume N) : ordre', sorted(res, key=lambda n: res[n]['N']), {n: res[n]['N'] for n in res})
