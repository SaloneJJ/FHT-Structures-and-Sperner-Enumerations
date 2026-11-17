"""V359 : recalcul indépendant de la torsion combinatoire (convention de l'article IV / V100 / V354 :
c = ½ Σ_{q ≥ −1} (−1)^q ln det′(d_q* d_q), complexe augmenté, poids unité), par SVD dense des cobords.
Contrôles : C3 (½ ln 3), simplexe (0), RP² à 6 sommets (torsion H¹ = 0, H² = ℤ/2 : c = ± ln 2).
Aussi : Σ(−1)^q ln det′Δ_q (superdéterminant, doit être 0) ; Betti mod 2, 3, 5, 7 (absence de p-torsion).
Usage : python3 V359_torsion.py <dossier JSON>   (numpy)."""
import json, sys, itertools as it, numpy as np
from collections import defaultdict
D = sys.argv[1] if len(sys.argv) > 1 else '.'
def load(f):
    d = json.load(open(f"{D}/{f}")); d = d['facettes'] if isinstance(d, dict) else d
    return [frozenset(x) for x in d]
def cells(F):
    S = set()
    for f in F:
        s = sorted(f)
        for k in range(1, len(s)+1): S.update(it.combinations(s, k))
    by = defaultdict(list)
    for c in S: by[len(c)-1].append(c)
    return {d: sorted(L) for d, L in by.items()}
def cobords(F):
    by = cells(F); top = max(by); idx = {d: {c: i for i, c in enumerate(by[d])} for d in by}
    M = {-1: np.ones((len(by[0]), 1))}
    for q in range(0, top):
        A = np.zeros((len(by[q+1]), len(by[q])))
        for j, c in enumerate(by[q+1]):
            for i in range(len(c)): A[j, idx[q][c[:i]+c[i+1:]]] = (-1)**i
        M[q] = A
    return M, by
def torsion(F, verbose=False):
    M, by = cobords(F); c = 0.0; lnd = {}; rk = {}
    for q, A in M.items():
        s = np.linalg.svd(A, compute_uv=False); s = s[s > 1e-9*max(1, s[0])]
        lnd[q] = 2*np.log(s).sum(); rk[q] = len(s); c += 0.5*(-1)**q*lnd[q]
    n = {d: len(L) for d, L in by.items()}; n[-1] = 1
    betti = [n[q] - rk.get(q, 0) - rk.get(q-1, 0) for q in range(0, max(by)+1)]
    # superdéterminant : ln det′Δ_q = lnd[q] + lnd[q−1] ; Σ(−1)^q (...) télescope
    sd = sum((-1)**q*(lnd.get(q, 0) + lnd.get(q-1, 0)) for q in range(-1, max(by)+1))
    cw = -0.5*sum((-1)**q*q*(lnd.get(q, 0) + lnd.get(q-1, 0)) for q in range(-1, max(by)+1))
    return c, betti, lnd, sd, cw
def tests():
    C3 = [frozenset(e) for e in ((0, 1), (1, 2), (0, 2))]
    simp = [frozenset(range(5))]
    RP2 = [frozenset(t) for t in ((0,1,2),(0,2,3),(0,3,4),(0,4,5),(0,5,1),(1,2,4),(2,3,5),(3,4,1),(4,5,2),(5,1,3))]
    for nm, F in (('C3', C3), ('simplexe Δ4', simp), ('RP² (6 sommets)', RP2)):
        c, b, _, sd, cw = torsion(F)
        print(f"  {nm:16s} : c = {c:+.10f} (e^c = {np.exp(c):.6f}) ; Betti réduits (ℚ) {b} ; superdét. {sd:+.1e} ; formule pondérée {cw:+.10f}")
    print(f"  référence : ½ ln 3 = {0.5*np.log(3):.10f} ; ln 2 = {np.log(2):.10f}")
print("== contrôles"); tests()
P_LIST = (2, 3, 5, 7)
def rank_modp(A, p):
    A = np.array(A, dtype=np.int64) % p; r = 0; rows, cols = A.shape
    for col in range(cols):
        piv = np.nonzero(A[r:, col])[0]
        if len(piv) == 0: continue
        i = r + piv[0]; A[[r, i]] = A[[i, r]]
        A[r] = (A[r]*pow(int(A[r, col]), p-2, p)) % p
        nz = np.nonzero(A[:, col])[0]; nz = nz[nz != r]
        A[nz] = (A[nz] - np.outer(A[nz, col], A[r])) % p
        r += 1
        if r == rows: break
    return r
# V32 (facettes jusqu'à 10 sommets) omis : trop coûteux en SVD dense sur la machine partagée
for nm, f in (('V37', 'V37_FHT_K26volmin.json'), ('V36', 'V36_FHT_K26vol.json')):
    F = load(f); c, b, lnd, sd, cw = torsion(F)
    print(f"== {nm} : c = {c:+.8f} ; e^c = {np.exp(c):.7f} ; e^2c = {np.exp(2*c):.7f} ; e^-c = {np.exp(-c):.4f}")
    print(f"     Betti réduits (ℚ, via SVD) {b} ; superdét. Σ(−1)^q ln det′Δ_q = {sd:+.1e} ; −½Σ(−1)^q q ln det′Δ_q = {cw:+.8f}")
    print("     ln det′(d_q* d_q) : " + " ".join(f"q={q}:{v:.3f}" for q, v in sorted(lnd.items())))
    if nm == 'V37':
        M, by = cobords(F)
        for p in P_LIST:
            rk = {q: rank_modp(A, p) for q, A in M.items()}
            n = {d: len(L) for d, L in by.items()}
            bp = [n[q] - rk.get(q, 0) - rk.get(q-1, 0) for q in range(0, max(by)+1)]
            print(f"     Betti réduits mod {p} : {bp}")
