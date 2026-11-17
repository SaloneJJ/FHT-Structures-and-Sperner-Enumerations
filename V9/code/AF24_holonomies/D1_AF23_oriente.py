"""D1 -- AF23 avec atlas oriente (Gamma_4 pseudo-scalaire). Criteres : 02_theoreme_holonomies/CRITERES_D1_AF23_oriente.md (ecrits avant calcul).
Usage : python3 D1_AF23_oriente.py 2 3 | tee D1.log   (N = 2 : ~2 min ; N = 3 : > 10 min)"""
import sys, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
import os
HERE = os.path.dirname(os.path.abspath(__file__))
for p in (os.path.join(HERE, '..', 'AF23_dirac_S3'), os.path.join(HERE, '..', 'V352')): sys.path.insert(0, p)
import AF23_dirac_natif_S3 as A
o = {(0,0):1,(1,1):1,(2,0):1,(3,1):1,(0,1):-1,(1,0):-1,(2,1):-1,(3,0):-1}
def ops(N, mu):
    A.MU = mu; G4 = A.G4
    H, cells = A.build(N, 'S3')
    A.G4 = 0 * G4; H0, _ = A.build(N, 'S3'); A.G4 = G4
    HG = (H - H0).tocsr()
    sgn = np.repeat([o[c[1]] for c in cells], 8)
    Hor = (H0 + sp.diags(sgn.astype(float)) @ HG).tocsr()
    assert abs(Hor - Hor.conj().T).max() < 1e-12
    # controle : O = Id redonne H
    assert abs((H0 + HG) - H).max() < 1e-12
    # HG interne aux facettes : sgn HG = HG sgn
    assert abs(sp.diags(sgn.astype(float)) @ HG - HG @ sp.diags(sgn.astype(float))).max() < 1e-12
    return H, Hor
def small(H, k=24):
    n = H.shape[0]
    if n <= 6000: E = np.linalg.eigvalsh(H.toarray())
    else: E = sla.eigsh(H.tocsc(), k=k, sigma=1e-6, which='LM', return_eigenvectors=False)
    return np.sort(np.abs(E))[:k]
for N in map(int, sys.argv[1:]):
    t = time.time(); res = {}
    for mu in (1.0, 0.0):
        H, Hor = ops(N, mu)
        res[('brut', mu)] = small(H); res[('oriente', mu)] = small(Hor)
    for key, v in res.items(): print(f'N = {N} {key[0]:8s} mu = {key[1]} : |E| min = {v[0]:.4f} ; 12 plus petites :', np.round(v[:12], 4).tolist(), flush=True)
    m0b, m0o = res[('brut', 0.0)][0], res[('oriente', 0.0)][0]
    print(f'   rapports : brut {res[("brut",1.0)][0]/m0b:.3f} ; oriente {res[("oriente",1.0)][0]/m0o:.3f}  ({time.time()-t:.0f} s)', flush=True)
