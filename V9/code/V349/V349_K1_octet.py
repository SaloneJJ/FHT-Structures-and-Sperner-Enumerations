"""V349 K1 : parités (Γ4, FΓ5, S) sous σ_d, inversion, C4 avec le Clifford du code de C3' (V335/V340) ; identité F = ±Γ4 S.
Critères : 01_synthese/V349 §2 (écrits avant calcul). numpy seul.   JJ SALONE, FHT 2026"""
import numpy as np, itertools as it
from rvA_common import lift, classify, H as Hbbh, S as Ssub
sxm = np.array([[0, 1], [1, 0]], complex); sym = np.array([[0, -1j], [1j, 0]]); szm = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
Xg = [k3(sxm, I2, I2), k3(szm, sxm, I2), k3(szm, szm, sxm)]; Yg = [k3(sym, I2, I2), k3(szm, sym, I2), k3(szm, szm, sym)]
G5 = sum(Xg)/np.sqrt(3); Am = (Xg[0] - Xg[1])/np.sqrt(2); Bm = (Xg[0] + Xg[1] - 2*Xg[2])/np.sqrt(6)
G4 = G5@Yg[0]@Yg[1]@Yg[2]; F = 1j*Am@Bm; FG = F@G5
kk = np.array([0.3, -1.1, 2.0]); assert np.allclose(sum((0.4 + np.cos(kk[d]))*Xg[d] + np.sin(kk[d])*Yg[d] for d in range(3)), Hbbh(kk))
def sgn(U, X):
    Y = U@X@U.conj().T
    return '+' if np.allclose(Y, X) else ('-' if np.allclose(Y, -X) else '?')
for nm, sg in (('F = -G4 S', -1), ('F = +G4 S', +1)): print(nm, np.allclose(F, sg*G4@Ssub))
print('[F,S]=0', np.allclose(F@Ssub, Ssub@F), ' [G4,S]=0', np.allclose(G4@Ssub, Ssub@G4), ' F* = -F', np.allclose(F.conj(), -F))
els = {'inversion': -np.eye(3),
       'C4 (z)': np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.]]),
       'sigma_d (x<->y)': np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1.]]),
       'S4 (z)': np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1.]]),
       'C3 (xyz)': np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0.]]),
       'C2 (z)': np.diag([-1, -1, 1.])}
print(f"{'element':18s} det garde_T+  G5 G4 FG S F")
for nm, g in els.items():
    U, ok = lift(g)
    det, keeps = classify(g)
    print(f"{nm:18s} {det:+d}  {str(keeps):6s}  ", ' '.join(sgn(U, X) for X in (G5, G4, FG, Ssub, F)))
# tous les éléments de O_h : tableau des parités par classe
from rvA_common import ohgroup
cnt = {}
for g in ohgroup():
    U, ok = lift(g); det, keeps = classify(g)
    key = (det, keeps, sgn(U, G4), sgn(U, FG), sgn(U, Ssub)); cnt[key] = cnt.get(key, 0) + 1
print('O_h complet (det, garde les tétraèdres, G4, FG, S) : nombre')
for k, v in sorted(cnt.items(), key=str): print('  ', k, v)
