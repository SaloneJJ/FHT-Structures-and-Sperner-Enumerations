"""V357 (relecture de V342 M3 / V343 R2) : gap des survivants avec l'insertion, sur l'arbre N = 2 (constructeur numpy de V352, sans flux, phase de Bloch k0).
Pour quelques k0, plus petite |E| des états de poids > 0,5 sur la couche fine ; test E(ε_D)² − E(0)² = ε_D² (O = A·Y123) et E(ε_D) = E(0) (O = Γ5)."""
import sys, numpy as np
import os; sys.path[:0] = [os.path.join(os.path.dirname(os.path.abspath(__file__)), d) for d in ('v352run', '../V352')]; import V352_arbre as M
N = 2; dep = M.depth_of(N)
def emin(k0, eD, O):
    H = M.Hins(N, k0, eD, O) if eD > 0 else M.Hins(N, k0, 1e-14, 'AY')
    E, V = np.linalg.eigh(H); P = (np.abs(V) ** 2).reshape(2, -1, 8, V.shape[1]).sum(axis=(0, 2))
    sel = P[dep == N].sum(0) > 0.5
    return np.sort(np.abs(E[sel]))[:2]
for k0 in [(0.4, 0, 0), (0.8, 0.3, 0), (2 * np.pi / 8, 2 * np.pi / 8, 2 * np.pi / 8)]:
    e0 = emin(k0, 0, 'AY')[0]
    out = []
    for eD in (0.05, 0.1, 0.3):
        ea = emin(k0, eD, 'AY')[0]; eg = emin(k0, eD, 'G5')[0]
        out.append(f'ε_D={eD}: √(E²−E0²)={np.sqrt(max(ea**2 - e0**2, 0)):.4f} ; Γ5 : E−E0={eg - e0:+.1e}')
    print(f'k0 = {np.round(k0, 3)} : E0 = {e0:.4f} | ' + ' | '.join(out))
