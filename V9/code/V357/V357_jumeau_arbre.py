"""V357 (relecture de V341) : l'arbre porte-t-il un Weyl pour ε < 1 − 1/√3 (jumeau) ? Constructeur numpy de V352 (sans flux), t_f = ε, μ = r = t_c = 1.
Pour k0 petit (racine), états de basse énergie : |E|, poids sur la couche fine, ⟨F⟩ (goût), ⟨chiralité⟩ estimée par le signe de la vitesse."""
import sys, numpy as np
import os; sys.path[:0] = [os.path.join(os.path.dirname(os.path.abspath(__file__)), d) for d in ('v352run', '../V352')]; import V352_arbre as M
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2
dep = M.depth_of(N); Nc = len(dep)
Fop = np.kron(np.eye(Nc), 1j * M.Am @ M.Bm)
for eps in [float(x) for x in (sys.argv[2:] or ['1.0', '0.8', '0.6', '0.2', '0.1', '0.0'])]:
    out = []
    for kx in (0.0, 0.4, 0.8):
        k0 = (kx, 0., 0.)
        H = M.Hnat(N, k0, +1, mu=1.0, tf=eps); E, V = np.linalg.eigh(H)
        P = (np.abs(V) ** 2).reshape(Nc, 8, -1).sum(1)
        fine = P[dep == N].sum(0); root = P[dep == 0].sum(0)
        i = np.argsort(np.abs(E))[:6]
        desc = [f'{E[j]:+.3f}(fin {fine[j]:.2f}, rac {root[j]:.2f}, F {np.real(np.vdot(V[:, j], Fop @ V[:, j])):+.2f})' for j in i[:4]]
        out.append(f'k0x={kx}: ' + ' '.join(desc))
    print(f'ε = {eps} :\n   ' + '\n   '.join(out), flush=True)
