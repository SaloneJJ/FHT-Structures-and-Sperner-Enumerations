"""Doublets exactement degeneres (q = 2) : l estimateur sum_i w_i s_i depend-il de la base dans le doublet ?
On teste si D = dH/dtheta et P_n restreints au doublet sont proportionnels a l identite. Usage : python3 rv3_degen.py N q th"""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
N, q, th = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]); A = list(sys.argv); g = C.setup(N, q); Nc, depth = g['Nc'], g['depth']
if N >= 4:
    z = np.load(f'/home/claude/fht/relecture3/fourier_N4_q{int(q)}.npz'); fast = C.Fast(z['ks'], z['Cf'], z['keys'], int(z['n4']))
else:
    ks, Cf, keys, n4, *_ = C.fourier_fit(kmax_factor=1.0); fast = C.Fast(ks, Cf, keys, n4)
H = fast.H(th); D = fast.H(th, 1); Hc = H.tocsc()
lu = sla.splu(Hc - 1e-5 * sp.identity(Hc.shape[0], format='csc'), permc_spec='MMD_AT_PLUS_A')
E, V = sla.eigsh(H, k=12, sigma=1e-5, which='LM', tol=1e-13, OPinv=sla.LinearOperator(H.shape, matvec=lu.solve, dtype=complex))
o = np.argsort(E); E, V = E[o], V[:, o]
print(f'N = {N}, q = {q}, theta = {th} : E = {np.round(E, 6)}')
i = 0
while i < len(E) - 1:
    if abs(E[i + 1] - E[i]) < 1e-3:
        W = V[:, i:i + 2]; Dm = W.conj().T @ (D @ W); devD = np.abs(Dm - np.trace(Dm) / 2 * np.eye(2)).max()
        devP = 0
        for n in range(N + 1):
            msk = np.repeat(depth == n, 4); Pm = W[msk].conj().T @ W[msk]; devP = max(devP, np.abs(Pm - np.trace(Pm) / 2 * np.eye(2)).max())
        print(f'  doublet E = {E[i]:+.6f} (ecart {E[i+1]-E[i]:.1e}) : |D - tr/2| = {devD:.1e} (pente {np.real(np.trace(Dm))/2:+.4f}) ; max_n |P_n - tr/2| = {devP:.1e}')
        i += 2
    else: i += 1
