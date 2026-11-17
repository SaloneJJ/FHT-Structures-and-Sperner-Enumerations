"""Balayage fin en theta de l integrande du flot fenetre (M1, fenetres eps et eps/2) + etats bas (E, pente, profil).
H et dH/dtheta : reconstruction de Fourier independante (rv3_core), verifiee contre V162. Diagonalisation creuse (shift-invert),
K etats les plus proches de 0. Sauvegarde npz par theta.
Usage : python3 rv3_scan.py N q th_min th_max npts K tag [fourier_npz]"""
import sys, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
N, q, a, b, npts, K, tag = int(sys.argv[1]), float(sys.argv[2]), float(eval(sys.argv[3], {'pi': np.pi})), float(eval(sys.argv[4], {'pi': np.pi})), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
ARGV = list(sys.argv); t0 = time.time(); g = C.setup(N, q); L = g["L"]
if len(ARGV) > 8:
    z = np.load(ARGV[8]); fast = C.Fast(z['ks'], z['Cf'], z['keys'], int(z['n4']))
else:
    ks, Cf, keys, n4, *_ = C.fourier_fit(kmax_factor=1.0); fast = C.Fast(ks, Cf, keys, n4)
EPS = 0.25 * 2 * np.pi / L
ths = np.linspace(a, b, npts)
out = dict(ths=ths, E=np.full((npts, K), np.nan), S=np.full((npts, K), np.nan), W=np.full((npts, N + 1, K), np.nan), eps=EPS, N=N, q=q, L=L)
for i, th in enumerate(ths):
    H = fast.H(th); D = fast.H(th, 1)
    Hc = H.tocsc(); lu = sla.splu(Hc - 1e-5 * sp.identity(Hc.shape[0], format='csc'), permc_spec='MMD_AT_PLUS_A')
    E, V = sla.eigsh(H, k=K, sigma=1e-5, which='LM', tol=1e-12, OPinv=sla.LinearOperator(H.shape, matvec=lu.solve, dtype=complex)); del lu
    o = np.argsort(E); E, V = E[o], V[:, o]
    out['E'][i] = E; out['S'][i] = np.real(np.einsum('ij,ij->j', V.conj(), D @ V)); out['W'][i] = C.profs(V)
    if np.sum(np.abs(E) < EPS) >= K - 1: print(f'  ATTENTION fenetre saturee a theta = {th}', flush=True)
    if i % 20 == 0 or N >= 4: print(f'  {i}/{npts} theta = {th:.5f}  t = {time.time() - t0:.0f} s', flush=True); np.savez(f'/home/claude/fht/relecture3/scan_{tag}.npz', **out)
np.savez(f'/home/claude/fht/relecture3/scan_{tag}.npz', **out)
print(f'  fini, {time.time() - t0:.0f} s')
