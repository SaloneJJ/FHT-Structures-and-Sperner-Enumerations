import sys, time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
g = C.setup(4, 2.0); z = np.load('fourier_N4_q2.npz'); fast = C.Fast(z['ks'], z['Cf'], z['keys'], int(z['n4']))
H = fast.H(0.3).tocsc(); n = H.shape[0]
for spec in ['MMD_AT_PLUS_A', 'COLAMD']:
    t = time.time(); lu = sla.splu(H - 1e-5 * sp.identity(n, format='csc'), permc_spec=spec); t1 = time.time() - t
    print(spec, f'factorisation {t1:.1f} s, nnz(L+U) {lu.L.nnz + lu.U.nnz:.3e}', flush=True)
    OP = sla.LinearOperator(H.shape, matvec=lu.solve, dtype=complex)
    t = time.time(); E, V = sla.eigsh(H, k=24, sigma=1e-5, which='LM', OPinv=OP); print(f'  eigsh {time.time()-t:.1f} s, E = {np.sort(np.abs(E))[:4]}', flush=True)
