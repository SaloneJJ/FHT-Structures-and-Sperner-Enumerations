"""Point 3 : relation spectrale entre H(theta) et H(-theta) (N = 3, diagonalisation complete), puis recherche de l'operateur.
Usage : python3 rv3_sym.py N q"""
import sys, itertools, numpy as np, scipy.sparse as sp
sys.path.insert(0, '/home/claude/fht/relecture3'); import rv3_core as C
N, q = int(sys.argv[1]), float(sys.argv[2]); g = C.setup(N, q); Nc = g['Nc']
for th in (0.37, 1.9):
    Ep, Vp = np.linalg.eigh(C.Href(th).toarray()); Em, Vm = np.linalg.eigh(C.Href(-th).toarray())
    print(f'theta = {th} : |spec H(-th) - spec H(th)| {np.abs(np.sort(Em) - np.sort(Ep)).max():.1e} ; |spec H(-th) - spec(-H(th))| {np.abs(np.sort(Em) - np.sort(-Ep)).max():.1e}')
# recherche : U H(th)^* U^+ = s H(-th), U = P_spatial (reflexions des coordonnees, niveau par niveau) x produit de Gamma, a une jauge diagonale pres
GAM = g['GAM']; sides = g['sides']; lin = g['lin']; grid = g['grid']
def perm(refl):
    p = np.zeros(Nc, int)
    for n in range(N + 1):
        s = sides[n]; I, J, K = grid(s); src = lin(n, I, J, K).ravel()
        I2, J2, K2 = [(s - 1 - A) if f else A for A, f in zip((I, J, K), refl)]
        p[src] = lin(n, I2, J2, K2).ravel()
    return p
gprods = {}
for sub in itertools.product((0, 1), repeat=5):
    M = np.eye(4, dtype=complex)
    for j, b in enumerate(sub):
        if b: M = M @ GAM[j]
    gprods[sub] = M
th = 0.37; Hp = C.Href(th).toarray(); Hm = C.Href(-th).toarray()
found = []
for conj in (0, 1):
    H1 = Hp.conj() if conj else Hp
    for refl in itertools.product((0, 1), repeat=3):
        p = perm(refl); P = sp.csr_matrix((np.ones(Nc), (p, np.arange(Nc))), shape=(Nc, Nc))
        for sub, gm in gprods.items():
            U = sp.kron(P, sp.csr_matrix(gm)).toarray()
            X = U @ H1 @ U.conj().T
            for s in (1, -1):
                # egalite des modules (independante d une jauge diagonale)
                if np.abs(np.abs(X) - np.abs(s * Hm)).max() < 1e-10:
                    # jauge diagonale : X = Lam (s Hm) Lam^* ; on la resout sur un arbre couvrant par propagation
                    T = s * Hm; n4 = X.shape[0]; lam = np.full(n4, np.nan + 0j); lam[0] = 1; stack = [0]
                    while stack:
                        i = stack.pop()
                        for j in np.where(np.abs(T[i]) > 1e-12)[0]:
                            if np.isnan(lam[j]):
                                # X_ij = lam_i T_ij conj(lam_j)
                                lam[j] = np.conj(X[i, j] / (lam[i] * T[i, j])); stack.append(j)
                    ok = not np.isnan(lam).any()
                    err = np.abs(np.diag(lam) @ T @ np.diag(lam.conj()) - X).max() if ok else np.inf
                    found.append((conj, refl, sub, s, err))
                    print(f'  candidat : conj {conj}, reflexions {refl}, Gamma {sub}, signe {s:+d} : modules egaux ; avec jauge diagonale, ecart {err:.1e}')
print('aucun candidat' if not found else f'{len(found)} candidats')
