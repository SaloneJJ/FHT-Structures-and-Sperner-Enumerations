# V193 — Q1 (saut R–R effectif, cube tordu ou non) et Q2 (réseau cubique à flux π dimérisé). Critères : V192.
import numpy as np, itertools as it
exec(open("V180_FHT_spin_cube_monopole.py").read().split("# ---------- octet (cube) ----------")[0])
X = [np.array(c) for c in it.product([0,1], repeat=3)]          # coins du cube (bits)
idx = {tuple(x): i for i, x in enumerate(X)}
def eta(d, x): return [1, (-1)**x[0], (-1)**(x[0]+x[1])][d]    # phases de Kawamoto–Smit
def cube(twist):
    H = np.zeros((8, 8))
    for x in X:
        for d in range(3):
            if x[d] == 0:
                y = x.copy(); y[d] = 1; s = eta(d, x) if twist else 1
                H[idx[tuple(x)], idx[tuple(y)]] = H[idx[tuple(y)], idx[tuple(x)]] = s
    return H
R = [i for i, x in enumerate(X) if x.sum() % 2 == 0]; Lc = [i for i in range(8) if i not in R]
print("=== Q1 : saut R–R effectif au second ordre (A² restreint à R) ===")
for lab, tw in [("non tordu", False), ("tordu (flux π)", True)]:
    A = cube(tw); M = (A @ A)[np.ix_(R, R)]
    print(" %-15s A²|_R =" % lab, M.astype(int).tolist(), " valeurs propres :", np.round(np.linalg.eigvalsh(M), 6).tolist())
# tétraèdre R tordu π/2 (twist de A_F) : V180
cubep = [2.*x - 1 for x in X]; R4 = [cubep[i] for i in R]
ted = list(it.combinations(range(4), 2)); tf = [outward_cycle(R4, list(f)) for f in it.combinations(range(4), 3)]
HT = build(R4, ted, tf, flux_list(4, np.pi/2, 1)); print(" tétraèdre R, π/2 par face : valeurs propres", np.round(np.linalg.eigvalsh(HT), 4).tolist())
HT0 = build(R4, ted, tf, flux_list(4, 0, 0)); print(" tétraèdre R, sans flux     : valeurs propres", np.round(np.linalg.eigvalsh(HT0), 4).tolist())
print("=== Q2 : réseau cubique à flux π dimérisé (maille = octet) ===")
def Hk(k, tf_, tc_):
    H = np.zeros((8, 8), complex)
    for x in X:
        for d in range(3):
            y = x.copy(); y[d] = 1 - x[d]
            if x[d] == 0:   # lien frère (intra-octet), de x vers x + e_d
                H[idx[tuple(y)], idx[tuple(x)]] += tf_*eta(d, x)
            else:           # lien cousin : de x (x_d = 1) vers la maille suivante, coin x − e_d
                H[idx[tuple(y)], idx[tuple(x)]] += tc_*eta(d, x)*np.exp(1j*k[d])
    return H + H.conj().T
K = np.array([np.pi]*3)
Hc = cube(True)
for tf_, tc_ in [(1.0, 1.0), (1.0, 0.8), (0.8, 1.0), (1.0, 0.5)]:
    H0 = Hk(K, tf_, tc_)
    ok_a = np.allclose(H0, (tf_ - tc_)*Hc)
    w0 = np.linalg.eigvalsh(H0)
    # dispersion le long de plusieurs directions
    out = []
    for dirv in [np.array([1,0,0]), np.array([1,1,0])/np.sqrt(2), np.array([1,1,1])/np.sqrt(3), np.array([0.3,-0.8,0.52])/np.linalg.norm([0.3,-0.8,0.52])]:
        q = 1e-3; w = np.linalg.eigvalsh(Hk(K + q*dirv, tf_, tc_)); Em = np.min(np.abs(w))
        v2 = (Em**2 - 3*(tf_-tc_)**2)/q**2
        out.append(round(np.sqrt(max(v2, 0)), 4))
    degs = np.round(np.linalg.eigvalsh(Hk(K + 1e-2*np.array([0.3,-0.8,0.52]), tf_, tc_)), 6)
    print(" t_f=%.2f t_c=%.2f : H(K) = (t_f − t_c) H_cube ? %s ; E(K) = %s ; v le long de 4 directions = %s ; dégénérescences près de K : %s"
          % (tf_, tc_, ok_a, sorted(set(np.round(w0, 6))), out, [int(np.sum(np.isclose(degs, e))) for e in sorted(set(degs))]))
# structure du terme linéaire
tf_, tc_ = 1.0, 0.8; h = 1e-6
Xd = [ (Hk(K + h*e, tf_, tc_) - Hk(K - h*e, tf_, tc_))/(2*h) for e in np.eye(3) ]
Xd = [ (x + x.conj().T)/2 for x in Xd ]
print(" ∂_i H : {∂_i, ∂_j} ∝ δ_ij ?", [np.allclose(Xd[i]@Xd[j] + Xd[j]@Xd[i], 0, atol=1e-6) for i, j in [(0,1),(0,2),(1,2)]],
      "; (∂_i H)² =", [round(np.real(np.trace(x@x))/8, 4) for x in Xd],
      "; {∂_i H, H_cube} = 0 ?", [np.allclose(x@Hc + Hc@x, 0, atol=1e-6) for x in Xd])
# covariance : ∂_i H dans l'espace T1 des relevés du cube tordu (jauge de Kawamoto–Smit)
Us = [lift(Hc.astype(complex), cubep, Rm)[0] for Rm in ROT]
rows = [np.kron(np.eye(3), np.kron(U, U.conj())) - np.kron(Rm.T, np.eye(64)) for U, Rm in zip(Us, ROT)]
A = np.vstack(rows); u, s, vt = np.linalg.svd(A); T1 = vt[np.sum(s > 1e-9):].conj()
v = np.concatenate([x.ravel() for x in Xd]); P = T1.T @ np.linalg.lstsq(T1.T, v, rcond=None)[0]
print(" dim T1 = %d ; ∂H dans T1 : résidu relatif = %.2e" % (len(T1), np.linalg.norm(v - P)/np.linalg.norm(v)))
# la masse est-elle singulet de goût ? commutant de l'algèbre engendrée par {∂_i H, H_cube}
G = Xd + [Hc/np.sqrt(3)]
rowsc = [np.kron(g, np.eye(8)) - np.kron(np.eye(8), g.T) for g in G]
uc, sc, vtc = np.linalg.svd(np.vstack(rowsc)); com = vtc[np.sum(sc > 1e-6):]
print(" commutant de {∂_i H, H_cube} : dimension %d (4 attendu : M₂ du goût)" % len(com))
