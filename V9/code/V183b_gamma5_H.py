# V183b — canonicité : existe-t-il un Cl(5) covariant avec Γ5 = H/√3 (le cube tordu lui-même) ? Dimension de l'ensemble des solutions.
import numpy as np
exec(open("V183c_controle_commutant.py").read().split("rng = np.random.default_rng(1)")[0].replace('print("A1 herm', '(lambda *a: None)("A1 herm').replace('print("A1 commutatif', '(lambda *a: None)("').replace('print("dim commutant de A1h[0]', '(lambda *a: None)("'))
G5 = H/np.sqrt(3)
def mats2(p):
    G4 = sum(c*X for c, X in zip(p[:4], A1h)); Gv = sum(c*X for c, X in zip(p[4:], T1h)); return [Gv[0], Gv[1], Gv[2], G4, G5]
def cost2(p):
    G = mats2(p); return sum(np.sum(np.abs(G[i]@G[j]+G[j]@G[i]-(2*I8 if i==j else 0))**2) for i in range(5) for j in range(i,5))
rng = np.random.default_rng(7); sols = []
for t in range(40):
    r = minimize(cost2, rng.normal(size=12), method="BFGS", options=dict(gtol=1e-12, maxiter=4000))
    if r.fun < 1e-13: sols.append(r.x)
print("Γ5 = H/√3 : solutions %d / 40" % len(sols))
if sols:
    G = mats2(sols[0]); G6 = G[0]@G[1]@G[2]@G[3]@G[4]
    print(" Γ6 =", np.round(G6[0,0], 6), "scalaire :", np.allclose(G6, G6[0,0]*I8, atol=1e-6))
    # dimension de la variété des solutions : rang du jacobien des contraintes au point
    from scipy.optimize import approx_fprime
    def res(p):
        G = mats2(p); out = []
        for i in range(5):
            for j in range(i, 5):
                A = G[i]@G[j]+G[j]@G[i]-(2*I8 if i==j else 0); out += list(A.real.ravel()) + list(A.imag.ravel())
        return np.array(out)
    J = np.array([(res(sols[0] + 1e-6*e) - res(sols[0] - 1e-6*e))/2e-6 for e in np.eye(12)]).T
    s = np.linalg.svd(J, compute_uv=False); print(" dimension locale de l'ensemble des solutions (Γ5 fixé) :", 12 - np.sum(s > 1e-5*s[0]))
    # les γ_i trouvés sont-ils liés aux sauts du cube ? recouvrement avec γ_i(sauts)
    gam = []
    for ax in range(3):
        g = np.zeros_like(H)
        for a, b in ced:
            if abs(cube[a][ax]-cube[b][ax]) > 1: g[a, b] = H[a, b]; g[b, a] = H[b, a]
        gam.append(g)
    print(" |<γ_i(sol), saut_i>|/8 :", [round(abs(np.trace(G[i] @ gam[i]))/8, 4) for i in range(3)])
