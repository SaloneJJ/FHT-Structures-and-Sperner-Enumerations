"""V359 : contrôle indépendant de V353 (d) — dépendance en φ de l'énergie du vide et d'autres fonctionnelles spectrales.
Construction propre : Clifford Cl(6) en 8 dimensions par une autre affectation de Pauli, puis conjugaison unitaire aléatoire
(rien n'est importé de V353). Modèle : H_s(k) = Σ_d [−sin k_d Y_d + s r (1 − cos k_d) X_d] + μ F Γ5 (+ m0 s Γ5 en option),
deux copies s = ±1 (γ_F = σ_z), couplées par ε O_φ ⊗ σ_x, O_φ = (cos φ A + sin φ B) Y1Y2Y3 ; grille k NK³ non décalée.
Mesures : harmoniques c2, c4, c6 de E_vac(φ) = Σ_{E<0} E et de S_t(φ) = Σ_k Tr e^{−t H²} ; signe de c6 (minimum à 0 ou 30 mod 60°).
Usage : python3 V359_phi.py [NK]   (numpy)."""
import sys, itertools as it, numpy as np
NK = int(sys.argv[1]) if len(sys.argv) > 1 else 10
s0 = np.eye(2); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1]).astype(complex)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
# autre affectation (ordre inversé des facteurs, rôles de σx et σy échangés)
Xr = [k3(s0, s0, sy), k3(s0, sy, sz), k3(sy, sz, sz)]; Yr = [k3(s0, s0, sx), k3(s0, sx, sz), k3(sx, sz, sz)]
rng = np.random.default_rng(359); Q, _ = np.linalg.qr(rng.normal(size=(8, 8)) + 1j*rng.normal(size=(8, 8)))
X = [Q@x@Q.conj().T for x in Xr]; Y = [Q@y@Q.conj().T for y in Yr]
G = X + Y
err = max(np.abs(a@b + b@a - 2*(i == j)*np.eye(8)).max() for i, a in enumerate(G) for j, b in enumerate(G))
print(f"Clifford : max |{{Γ_i, Γ_j}} − 2δ_ij| = {err:.1e}")
G5 = sum(X)/np.sqrt(3); A = (X[0]-X[1])/np.sqrt(2); B = (X[0]+X[1]-2*X[2])/np.sqrt(6); Y123 = Y[0]@Y[1]@Y[2]
F = 1j*A@B; FG = F@G5
O = lambda p: (np.cos(p)*A + np.sin(p)*B)@Y123
U = lambda th: np.cos(th)*np.eye(8) + 1j*np.sin(th)*F
print(f"F² = 1 : {np.abs(F@F - np.eye(8)).max():.1e} ; [F, G5] {np.abs(F@G5-G5@F).max():.1e} ; [F, Y123] {np.abs(F@Y123-Y123@F).max():.1e} ; "
      f"U_θ O_φ U_θ† − O_(φ+2θ) : {np.abs(U(.3)@O(.2)@U(.3).conj().T - O(.8)).max():.1e}")
ks = 2*np.pi*np.arange(NK)/NK; KK = np.array(list(it.product(ks, ks, ks))); c = np.cos(KK); sn = np.sin(KK)
def Hs(s, mu, r, m0):
    H = np.zeros((len(KK), 8, 8), complex)
    for d in range(3): H += -sn[:, d, None, None]*Y[d] + s*r*(1 - c[:, d, None, None])*X[d]
    return H + mu*FG + s*m0*G5
PH = np.deg2rad(np.arange(0, 180, 5.0))
def run(mu, r, m0, eps, ts=(0.3, 1.0, 3.0)):
    Hp, Hm = Hs(+1, mu, r, m0), Hs(-1, mu, r, m0); ev, st = [], {t: [] for t in ts}
    for p in PH:
        H = np.zeros((len(KK), 16, 16), complex); H[:, :8, :8] = Hp; H[:, 8:, 8:] = Hm
        H[:, :8, 8:] = eps*O(p); H[:, 8:, :8] = eps*O(p).conj().T
        E = np.linalg.eigvalsh(H); ev.append(E[E < 0].sum())
        for t in ts: st[t].append(np.exp(-t*E**2).sum())
    def h(v):
        v = np.array(v) - np.mean(v); return {n: 2*np.mean(v*np.exp(-1j*n*PH)) for n in (2, 4, 6)}
    return h(ev), {t: h(v) for t, v in st.items()}, np.array(ev)
print(f"\nGrille {NK}³ non décalée. c_n = coefficient complexe de e^(inφ). « min » : position du minimum de la fonctionnelle (mod 60°).")
def loc(c6): return "0°" if np.real(c6) < 0 else "30°"
rows = [(1.0, 1, 0), (0.5, 1, 0), (1.5, 1, 0), (-1.0, 1, 0), (2.5, 1, 0), (1.0, -1, 0), (1.0, 1, -0.5), (1.0, 1, 0.5), (0.0, 1, 0)]
for (mu, r, m0) in rows:
    out = []
    for eps in (0.1, 0.2, 0.4):
        he, hs, _ = run(mu, r, m0, eps); out.append((eps, he, hs))
    a1, a2 = abs(out[0][1][6]), abs(out[1][1][6])
    p = np.log(a2/a1)/np.log(2) if a1 > 0 else np.nan
    print(f"μ={mu:+.1f} r={r:+d} m0={m0:+.1f} : |c2|,|c4| (ε=0,4) = {abs(out[2][1][2]):.0e}, {abs(out[2][1][4]):.0e} ; c6(E_vac) ε=0,1/0,2/0,4 = "
          + " / ".join(f"{np.real(o[1][6]):+.2e}" for o in out) + f" (Im ≤ {max(abs(np.imag(o[1][6])) for o in out):.0e}) ; p = {p:.2f} ; min E_vac à {loc(out[2][1][6])}")
    print("      fonctionnelles Tr e^(−tH²) (ε = 0,4) : " + " ; ".join(f"t={t}: c6 = {np.real(h[6]):+.2e} → min à {loc(h[6])}" for t, h in out[2][2].items()))
