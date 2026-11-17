"""V359 : contrôle indépendant des parités de V349 K1 / V333 : (Γ4, FΓ5, S) sous les relevés des 48 éléments de O_h
sur l'octet BBH (forme séparable de l'article VII), F = −Γ4 S ; puis argument de symétrie pour le flot (ω, μ).
Code propre (pas d'import de rvA_common). Usage : python3 V359_parites.py   (numpy)."""
import numpy as np, itertools as it
tf, tc = 0.4, 1.0
sz = np.diag([1., -1]).astype(complex); I2 = np.eye(2)
def h(k): return np.array([[0, tf + tc*np.exp(-1j*k)], [tf + tc*np.exp(1j*k), 0]])
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
def H(k): return k3(h(k[0]), I2, I2) + k3(sz, h(k[1]), I2) + k3(sz, sz, h(k[2]))
SITES = list(it.product((0, 1), repeat=3)); POS = np.array(SITES) - 0.5
S = np.diag([(-1)**sum(s) for s in SITES]).astype(complex)
K = np.full(3, np.pi); n2 = lambda M: M/np.sqrt(np.real(np.trace(M@M))/8)
G5 = n2(H(K)); d = 1e-4
al = [n2((H(K + d*e) - H(K - d*e))/(2*d)) for e in np.eye(3)]
cs = [n2((H(K + d*e) + H(K - d*e) - 2*H(K))/d**2) for e in np.eye(3)]
A = n2((cs[0] - cs[1])/np.sqrt(2)); B = n2((cs[0] + cs[1] - 2*cs[2])/np.sqrt(6))
G4 = G5@al[0]@al[1]@al[2]; F = 1j*A@B; FG5 = F@G5
for sg in (+1, -1):
    print(f"F − ({sg:+d}) Γ4 S : {np.abs(F - sg*G4@S).max():.1e}")
print(f"[F,S] {np.abs(F@S-S@F).max():.1e} ; [Γ4,S] {np.abs(G4@S-S@G4).max():.1e} ; F* + F (base des sites) {np.abs(F.conj()+F).max():.1e}")
def lift(g):
    perm = [int(np.argmin(np.linalg.norm(POS - g@p, axis=1))) for p in POS]
    P = np.zeros((8, 8)); P[perm, range(8)] = 1
    # phases diagonales : U = D P avec U H(k) U† = H(gk) ; on résout D par moindres carrés sur des k aléatoires
    rng = np.random.default_rng(1); M = []
    for _ in range(5):
        k = rng.uniform(-np.pi, np.pi, 3); Hk = P@H(k)@P.T; Hg = H(g@k)
        # D Hk D† = Hg  ⇔  d_i Hk_ij conj(d_j) = Hg_ij  ⇔  d_i Hk_ij − Hg_ij d_j = 0 (|d| = 1)
        for i in range(8):
            for j in range(8):
                if abs(Hk[i, j]) > 1e-12 or abs(Hg[i, j]) > 1e-12:
                    row = np.zeros(8, complex); row[i] += Hk[i, j]; row[j] -= Hg[i, j]; M.append(row)
    u, s, vh = np.linalg.svd(np.array(M)); dvec = vh[-1].conj(); dvec /= dvec[np.argmax(abs(dvec))]
    U = np.diag(dvec)@P
    k = np.array([0.3, -1.1, 2.0]); ok = np.abs(U@H(k)@U.conj().T - H(g@k)).max() < 1e-8 and s[-2] > 1e-6
    return U, ok
def sgn(U, O):
    R = U@O@U.conj().T
    for s in (+1, -1):
        if np.abs(R - s*O).max() < 1e-6: return s
    return 0
tab = {}
for p in it.permutations(range(3)):
    for e in it.product((1, -1), repeat=3):
        g = np.diag(e)@np.eye(3)[list(p)]; U, ok = lift(g)
        key = (sgn(U, G4), sgn(U, FG5), sgn(U, S)); det = round(np.linalg.det(g))
        tab.setdefault(key, []).append((det, ok))
for key, L in sorted(tab.items(), reverse=True):
    print(f"(Γ4, FΓ5, S) = {key} : {len(L)} éléments ; det = {sorted(set(x[0] for x in L))} ; relevé vérifié : {all(x[1] for x in L)}")
print("Lecture : un élément qui renverse Γ4 en gardant S existe (−,−,+) : une donnée de type S (marque L/R) ne fixe pas ω ;\n"
      "un élément (+,−,−) renverse μ sans toucher Γ4 : en l'absence de terme en S, flot(ω, μ) = flot(ω, −μ) exactement.")
