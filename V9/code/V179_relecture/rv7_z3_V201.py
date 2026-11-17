# Relecture V200–V201 : cristal alterné de monopôles ℤ3 (construction indépendante, jauge axiale), modes nuls des « quarks ».
import numpy as np, itertools as it
from rv_common import ROTS, HX, HY, lift_diag, commutator_class, level_structure

def ks(d, r): return [1, (-1) ** r[0], (-1) ** (r[0] + r[1])][d]
def q_cube(c): return 1 if sum(c) % 2 == 0 else -1          # charge ℤ3 sortante par face (unités 2π/3)

def flux_plus(s, k, N):
    """flux ℤ3 le long de +e_k à travers la plaquette de coin s normale à e_k."""
    if s[k] >= 1:
        c = list(s); c[k] -= 1; return q_cube(c)                # sortant du cube situé du côté −e_k
    return -q_cube(s)                                           # entrant dans le cube situé du côté +e_k
def F(s, d1, d2, N):
    """circulation orientée e_d1 puis e_d2 (d1<d2) ; normale e_d1×e_d2 = ±e_k"""
    k = 3 - d1 - d2; orient = {(0, 1): 1, (1, 2): 1, (0, 2): -1}[(d1, d2)]
    return orient * flux_plus(s, k, N)

def z3_axial(N):
    """a[(s,d)] ∈ ℤ3 en jauge axiale : a_x = 0 ; a_y, a_z intégrés le long de x ; condition initiale sur le plan x=0."""
    a = {}
    for s in it.product(range(N), repeat=3):
        for d in range(3):
            if s[d] + 1 < N: a[(s, d)] = 0
    # plan x = 0 : a_y(0,y,z) = 0 ; a_z(0,y+1,z) = a_z(0,y,z) + F_yz(0,y,z) ; a_z(0,0,z) = 0
    for z in range(N - 1):
        for y in range(N - 1):
            a[((0, y + 1, z), 2)] = (a[((0, y, z), 2)] + F((0, y, z), 1, 2, N)) % 3
    # propagation en x : a_y(x+1) = a_y(x) + F_xy(x) ; a_z(x+1) = a_z(x) + F_xz(x)
    for x in range(N - 1):
        for y in range(N):
            for z in range(N):
                if y + 1 < N: a[((x + 1, y, z), 1)] = (a[((x, y, z), 1)] + F((x, y, z), 0, 1, N)) % 3
                if z + 1 < N: a[((x + 1, y, z), 2)] = (a[((x, y, z), 2)] + F((x, y, z), 0, 2, N)) % 3
    # contrôle : toutes les plaquettes
    bad = 0
    for s in it.product(range(N), repeat=3):
        for d1, d2 in [(0, 1), (0, 2), (1, 2)]:
            if s[d1] + 1 < N and s[d2] + 1 < N:
                e1 = list(s); e1[d1] += 1; e2 = list(s); e2[d2] += 1
                circ = a[(s, d1)] + a[(tuple(e1), d2)] - a[(tuple(e2), d1)] - a[(s, d2)]
                bad += (circ - F(s, d1, d2, N)) % 3 != 0
    assert bad == 0, bad
    return a

def H_cluster(N, tfun, phi):
    a = z3_axial(N); sites = list(it.product(range(N), repeat=3)); idx = {s: i for i, s in enumerate(sites)}
    H = np.zeros((N ** 3,) * 2, complex)
    for (s, d), av in a.items():
        u = list(s); u[d] += 1
        v = tfun(s, d) * ks(d, s) * np.exp(2j * np.pi * phi * av / 3)
        H[idx[tuple(u)], idx[s]] = v; H[idx[s], idx[tuple(u)]] = np.conj(v)
    return H, sites

print("=== octet isolé (N=2) ===")
for phi in [0, 1, -1]:
    H, sites = H_cluster(2, lambda s, d: 1.0, phi); pos = [2. * np.array(s) - 1 for s in sites]
    Us = [lift_diag(H, pos, R)[0] for R in ROTS]
    lam, sc = commutator_class(lift_diag(H, pos, HX)[0], lift_diag(H, pos, HY)[0])
    print(" φ=%+d : relevés 24/24 %s ; [Ux,Uy] = %+.3f ; niveaux %s ; spectre de H² %s"
          % (phi, all(U is not None for U in Us), lam.real, level_structure(H, Us), np.round(np.linalg.eigvalsh(H @ H), 4).tolist()))

print("=== amas dimérisés : modes |E| < 1e-8 (φ = 1, quarks) ===")
for eps in [0.3, 0.5, 1.0]:
    out = []
    for n in [1, 2, 3, 4]:
        H, sites = H_cluster(2 * n, lambda s, d: eps if s[d] % 2 == 0 else 1.0, 1)
        w, V = np.linalg.eigh(H); z = np.abs(w) < 1e-8; out.append(int(z.sum()))
        if n == 3 and eps == 0.5:
            N = 2 * n; corner = np.array([all(c in (0, 1, N - 2, N - 1) for c in s) for s in sites])
            P = V[:, z]; wz = np.real(np.sum(np.abs(P) ** 2 * corner[:, None])) / P.shape[1]
            frac_corner_sites = corner.mean()
    print(" ε=%.1f : nombre de modes nuls pour n = 1..4 : %s" % (eps, out))
print(" n=3, ε=0.5 : poids moyen des modes nuls sur les octets de coin = %.3f (fraction de sites dans ces octets = %.3f)" % (wz, frac_corner_sites))
# leptons pour comparaison
H, _ = H_cluster(6, lambda s, d: 0.5 if s[d] % 2 == 0 else 1.0, 0)
print(" leptons (φ=0), n=3, ε=0.5 : 10 plus petits |E| =", np.round(np.sort(np.abs(np.linalg.eigvalsh(H)))[:10], 4).tolist())
# robustesse : rotation sur le spectre (invariance de jauge) — spectre identique avec une jauge ℤ3 aléatoire
H1, _ = H_cluster(4, lambda s, d: 0.5 if s[d] % 2 == 0 else 1.0, 1)
D = np.diag(np.exp(2j * np.pi * np.random.default_rng(1).integers(0, 3, 64) / 3))
print(" jauge ℤ3 aléatoire : spectre inchangé :", np.allclose(np.linalg.eigvalsh(D @ H1 @ D.conj().T), np.linalg.eigvalsh(H1)))

# lemme : cube à sauts unitaires, H² ∝ 1 ⇔ flux π sur chaque face
rng = np.random.default_rng(3)
X = list(it.product([0, 1], repeat=3)); edges = [(i, j) for i, j in it.combinations(range(8), 2) if sum(abs(np.array(X[i]) - np.array(X[j]))) == 1]
cnt = 0; okc = 0
for _ in range(2000):
    th = rng.choice([0, np.pi, 2 * np.pi / 3, -2 * np.pi / 3, np.pi / 3, -np.pi / 3], len(edges))
    Hc = np.zeros((8, 8), complex)
    for (i, j), t in zip(edges, th): Hc[i, j] = np.exp(1j * t); Hc[j, i] = np.exp(-1j * t)
    prop = np.allclose(Hc @ Hc, 3 * np.eye(8))
    allpi = True
    for ax in range(3):
        for v in range(2):
            f = [i for i in range(8) if X[i][ax] == v]; a0 = f[0]
            nb = [i for i in f if sum(abs(np.array(X[i]) - np.array(X[a0]))) == 1]; op = [i for i in f if sum(abs(np.array(X[i]) - np.array(X[a0]))) == 2][0]
            ph = Hc[a0, nb[0]] * Hc[nb[0], op] * Hc[op, nb[1]] * Hc[nb[1], a0]
            allpi &= np.isclose(ph, -1)
    cnt += prop; okc += (prop == allpi)
print(" lemme H²=3·1 ⇔ flux π sur les 6 faces : vérifié sur %d/2000 cubes aléatoires (%d avec H² ∝ 1)" % (okc, cnt))
