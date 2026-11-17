# Relecture V193 : réseau cubique à flux π dimérisé (maille = octet). Construction depuis l'espace réel (tore L=4 sites/côté).
import numpy as np, itertools as it
from rv_common import ROTS, lift_diag
rng = np.random.default_rng(193)

def ks(d, r): return [1, (-1) ** r[0], (-1) ** (r[0] + r[1])][d]
def bloch(k, tf, tc):
    """H(k) 8×8 obtenu par transformée de Fourier des liens d'un réseau infini (r = 2R + x)."""
    H = np.zeros((8, 8), complex); xs = list(it.product([0, 1], repeat=3))
    for a, x in enumerate(xs):
        for d in range(3):
            r = np.array(x); r2 = r.copy(); r2[d] += 1          # lien r -> r + e_d
            t = tf if r[d] % 2 == 0 else tc                     # dans l'octet si r_d pair
            R2 = r2 // 2; x2 = tuple(r2 % 2); b = xs.index(x2)
            val = t * ks(d, r) * np.exp(1j * np.dot(k, R2))     # phase de Bloch e^{ik·R}
            H[b, a] += val; H[a, b] += np.conj(val)
    return H
xs = list(it.product([0, 1], repeat=3)); pos = [2. * np.array(x) - 1 for x in xs]
Hcube = bloch(np.zeros(3), 1.0, 0.0)
K = np.array([np.pi] * 3)

# (a) H(K) et H(0)
for tf, tc in [(1.0, 0.7), (0.4, 1.0)]:
    print("t_f=%.1f t_c=%.1f : H(K) = (t_f−t_c)H_cube : %s ; H(0) = (t_f+t_c)H_cube : %s"
          % (tf, tc, np.allclose(bloch(K, tf, tc), (tf - tc) * Hcube), np.allclose(bloch(np.zeros(3), tf, tc), (tf + tc) * Hcube)))

# (b) spectre exact : E² = Σ_d (t_f² + t_c² + 2 t_f t_c cos k_d), 4 fois dégénéré partout
err = 0
for _ in range(200):
    tf, tc = rng.random(2) * 2; k = rng.uniform(-np.pi, np.pi, 3)
    w = np.linalg.eigvalsh(bloch(k, tf, tc)); E2 = np.sum(tf ** 2 + tc ** 2 + 2 * tf * tc * np.cos(k))
    err = max(err, np.max(np.abs(np.sort(w) - np.array([-1] * 4 + [1] * 4) * np.sqrt(E2))))
print("E² = Σ_d (t_f²+t_c²+2t_f t_c cos k_d), 4+4 dégénéré (200 points aléatoires) : écart max %.1e" % err)
print("  => en K+q : E² = 3(t_f−t_c)² + 2 t_f t_c Σ(1−cos q_d) ≈ 3(t_f−t_c)² + t_f t_c |q|² : v² = t_f t_c ; le gap ne se ferme qu'à t_f = t_c")

# (c) décomposition exacte pour q ARBITRAIRE (pas seulement petit)
G5 = Hcube / np.sqrt(3)
gam = [bloch(np.eye(3)[d] * 0, 1.0, 0.0) for d in range(3)]
# sauts du cube par direction
def hop_dir(d):
    g = np.zeros((8, 8), complex)
    for a, x in enumerate(xs):
        if x[d] == 0:
            y = list(x); y[d] = 1; b = xs.index(tuple(y)); g[a, b] = g[b, a] = ks(d, x)
    return g
gam = [hop_dir(d) for d in range(3)]
alpha = []
for d in range(3):
    e = np.eye(3)[d]; h = 1e-6
    alpha.append((bloch(K + h * e, 0, 1.0) - bloch(K - h * e, 0, 1.0)) / (2 * h))
gperp = [g - G5 / np.sqrt(3) for g in gam]
worst = 0
for _ in range(100):
    tf, tc = rng.random(2) * 2; q = rng.uniform(-np.pi, np.pi, 3)
    m = np.sqrt(3) * (tf - tc); rW = tc / np.sqrt(3)
    pred = (m + rW * np.sum(1 - np.cos(q))) * G5 + tc * sum(np.sin(q[d]) * alpha[d] for d in range(3)) + tc * sum((1 - np.cos(q[d])) * gperp[d] for d in range(3))
    worst = max(worst, np.abs(bloch(K + q, tf, tc) - pred).max())
print("décomposition H(K+q) = [m + r_W Σ(1−cos q)]Γ5 + t_c Σ sin q α + t_c Σ(1−cos q) γ⊥, q arbitraire : écart max %.1e" % worst)
six = gam + [a for a in alpha]
print("les six matrices {γ_d, α_d} anticommutent deux à deux (forme BBH) :",
      all(np.allclose(six[i] @ six[j] + six[j] @ six[i], 2 * np.eye(8) * (i == j)) for i in range(6) for j in range(6)))

# (d) covariance du réseau complet : U H(K+q) U† = H(K + R q) avec les relevés de l'octet ?
Us = [lift_diag(Hcube, pos, R)[0] for R in ROTS]
ok = True; okshift = True
for U, R in zip(Us, ROTS):
    for _ in range(5):
        q = rng.uniform(-np.pi, np.pi, 3); tf, tc = 0.6, 1.0
        ok &= np.allclose(U @ bloch(K + q, tf, tc) @ U.conj().T, bloch(K + R @ q, tf, tc), atol=1e-9)
print("covariance du réseau complet U_g H(K+q) U_g† = H(K+gq), 24 rotations, relevés de l'octet seuls :", ok)

# (e) saut R–R effectif au second ordre
def cube_hops(twist):
    A = np.zeros((8, 8))
    for a, x in enumerate(xs):
        for d in range(3):
            if x[d] == 0:
                y = list(x); y[d] = 1; b = xs.index(tuple(y)); A[a, b] = A[b, a] = ks(d, x) if twist else 1
    return A
Rset = [i for i, x in enumerate(xs) if sum(x) % 2 == 0]
for tw in (False, True):
    M = (cube_hops(tw) @ cube_hops(tw))[np.ix_(Rset, Rset)]
    print("A²|_R %s : %s ; valeurs propres %s" % ("tordu" if tw else "non tordu", M.astype(int).tolist(), np.round(np.linalg.eigvalsh(M), 6).tolist()))

# (f) amas ouvert : comptage des états dans le gap, pour tester la « fenêtre » 1/3 < t_f/t_c < 1
def cluster(N, tf, tc):
    sites = list(it.product(range(N), repeat=3)); idx = {s: i for i, s in enumerate(sites)}
    H = np.zeros((N ** 3,) * 2)
    for s in sites:
        for d in range(3):
            if s[d] + 1 < N:
                u = list(s); u[d] += 1; t = tf if s[d] % 2 == 0 else tc
                H[idx[s], idx[tuple(u)]] = H[idx[tuple(u)], idx[s]] = t * ks(d, s)
    return H
print("amas 8×8×8 (n = 4 octets/côté), t_c = 1 : nombre d'états avec |E| < 0,5·gap massif √3|t_f−t_c|")
for eps in [0.1, 0.2, 0.3, 0.5, 0.8, 1.25, 2.0]:
    w = np.abs(np.linalg.eigvalsh(cluster(8, eps, 1.0))); gap = np.sqrt(3) * abs(eps - 1)
    print("  t_f/t_c = %.2f : %d états (min |E| = %.2e)" % (eps, np.sum(w < 0.5 * gap), w.min()))
