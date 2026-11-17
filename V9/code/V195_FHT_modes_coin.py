# V195 — E1 : états de coin du Dirac natif (réseau à flux π dimérisé, amas ouvert) ; E2 : ondelettes d'un nœud sous les rotations. Critères : V194.
import numpy as np, itertools as it
exec(open("V180_FHT_spin_cube_monopole.py").read().split("# ---------- octet (cube) ----------")[0])
def eta(d, x): return [1, (-1)**x[0], (-1)**(x[0]+x[1])][d]
def cluster(n, tf, tc):
    N = 2*n; sites = list(it.product(range(N), repeat=3)); idx = {s: i for i, s in enumerate(sites)}
    H = np.zeros((len(sites),)*2)
    for s in sites:
        for d in range(3):
            if s[d] + 1 < N:
                t = tf if s[d] % 2 == 0 else tc
                u = list(s); u[d] += 1; H[idx[s], idx[tuple(u)]] = H[idx[tuple(u)], idx[s]] = t*eta(d, s)
    return H, sites
print("=== E1 : états de coin ===")
for ratio in [0.5, 0.3, 2.0]:
    for n in [2, 3, 4]:
        H, sites = cluster(n, ratio, 1.0); w, V = np.linalg.eigh(H); o = np.argsort(np.abs(w)); N = 2*n
        e8 = np.abs(w[o[:8]]); e9 = np.abs(w[o[8]])
        corner_oct = np.array([all(c in (0, 1) or c in (N-2, N-1) for c in s) for s in sites])
        wc = np.sum(np.abs(V[:, o[:8]])**2 * corner_oct[:, None], axis=0)
        par = np.array([sum(s) % 2 for s in sites])
        pe = np.sum(np.abs(V[:, o[:8]])**2 * (par == 0)[:, None], axis=0)
        print(" t_f/t_c=%.1f n=%d (%3d sites) : 8 plus petits |E| max = %.2e, 9e = %.3f ; poids aux octets de coin (min) = %.3f ; poids pair par état = %s"
              % (ratio, n, len(sites), e8.max(), e9, wc.min(), np.round(pe, 2).tolist()))
# loi de décroissance avec n à ratio fixé
print("=== E1-b : |E| des modes de coin contre n ===")
for ratio in [0.5, 0.7]:
    out = []
    for n in [2, 3, 4, 5]:
        H, _ = cluster(n, ratio, 1.0); w = np.linalg.eigvalsh(H); out.append(np.sort(np.abs(w))[7])
    sl = [np.log(out[k+1]/out[k])/np.log(ratio) for k in range(len(out)-1)]
    print(" t_f/t_c=%.1f : max |E| des 8 états, n = 2..5 : %s ; exposant a (|E| ∝ ε^{a n}) entre n successifs : %s" % (ratio, ["%.2e" % e for e in out], np.round(sl, 2).tolist()))
# chiralité de sous-réseau : dans quel coin vit chaque état (base localisée : diagonaliser la projection sur les 8 états du coin)
H, sites = cluster(3, 0.5, 1.0); w, V = np.linalg.eigh(H); o = np.argsort(np.abs(w)); P = V[:, o[:8]]
N = 6; corners = [tuple(c) for c in it.product([0, N-1], repeat=3)]
print("=== E1-c : chaque coin porte-t-il un état, sur quel sous-réseau ? (n = 3, ratio 0,5) ===")
for c in corners:
    mask = np.array([max(abs(s[k]-c[k]) for k in range(3)) <= 1 for s in sites])   # octet du coin
    Pc = P[mask]; ev = np.linalg.eigvalsh(Pc.T @ Pc)
    # sous-réseau du site de coin lui-même
    print("  coin %s (parité du site de coin %d) : poids max d'un état projeté sur l'octet du coin = %.3f" % (c, sum(c) % 2, ev.max()))
print("=== E2 : ondelettes d'un nœud (7 = 8 − moyenne) sous O ===")
cube = [np.array(c, float)*2-1 for c in it.product([0,1], repeat=3)]
key = {tuple(p): i for i, p in enumerate(cube)}
cls = {}
for R in ROT:
    P_ = np.zeros((8, 8))
    for i, p in enumerate(cube): P_[key[tuple(R @ p)], i] = 1
    chi = np.trace(P_) - 1; cls.setdefault(round(np.trace(R)), set()).add(round(chi))
print("  caractère des 7 ondelettes par classe (trace de R : 3 = E, 0 = C3, −1 = C2 et C2', 1 = C4) :", {k: sorted(v) for k, v in cls.items()})
print("  attendu A2 + T1 + T2 : E 7, C3 1, C2 face −1, C4 −1, C2' −1")
ced = [(a, b) for a, b in it.combinations(range(8), 2) if np.sum(np.abs(cube[a]-cube[b])) == 2]
cfaces = [outward_cycle(cube, [i for i in range(8) if cube[i][ax] == s]) for ax in range(3) for s in (1, -1)]
Ht = build(cube, ced, cfaces, flux_list(6, np.pi, 3)); Us = [lift(Ht, cube, R)[0] for R in ROT]
u = np.ones(8)/np.sqrt(8)
print("  tordu : la moyenne est-elle invariante (à une phase près) pour les 24 relevés ?", all(abs(abs(np.vdot(u, U @ u)) - 1) < 1e-9 for U in Us))
