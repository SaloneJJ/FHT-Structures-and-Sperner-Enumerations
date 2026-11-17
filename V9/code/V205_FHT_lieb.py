# V205 — test de Lieb : le flux π minimise-t-il l'énergie des fermions à demi-remplissage ? Critères : V204.
import numpy as np, itertools as it
exec(open("V180_FHT_spin_cube_monopole.py").read().split("# ---------- octet (cube) ----------")[0])
rng = np.random.default_rng(205)
def egs(H):
    w = np.linalg.eigvalsh(H); return np.sum(w[:len(w)//2])
print("=== L1 : cube seul, toutes les configurations ℤ/2 (2^12 signes d'arêtes) ===")
cube = [np.array(c, float)*2-1 for c in it.product([0,1], repeat=3)]
ced = [(a, b) for a, b in it.combinations(range(8), 2) if np.sum(np.abs(cube[a]-cube[b])) == 2]
faces = [[i for i in range(8) if cube[i][ax] == s] for ax in range(3) for s in (1, -1)]
res = {}
for signs in it.product([1, -1], repeat=12):
    H = np.zeros((8, 8))
    for (a, b), s in zip(ced, signs): H[a, b] = H[b, a] = s
    npi = 0
    for f in faces:
        p = 1
        for (a, b), s in zip(ced, signs):
            if a in f and b in f: p *= s
        npi += (p == -1)
    e = round(egs(H), 6); res.setdefault(npi, set()).add(e)
for k in sorted(res): print("  %d faces à π : énergie(s) du fondamental = %s" % (k, sorted(res[k])))
cfaces = [outward_cycle(cube, f) for f in faces]
print("  flux U(1) uniforme sortant :", [(m, round(egs(build(cube, ced, cfaces, flux_list(6, 2*np.pi*m/6, m))), 6)) for m in range(6)])
def cubic(N, sign_fn, periodic=False):
    sites = list(it.product(range(N), repeat=3)); idx = {s: i for i, s in enumerate(sites)}
    H = np.zeros((N**3,)*2)
    for s in sites:
        for d in range(3):
            u = list(s); u[d] += 1
            if u[d] >= N:
                if not periodic: continue
                u[d] = 0
            H[idx[s], idx[tuple(u)]] += sign_fn(s, d); H[idx[tuple(u)], idx[s]] += sign_fn(s, d)
    return H
ks = lambda s, d: [1, (-1)**s[0], (-1)**(s[0]+s[1])][d]
print("=== L2 : amas ouverts ===")
for N in [4, 6]:
    e0 = egs(cubic(N, lambda s, d: 1)); epi = egs(cubic(N, ks))
    er = []
    for _ in range(200):
        R = {}
        def rs(s, d):
            k = (s, d)
            if k not in R: R[k] = rng.choice([1, -1])
            return R[k]
        er.append(egs(cubic(N, rs)))
    er = np.array(er)
    print("  %d³ : E(0) = %.4f ; E(π) = %.4f ; aléatoires : min %.4f, médiane %.4f ; nombre < E(π) : %d / 200 (par site : E(π)/N³ = %.4f, E(0)/N³ = %.4f)"
          % (N, e0, epi, er.min(), np.median(er), int(np.sum(er < epi - 1e-9)), epi/N**3, e0/N**3))
print("=== L3 : réseau périodique 8³ ===")
e0 = egs(cubic(8, lambda s, d: 1, True)); epi = egs(cubic(8, ks, True))
print("  E(0)/site = %.5f ; E(π)/site = %.5f" % (e0/512, epi/512))
print("=== L4 : coquille du germe K26 (26 sites, 24 carrés) ===")
shell = [np.array(p, float) for p in it.product([-1,0,1], repeat=3) if any(p)]
sed = [(a, b) for a, b in it.combinations(range(26), 2) if np.isclose(np.linalg.norm(shell[a]-shell[b]), 1)]
sq = []
for ax in range(3):
    for s in (1, -1):
        for u0, v0 in it.product([-1, 0], repeat=2):
            o = [k for k in range(3) if k != ax]; idq = []
            for du, dv in [(0,0),(1,0),(1,1),(0,1)]:
                p = np.zeros(3); p[ax] = s; p[o[0]] = u0+du; p[o[1]] = v0+dv
                idq.append(next(i for i, q in enumerate(shell) if np.allclose(q, p)))
            sq.append(outward_cycle(shell, idq))
for lab, phi, m in [("0", 0, 0), ("π par carré", np.pi, 12)]:
    H = build(shell, sed, sq, flux_list(24, phi, m)); w = np.linalg.eigvalsh(H)
    print("  %-12s : énergie des 13 plus bas niveaux = %.4f" % (lab, np.sum(w[:13])))
