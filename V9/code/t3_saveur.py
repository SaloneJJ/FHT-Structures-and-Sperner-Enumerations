"""V68 : dictionnaire residuels (G_e, G_nu) de S4_gen (rep 3 standard) -> colonnes PMNS fixees.
G_e abelien a spectre conjoint non degenere -> U_e ; G_nu : colonnes = vecteurs propres conjoints non degeneres.
Compare a TBM (2/3,1/6,1/6 | 1/3,1/3,1/3 | 0,1/2,1/2), TM1 (col1 TBM), TM2 (col2 TBM). JJ SALONE, FHT 2026"""
import itertools, numpy as np
S4 = list(itertools.permutations(range(4)))
def mul(a, b): return tuple(a[b[i]] for i in range(4))
def gen(gs):
    G = {tuple(range(4))}; fr = list(G)
    while fr:
        n = []
        for g in fr:
            for s in gs:
                h = mul(s, g)
                if h not in G: G.add(h); n.append(h)
        fr = n
    return sorted(G)
Bz = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float)
Bz = (Bz.T / np.linalg.norm(Bz, axis=1)).T  # base orthonormale de la rep standard
def rho(g):
    P = np.zeros((4, 4))
    for i in range(4): P[g[i], i] = 1
    return Bz @ P @ Bz.T
def joint_eig(G):
    """vecteurs propres conjoints ; renvoie liste des sous-espaces (matrices de colonnes)"""
    rng = np.random.default_rng(0)
    M = sum(rng.normal() * rho(g) for g in G) ; M = M + M.T if np.allclose(M, M.T) else M
    w, U = np.linalg.eig(M)
    spaces = []
    for val in np.unique(np.round(w, 8)):
        cols = U[:, np.abs(w - val) < 1e-6]
        spaces.append(cols)
    return spaces
def label(G):
    n = len(G); odd = sum(1 for g in G if sum(1 for i in range(4) for j in range(i+1, 4) if g[i] > g[j]) % 2)
    fixed = [i for i in range(4) if all(g[i] == i for g in G)]
    return f'|G|={n} impairs={odd} fixe={fixed} gens={G[1:3]}'
# sous-groupes abeliens (a conjugaison pres pour G_e, tous pour G_nu)
subs = {}
for k in range(1, 3):
    for gs in itertools.combinations(S4, k):
        G = tuple(gen(gs))
        if all(mul(a, b) == mul(b, a) for a in G for b in G): subs[G] = 1
subs = list(subs)
TBM = np.array([[2/3, 1/3, 0], [1/6, 1/3, 1/2], [1/6, 1/3, 1/2]])
def colname(c):
    c = np.sort(np.round(c, 4))
    for nm, t in [('TBM1/TM1 (2/3,1/6,1/6)', [1/6, 1/6, 2/3]), ('TBM2/TM2 (1/3,1/3,1/3)', [1/3, 1/3, 1/3]), ('TBM3 (0,1/2,1/2)', [0, .5, .5])]:
        if np.allclose(c, t, atol=1e-3): return nm
    return str(c)
C3 = tuple(gen([(0, 2, 3, 1)]))   # 3-cycle fixant 0
for Ge in [C3]:
    Se = joint_eig(Ge)
    if any(s.shape[1] > 1 for s in Se): continue
    Ue = np.hstack(Se)
    print('G_e =', label(Ge))
    seen = set()
    for Gn in subs:
        if len(Gn) == 1: continue
        Sn = joint_eig(Gn)
        cols = [s[:, 0] for s in Sn if s.shape[1] == 1]
        out = sorted(colname(np.abs(Ue.conj().T @ v) ** 2) for v in cols)
        key = (len(Gn), tuple(out), tuple(i for i in range(4) if all(g[i] == i for g in Gn)))
        if key in seen: continue
        seen.add(key)
        print('   G_nu', label(list(Gn)), '->', out if out else 'aucune colonne fixee')

# ---- regle FHT : stabilisateurs tordus (sgn) des chaines orientees de K4 (partie II, rep 3' = sgn x 3)
def sgn(g): return (-1) ** sum(1 for i in range(4) for j in range(i+1, 4) if g[i] > g[j])
def chain_sign(g, cell):
    img = [g[i] for i in cell]
    if set(img) != set(cell): return 0
    perm = [cell.index(x) for x in img]
    return (-1) ** sum(1 for i in range(len(perm)) for j in range(i+1, len(perm)) if perm[i] > perm[j])
def tw_stab(cell): return [g for g in S4 if chain_sign(g, list(cell)) * sgn(g) == 1]
print('\nRegle : residuel = stabilisateur tordu de la cellule marquee de K4 (sommet -> e, arete -> nu)')
Ge = tw_stab([0]); print('sommet g0 -> G_e =', label(Ge))
Ue = np.hstack(joint_eig(Ge))
for E in [(0, 1), (0, 2), (1, 2), (2, 3)]:
    Gn = tw_stab(list(E))
    cols = [s[:, 0] for s in joint_eig(Gn) if s.shape[1] == 1]
    print(f'arete {E} ({"drapeau" if 0 in E else "non incidente"}) -> G_nu = {Gn} ->', sorted(colname(np.abs(Ue.conj().T @ v) ** 2) for v in cols))
s13 = 0.0220; print('TM1 : sin2 th12 =', round(1 - (2/3) / (1 - s13), 4), ' (NuFit-6.0 : 0.307 +0.012 -0.011)')
