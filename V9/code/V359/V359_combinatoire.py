"""V359 (relecture indépendante de V347-V350) : combinatoire de K26 recalculée sans importer le code relu.
Rangs sur GF(p) (p = 32003) par élimination ; homologie réduite ; zones induites ; caractère de S4_L sur H1(L-X) ;
cône de R_z ; aires des 3-cycles du 1-squelette ; secteurs I/II ; automorphismes impairs ; liens des coins.
Usage : python3 V359_combinatoire.py <dossier des JSON>   (numpy seul)."""
import json, sys, itertools as it, numpy as np
from collections import Counter, defaultdict
D = sys.argv[1] if len(sys.argv) > 1 else '.'
P = 32003
def load(f):
    d = json.load(open(f"{D}/{f}")); d = d['facettes'] if isinstance(d, dict) else d
    return [frozenset(x) for x in d]
REP = {'V37': 'V37_FHT_K26volmin.json', 'V36': 'V36_FHT_K26vol.json', 'V32': 'V32_FHT_K26etoile.json'}
def closure(F):
    S = set()
    for f in F:
        s = sorted(f)
        for k in range(1, len(s)+1):
            S.update(it.combinations(s, k))
    return S
def rank_modp(rows, ncols):
    """rows : liste de dict {col: val} ; rang sur GF(P) par élimination creuse."""
    piv = {}
    r = 0
    for row in rows:
        row = {c: v % P for c, v in row.items() if v % P}
        while row:
            c = min(row)
            if c in piv:
                prow = piv[c]; f = row[c]
                for cc, vv in prow.items():
                    row[cc] = (row.get(cc, 0) - f*vv) % P
                    if row[cc] == 0: del row[cc]
            else:
                inv = pow(row[c], P-2, P); row = {cc: vv*inv % P for cc, vv in row.items()}
                piv[c] = row; r += 1; break
    return r
def betti(F):
    """Betti réduits (degrés 0..dim) sur GF(P)."""
    S = closure(F); by = defaultdict(list)
    for c in S: by[len(c)-1].append(c)
    top = max(by); idx = {d: {c: i for i, c in enumerate(sorted(by[d]))} for d in by}
    rk = {0: 1}  # augmentation (réduite)
    for d in range(1, top+1):
        rows = []
        for c in by[d]:
            rows.append({idx[d-1][c[:i]+c[i+1:]]: (-1)**i for i in range(len(c))})
        rk[d] = rank_modp(rows, len(by[d-1]))
    return [len(by[d]) - rk[d] - rk.get(d+1, 0) for d in range(top+1)]
def induced(F, V):
    R = set(frozenset(f & V) for f in F if f & V)
    return [f for f in R if not any(f < g for g in R)]
def zone(v):
    if v.startswith('g'): return 'G'
    s, n = v.split('_', 1)
    if n == 'z' or n.startswith('x'): return s+'X'
    if n == 'f' or n.startswith('p'): return s+'P'
    return s+n
out = {}
print("== 1. Betti réduits (GF(32003)) des zones induites et de K26 entier")
for nm, f in REP.items():
    F = load(f); V = sorted(set().union(*F)); Z = defaultdict(set)
    for v in V: Z[zone(v)].add(v)
    if nm != 'V32':
        line = {z: betti(induced(F, Z[z])) for z in ('G', 'LX', 'RX', 'LP', 'RP')}
        print(f"  {nm} : {len(F)} facettes ; zones : " + " ; ".join(f"{z} {b}" for z, b in line.items()))
    print(f"  {nm} : Betti réduits de K entier : {betti(F)}")
F = load(REP['V37']); S37 = closure(F)
E = {c for c in S37 if len(c) == 2}; T = {c for c in S37 if len(c) == 3}; V = sorted(set().union(*F))
print("\n== 2. Caractère de S4_L sur H1(L-X induit) (K5 sans triangle dans V37 ?)")
LX = ['L_z', 'L_x0', 'L_x1', 'L_x2', 'L_x3']
eLX = [e for e in E if set(e) <= set(LX)]; tLX = [t for t in T if set(t) <= set(LX)]
print(f"  arêtes {len(eLX)}, triangles {len(tLX)}")
# cycles : noyau de ∂1 sur les arêtes de K5 ; action de S4 (permutation des x) avec signes d'orientation
edges = sorted(tuple(sorted(e)) for e in eLX); ei = {e: i for i, e in enumerate(edges)}
vi = {v: i for i, v in enumerate(LX)}
d1 = np.zeros((5, len(edges)))
for j, (a, b) in enumerate(edges): d1[vi[a], j] = -1; d1[vi[b], j] = 1
u, s, vh = np.linalg.svd(d1); Zb = vh[np.sum(s > 1e-9):].T      # base orthonormée des cycles (10 x 6)
def act(p):  # p : permutation de range(4) -> matrice sur les arêtes
    mp = {'L_z': 'L_z', **{f'L_x{i}': f'L_x{p[i]}' for i in range(4)}}
    M = np.zeros((len(edges), len(edges)))
    for j, (a, b) in enumerate(edges):
        a2, b2 = mp[a], mp[b]; e2 = tuple(sorted((a2, b2))); M[ei[e2], j] = 1 if (a2, b2) == e2 else -1
    return M
cls = {'1': (0, 1, 2, 3), '(12)': (1, 0, 2, 3), '(12)(34)': (1, 0, 3, 2), '(123)': (1, 2, 0, 3), '(1234)': (1, 2, 3, 0)}
chi = {k: round(float(np.trace(Zb.T @ act(p) @ Zb)), 6) for k, p in cls.items()}
print(f"  caractère sur H1 : {chi}   (3 ⊕ 3' = (6, 0, -2, 0, 0))")
print("  remarque : H1(K5) = C1(K4 orienté) = Λ²(perm4) = Λ²(1 ⊕ 3) = 3 ⊕ 3' pour TOUT S4 agissant sur K5 avec un sommet fixe")
print("\n== 3. Cône de R_z (V37)")
nb = {v for e in E for v in e if 'R_z' in e} - {'R_z'}
lk_edges = {tuple(sorted(set(t) - {'R_z'})) for t in T if 'R_z' in t}
other = {tuple(sorted(e)) for e in E if 'R_z' not in e}
print(f"  voisins de R_z : {len(nb)} / {len(V)-1} ; arêtes hors R_z : {len(other)} ; dans le lien : {len(other & lk_edges)} ; manquantes : {sorted(other - lk_edges)}")
for apex in V:
    nbv = {v for e in E for v in e if apex in e} - {apex}
    lk = {tuple(sorted(set(t) - {apex})) for t in T if apex in t}
    oth = {tuple(sorted(e)) for e in E if apex not in e}
    if len(nbv) >= 24: print(f"  sommet {apex} : voisins {len(nbv)}, arêtes du lien {len(oth & lk)}/{len(oth)}")
# 3-cycles du 1-squelette : remplis ou non ; aire minimale par recherche bornée (≤ 3 triangles)
adj = defaultdict(set)
for a, b in E: adj[a].add(b); adj[b].add(a)
cyc = [c for c in it.combinations(V, 3) if c[1] in adj[c[0]] and c[2] in adj[c[0]] and c[2] in adj[c[1]]]
Ts = {tuple(sorted(t)) for t in T}
unf = [c for c in cyc if c not in Ts]
def area3(c):
    """aire 3 ? on cherche un sommet v tel que les 3 triangles v-arête soient des triangles de K (cône) ;
    sinon on signale (l'aire 2 est exclue par parité)."""
    a, b, cc = c
    return [v for v in V if v not in c and all(tuple(sorted((v,)+e)) in Ts for e in ((a, b), (b, cc), (a, cc)))]
nocone = [c for c in unf if not area3(c)]
apexes = Counter(v for c in unf for v in area3(c))
print(f"  3-cycles du 1-squelette : {len(cyc)} ; remplis : {len(cyc)-len(unf)} ; non remplis : {len(unf)} ; non remplis sans cône d'un sommet : {len(nocone)}")
print(f"  sommets-cônes disponibles (nombre de 3-cycles non remplis qu'ils remplissent) : {apexes.most_common(6)}")
print(f"  3-cycles non remplis contenant l'arête L_w–R_f : {[c for c in unf if 'L_w' in c and 'R_f' in c]}")
print("  (parité : une 2-chaîne ℤ/2 de k triangles a un bord de longueur ≡ k mod 2 ⇒ aire d'un 3-cycle non rempli ≥ 3)")
print("\n== 4. Secteurs I / II (nombre de g par facette) et zones")
ng = Counter(sum(v.startswith('g') for v in f) for f in F); print(f"  facettes par nombre de g : {dict(sorted(ng.items()))}")
for k in (0, 1, 2):
    zz = Counter(zone(v) for f in F if sum(v.startswith('g') for v in f) == k for v in f if not v.startswith('g'))
    print(f"   {k} g : zones présentes {dict(sorted(zz.items()))}")
# vecteur de Fiedler du 1-squelette pondéré (multiplicité des arêtes dans les facettes)
W = Counter()
for f in F:
    for a, b in it.combinations(sorted(f), 2): W[(a, b)] += 1
ix = {v: i for i, v in enumerate(V)}; Aj = np.zeros((len(V),)*2)
for (a, b), w in W.items(): Aj[ix[a], ix[b]] = Aj[ix[b], ix[a]] = w
ev, U = np.linalg.eigh(np.diag(Aj.sum(1)) - Aj)
fi = U[:, 1]; pos = sorted(v for v in V if fi[ix[v]] > 0); neg = sorted(v for v in V if fi[ix[v]] <= 0)
print(f"  laplacien pondéré : λ2 = {ev[1]:.2f}, λ3 = {ev[2]:.2f} ; Fiedler + : {pos}\n                                          Fiedler − : {neg}")
# Fiedler en grandeur : séparation nette ?
print("  composantes de Fiedler (triées) : " + " ".join(f"{v}:{fi[ix[v]]:+.2f}" for v in sorted(V, key=lambda v: fi[ix[v]])))
print("\n== 5. Transposition x0 ↔ x1 de S4_L (et de S4_R) avec action induite sur les p : automorphisme ?")
for nm in ('V37', 'V36', 'V32'):
    FF = set(load(REP[nm])); res = []
    for side in 'LR':
        found = []
        for pp in it.permutations(range(3)):
            m = {f'{side}_x0': f'{side}_x1', f'{side}_x1': f'{side}_x0', **{f'{side}_p{i}': f'{side}_p{pp[i]}' for i in range(3)}}
            if {frozenset(m.get(v, v) for v in f) for f in FF} == FF: found.append(pp)
        res.append((side, found))
    print(f"  {nm} : " + " ; ".join(f"côté {s} : permutations des p qui conviennent {fd}" for s, fd in res))
print("\n== 6. Liens des coins (Betti réduits) et homologie des secteurs (degré de I contre II)")
for nm in ('V37', 'V36', 'V32'):
    FF = load(REP[nm])
    for v in ('R_x0', 'L_x0'):
        L = [f - {v} for f in FF if v in f]
        print(f"  {nm} : lk({v}) : {len(L)} facettes, Betti réduits {betti(L)}")
FI = [f for f in F if sum(v.startswith('g') for v in f) == 1]; FII = [f for f in F if sum(v.startswith('g') for v in f) == 2]
print(f"  V37 : sous-complexe des facettes à 1 g : Betti {betti(FI)} ; à 2 g : Betti {betti(FII)}")
