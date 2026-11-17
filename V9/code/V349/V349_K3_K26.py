"""V349 K3 : (a) permutations impaires de S4_L, S4_R dans Aut des représentants ; (b) orientation des faces du tétraèdre x (zone L-X induite) ;
(c) liens des coins L_x, R_x ; (d) [ajouté après le message de V350, avant tout calcul de K3] les boucles des faces x_i x_j x_k sont-elles
bordées dans K26 entier ?  Critères : 01_synthese/V349 §2.  numpy seul.   JJ SALONE, FHT 2026"""
import json, itertools as it, numpy as np
from collections import defaultdict
def load(f):
    d = json.load(open(f)); return [frozenset(x) for x in (d if isinstance(d, list) else d.get('facettes', d.get('facets')))]
# appariements : p0 = {01|23}, p1 = {02|13}, p2 = {03|12} (convention à vérifier : on teste les 6 bijections induites possibles)
PAIRS = [frozenset([frozenset([0, 1]), frozenset([2, 3])]), frozenset([frozenset([0, 2]), frozenset([1, 3])]), frozenset([frozenset([0, 3]), frozenset([1, 2])])]
def induced_p(pi, order):
    """action de pi (sur 0..3) sur les appariements, pour une étiquette order : p_i -> PAIRS[order[i]]"""
    lab = {PAIRS[order[i]]: i for i in range(3)}
    return {i: lab[frozenset(frozenset(pi[a] for a in pr) for pr in PAIRS[order[i]])] for i in range(3)}
def perm_vertex(side, pi, order):
    pp = induced_p(pi, order)
    def m(v):
        if not v.startswith(side + '_'): return v
        n = v[2:]
        if n.startswith('x'): return f'{side}_x{pi[int(n[1])]}'
        if n.startswith('p'): return f'{side}_p{pp[int(n[1])]}'
        return v
    return m
def is_aut(F, m): S = set(F); return all(frozenset(map(m, f)) in S for f in F)
def faces(F):
    S = set()
    for f in F:
        for k in range(1, len(f) + 1): S.update(it.combinations(sorted(f), k))
    return S
def betti(F):
    S = faces(F); bd = defaultdict(list)
    for c in S: bd[len(c) - 1].append(c)
    idx = {d: {c: i for i, c in enumerate(sorted(L))} for d, L in bd.items()}; rk = {0: 0}
    for d in sorted(bd):
        if d == 0: continue
        M = np.zeros((len(bd[d-1]), len(bd[d])))
        for j, c in enumerate(sorted(bd[d])):
            for i in range(len(c)): M[idx[d-1][c[:i] + c[i+1:]], j] = (-1)**i
        rk[d] = np.linalg.matrix_rank(M)
    top = max(bd); b = [len(bd[d]) - rk[d] - rk.get(d+1, 0) for d in range(top + 1)]; b[0] -= 1
    return b
reps = (('K26 V37', 'V37_FHT_K26volmin.json'), ('K26vol V36', 'V36_FHT_K26vol.json'), ('K26* V32', 'V32_FHT_K26etoile.json'))
tr = {0: 1, 1: 0, 2: 2, 3: 3}; cyc = {0: 1, 1: 2, 2: 0, 3: 3}; dbl = {0: 1, 1: 0, 2: 3, 3: 2}
for nm, fn in reps:
    F = load(fn); print(f'===== {nm} : {len(F)} facettes')
    # étiquetage des p : trouver un ordre pour lequel le 3-cycle (pair) est un automorphisme
    for side in 'LR':
        good = [o for o in it.permutations(range(3)) if is_aut(F, perm_vertex(side, cyc, o)) and is_aut(F, perm_vertex(side, dbl, o))]
        o = good[0] if good else None
        print(f'  côté {side} : étiquetages des p compatibles avec A4 : {len(good)} ; ', end='')
        if o is None: print('aucun'); continue
        print(f'transposition x0<->x1 (impaire) automorphisme : {is_aut(F, perm_vertex(side, tr, o))} ; 3-cycle : True ; double transposition : True')
    # (c) liens des coins
    for v in ('L_x0', 'R_x0', 'L_x2', 'R_x3'):
        L = [f - {v} for f in F if v in f]; print(f'  lien de {v} : {len(L)} facettes, Betti réduits {betti(L)}')
    # (d) faces x_i x_j x_k : bordées dans K26 entier ? (cycle triangle = bord d'un 2-simplexe ou bord dans le lien d'un cône)
    for side in 'LR':
        tri = [frozenset(f'{side}_x{i}' for i in t) for t in it.combinations(range(4), 3)]
        filled = [any(t <= f for f in F) for t in tri]
        cones = sorted({w for t in tri for f in F if t <= f or False for w in ()})
        apex = sorted({w for t in tri for f in F for w in f if all(frozenset({a, b, w}) <= g for a, b in it.combinations(sorted(t), 2) for g in [next((h for h in F if frozenset({a, b, w}) <= h), frozenset())])} - set().union(*tri))
        print(f'  côté {side} : triangles x_i x_j x_k faces de K26 : {filled} ; sommets w tels que les 3 triangles (x_a x_b w) existent pour chaque face : {apex[:8]}{"..." if len(apex) > 8 else ""} ({len(apex)})')

# (b) motif de flux « +π/2 sur chaque face sortante » du tétraèdre x et image par la transposition x0<->x1 (lemme « Opposite faces »)
def outward(perm_or_id):
    # face opposée à x_i, orientée sortante : (-1)^i [faces ordonnées] ; on renvoie l'orientation induite par pi comparée à la sortante
    out = []
    for i in range(4):
        fc = [j for j in range(4) if j != i]; img = [perm_or_id[j] for j in fc]
        # signe de la permutation qui ordonne img, et la face image est opposée à pi(i)
        sg = np.linalg.det(np.eye(3)[np.argsort(img)]); ii = perm_or_id[i]
        out.append(round((-1)**i * sg * (-1)**ii))
    return out
print('(b) image des 4 faces sortantes par x0<->x1 (+1 : sortante, -1 : rentrante) :', outward(tr), ' ; par le 3-cycle :', outward(cyc))
