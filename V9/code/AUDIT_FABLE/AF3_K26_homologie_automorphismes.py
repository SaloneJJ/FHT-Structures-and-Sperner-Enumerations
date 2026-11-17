#!/usr/bin/env python3
"""AUDIT FABLE -- A2, registre T3-T4 : K26 reconstruit a partir de la seule Definition de l'article (pieces A et B, I = A^L * B^R * G0,
II = B^L * A^R * K4), puis : nombre d'aretes (716), de simplexes (47 167), nombres de Betti reduits (12, 36, 36, 12 en degres 3-6),
ordre du groupe d'automorphismes (13 824 = 24^3), homologie des pieces. Meme controle sur les representants deposes dans 04_objets
(V37 = K26 « volume minimal », 542 facettes, 7 489 simplexes ; V36 ; V32).
Code independant : aucune ligne des scripts de CalculCodes n'est reutilisee. Homologie sur F_2 et sur F_p (p = 65521) par reduction de colonnes.
Automorphismes : chaine de stabilisateurs, recherche d'existence par retour arriere (elagage par types de sommets et co-degres de paires).
Usage : python3 AF3_K26_homologie_automorphismes.py [chemin de 04_objets]       (python3 pur ; ~2 a 4 min)
"""
import itertools, json, os, sys, time
from collections import Counter, defaultdict

# ---------------------------------------------------------------- construction d'apres la Definition
def piece_A(side):
    v = lambda n: f'{side}_{n}'
    Z = [v('z')] + [v(f'x{i}') for i in range(4)]
    P = {0: ((0, 1), (2, 3)), 1: ((0, 2), (1, 3)), 2: ((0, 3), (1, 2))}       # les trois appariements de {0,1,2,3}
    F = [frozenset(Z)]
    F += [frozenset([v('f'), v(f'p{i}'), v(f'p{j}')]) for i, j in itertools.combinations(range(3), 2)]
    F += [frozenset(v(f'p{i}') for i in range(3))]
    for k, (ab, cd) in P.items():
        for (a, b) in (ab, cd):
            F.append(frozenset([v(f'p{k}'), v('z'), v(f'x{a}'), v(f'x{b}')]))
    F.append(frozenset([v('q')]))
    return F

def piece_B(side):
    v = lambda n: f'{side}_{n}'
    Z = [v('z')] + [v(f'x{i}') for i in range(4)]
    return [frozenset(Z)] + [frozenset([v('w'), z]) for z in Z]

def join(*families):
    return [frozenset().union(*c) for c in itertools.product(*families)]

def maximal(faces):
    faces = sorted(set(faces), key=len, reverse=True); out = []
    for f in faces:
        if not any(f < g for g in out): out.append(f)
    return out

def K26_definition():
    G0 = [frozenset([f'g{i}']) for i in range(4)]
    K4 = [frozenset([f'g{i}', f'g{j}']) for i, j in itertools.combinations(range(4), 2)]
    I = join(piece_A('L'), piece_B('R'), G0); II = join(piece_B('L'), piece_A('R'), K4)
    return maximal(I + II)

# ---------------------------------------------------------------- homologie reduite sur F_p
def simplices(facets):
    S = set()
    for f in facets:
        f = sorted(f)
        for r in range(1, len(f) + 1):
            S.update(itertools.combinations(f, r))
    by = defaultdict(list)
    for s in S: by[len(s) - 1].append(s)
    return {d: sorted(v) for d, v in by.items()}

def rank_boundary(cols, p):
    """cols : liste de dict {ligne: coeff} ; rang sur F_p par reduction (pivot = plus grand indice de ligne)"""
    piv = {}; rank = 0
    for c in cols:
        c = dict(c)
        while c:
            l = max(c)
            if l not in piv: break
            pc, inv = piv[l]; fac = (c[l] * inv) % p
            for k, val in pc.items():
                nv = (c.get(k, 0) - fac * val) % p
                if nv: c[k] = nv
                else: c.pop(k, None)
        if c:
            l = max(c); piv[l] = (c, pow(c[l], p - 2, p)); rank += 1
    return rank

def betti(facets, p=2):
    S = simplices(facets); dmax = max(S); idx = {d: {s: i for i, s in enumerate(S[d])} for d in S}
    n = {d: len(S[d]) for d in S}; n[-1] = 1
    rk = {0: 1}                                      # bord augmente C_0 -> C_{-1} : rang 1
    for d in range(1, dmax + 1):
        cols = []
        for s in S[d]:
            col = {}
            for i in range(len(s)):
                col[idx[d - 1][s[:i] + s[i + 1:]]] = (1 if i % 2 == 0 else p - 1) % p
            cols.append(col)
        rk[d] = rank_boundary(cols, p)
    rk[dmax + 1] = 0
    return [n[d] - rk[d] - rk[d + 1] for d in range(0, dmax + 1)], sum(n[d] for d in S)

# ---------------------------------------------------------------- automorphismes
def aut_order(facets):
    V = sorted(set().union(*facets)); nv = len(V); ix = {v: i for i, v in enumerate(V)}
    F = [frozenset(ix[v] for v in f) for f in facets]; Fset = set(F)
    vtype = [tuple(sorted(Counter(len(f) for f in F if v in f).items())) for v in range(nv)]
    cod = [[tuple(sorted(Counter(len(f) for f in F if u in f and v in f).items())) for v in range(nv)] for u in range(nv)]
    def extend(img, order, pos):
        if pos == len(order):
            return all(frozenset(img[v] for v in f) in Fset for f in F)
        u = order[pos]; used = set(img.values())
        for y in range(nv):
            if y in used or vtype[y] != vtype[u]: continue
            if all(cod[u][w] == cod[y][img[w]] for w in img):
                img[u] = y
                if extend(img, order, pos + 1): return True
                del img[u]
        return False
    order_total = 1; fixed = {}
    base = sorted(range(nv), key=lambda v: (vtype[v], v))
    for b in base:
        if b in fixed: continue
        rest = [v for v in base if v not in fixed and v != b]
        orbit = 0
        for y in range(nv):
            if y in fixed.values() or vtype[y] != vtype[b]: continue
            if not all(cod[b][w] == cod[y][fixed[w]] for w in fixed): continue
            img = dict(fixed); img[b] = y
            if extend(img, rest, 0): orbit += 1
        order_total *= orbit; fixed[b] = b
    return order_total

def report(name, facets, expect):
    t0 = time.time()
    V = set().union(*facets)
    b2, ns = betti(facets, 2); bp, _ = betti(facets, 65521)
    a = aut_order(facets)
    line = f'{name:28s} sommets {len(V)}, aretes {len(facets)}, simplexes {ns}, Betti F2 {b2}, Betti F_65521 {bp}, |Aut| = {a}  ({time.time()-t0:.0f} s)'
    trim = lambda b: [x for i, x in enumerate(b) if i < 7 or x]            # zeros de queue ignores
    ok = all(expect.get(k) in (None, val) for k, val in (('edges', len(facets)), ('simplices', ns), ('betti', trim(bp)), ('aut', a))) and b2 == bp
    print(line, ' ->', 'CONFORME' if ok else 'ECART', flush=True)
    return ok

if __name__ == '__main__':
    obj = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '04_objets')
    print('Pieces (Proposition « Homology of the pieces ») :')
    bA, nA = betti(piece_A('L'), 65521); bB, nB = betti(piece_B('L'), 65521)
    print(f'  A : 10 sommets, {len(piece_A("L"))} faces maximales, Betti reduits {bA}  (attendu 1, 2, 1 en degres 0, 1, 2) ; |Aut(A)| = {aut_order(piece_A("L"))} (attendu 24)')
    print(f'  B : 6 sommets, {len(piece_B("L"))} faces maximales, Betti reduits {bB}  (attendu 4 en degre 1) ; |Aut(B)| = {aut_order(piece_B("L"))} (attendu 120 : S_5 sur Z)')
    B26 = [0, 0, 0, 12, 36, 36, 12]
    K = K26_definition()
    ok = report('K26 = K^(4) (Definition A2)', K, dict(edges=716, simplices=47167, betti=B26, aut=13824))
    for fn, exp in (('V37_FHT_K26volmin.json', dict(edges=542, simplices=7489, aut=13824)), ('V36_FHT_K26vol.json', dict(aut=13824)), ('V32_FHT_K26etoile.json', dict(aut=13824))):
        path = os.path.join(obj, fn)
        if not os.path.exists(path): print('  (absent :', path, ')'); continue
        d = json.load(open(path)); fac = d if isinstance(d, list) else d.get('facets', d.get('facettes'))
        fac = [frozenset(f) for f in fac]
        exp = dict(exp, betti=B26)
        ok &= report(fn, fac, exp)
    print('VERDICT :', 'conforme' if ok else 'ECART')
