"""FHT_21 / code_v30 : FHT hierarchiques (definition de V0_FHT_DefinitionAndEnumeration, A5 par composant).
Decomposition recursive depuis la realisation, et trois lectures homologiques :
  R0 : complexe engendre par la realisation (lecture plate actuelle appliquee telle quelle) ;
  R1 : substitution / join polyedral : |H| = U_{e^ in hote} join_{C in e^} |C|, |feuille| = point ;
  R2 : noeuds internes comme sommets : simplexes de l'hote sur les noeuds nu_C, plus cones nu_C * |C|_2."""
import itertools
from fht_tools import maximal, betti_z2

def decompose(V, E):
    """V : ensemble de feuilles, E : realisation (frozensets). Renvoie un arbre
    ('leaf', x) ou ('comp', V, E_top, [enfants]) ; E_top = aretes propres du composant racine."""
    V = frozenset(V)
    if len(V) == 1: return ('leaf', next(iter(V)))
    E = [frozenset(e) for e in E if frozenset(e) <= V]
    top = maximal(E)
    lower = [e for e in E if e not in top]
    # enfants = composantes connexes de (V, lower)
    parent = {x: x for x in V}
    def f(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for e in lower:
        l = list(e)
        for y in l[1:]: parent[f(y)] = f(l[0])
    blocks = {}
    for x in V: blocks.setdefault(f(x), set()).add(x)
    kids = [decompose(B, [e for e in lower if e <= B]) for B in blocks.values()]
    return ('comp', V, top, kids)

def leaves(t): return {t[1]} if t[0] == 'leaf' else set(t[1])
def depth(t): return 0 if t[0] == 'leaf' else 1 + max(depth(k) for k in t[3])
def host(t):
    """hote : enfants (indices) et aretes e^"""
    _, V, top, kids = t
    kl = [frozenset(leaves(k)) for k in kids]
    return kids, [frozenset(i for i, B in enumerate(kl) if B <= e) for e in top]

def R0(t):
    if t[0] == 'leaf': return [frozenset([t[1]])]
    out = []
    def rec(t):
        if t[0] == 'leaf': return
        out.extend(t[2]); [rec(k) for k in t[3]]
    rec(t); return maximal(out)

def join(cxs):
    out = [frozenset()]
    for C in cxs: out = [a | b for a in out for b in C]
    return out

def R1(t):
    if t[0] == 'leaf': return [frozenset([t[1]])]
    kids, he = host(t)
    sub = [R1(k) for k in kids]
    return maximal([f for e in he for f in join([sub[i] for i in e])])

def R2(t, tag=''):
    """noeuds internes nu_C (sauf racine) ; cone sur |C|_2 ; hote sur les noeuds"""
    if t[0] == 'leaf': return [frozenset([t[1]])], t[1]
    kids, he = host(t)
    facets, node = [], []
    for i, k in enumerate(kids):
        if k[0] == 'leaf': node.append(k[1]); continue
        sub, _ = R2(k, tag + '.%d' % i)
        nu = 'nu' + tag + '.%d' % i
        facets += [f | {nu} for f in sub]; node.append(nu)
    facets += [frozenset(node[i] for i in e) for e in he]
    return maximal(facets), None

def contracted_host(t):
    kids, he = host(t)
    return maximal([frozenset('c%d' % i for i in e) for e in he])
