"""AF21b -- Automorphismes par raffinement de couleurs (remplace aut_sgs de AF9 pour les gros groupes).

Diagnostic (14/11/2026) : AF21 ne terminait pas parce que aut_sgs (AF9) cherche les automorphismes par retour arriere
elague seulement par (type de sommet, co-degres).  Quand un candidat n'est pas dans l'orbite, la recherche doit epuiser
un arbre dont la taille croit avec |Aut| ; pour X_grand (|Aut| = 796 262 400) elle ne termine pas en temps raisonnable.
La version ci-dessous individualise un sommet puis raffine la coloration du graphe d'incidence sommets-facettes
(meme methode que fht_tools.aut_order, qui avait donne l'ordre attendu) : les branches mortes sont coupees tot.
Elle renvoie un ensemble fort de generateurs (un representant par point d'orbite non deja atteint) et l'ordre.
Les candidats deja dans l'orbite des generateurs trouves au meme niveau ne sont pas recherches (gain d'un facteur ~|orbite|).
"""
from collections import defaultdict, Counter


def aut_sgs_raffine(C):
    n = len(C.V); F = C.F; m = len(F); N = n + m
    adj = [[] for _ in range(N)]
    for j, f in enumerate(F):
        for v in f:
            adj[v].append(n + j); adj[n + j].append(v)
    Fset = C.Fset

    def refine_pair(c1, c2):
        while True:
            s1 = [(c1[x], tuple(sorted(c1[y] for y in adj[x]))) for x in range(N)]
            s2 = [(c2[x], tuple(sorted(c2[y] for y in adj[x]))) for x in range(N)]
            keys = sorted(set(s1) | set(s2)); k = {s: i for i, s in enumerate(keys)}
            n1 = [k[s] for s in s1]; n2 = [k[s] for s in s2]
            if len(set(n1)) == len(set(c1)) and len(set(n2)) == len(set(c2)):
                return n1, n2
            c1, c2 = n1, n2

    def search(c1, c2):
        """renvoie un automorphisme (liste) compatible avec les colorations, ou None"""
        c1, c2 = refine_pair(c1, c2)
        if Counter(c1) != Counter(c2): return None
        cells = defaultdict(list); cells2 = defaultdict(list)
        for x in range(n):
            cells[c1[x]].append(x); cells2[c2[x]].append(x)
        nontriv = [c for c in cells if len(cells[c]) > 1]
        if not nontriv:
            mp = {cells[c][0]: cells2[c][0] for c in cells}
            perm = [mp[v] for v in range(n)]
            return perm if all(frozenset(perm[v] for v in f) in Fset for f in F) else None
        c = min(nontriv, key=lambda c: (len(cells[c]), c)); v = cells[c][0]
        big = max(max(c1), max(c2)) + 1
        for w in cells2[c]:
            d1 = list(c1); d2 = list(c2); d1[v] = big; d2[w] = big
            r = search(d1, d2)
            if r is not None: return r
        return None

    def orbit(pt, gens):
        seen = {pt}; fr = [pt]
        while fr:
            nf = []
            for x in fr:
                for g in gens:
                    y = g[x]
                    if y not in seen: seen.add(y); nf.append(y)
            fr = nf
        return seen

    base0 = [0] * n + [1] * m
    col, _ = refine_pair(base0, base0)
    gens_all = []; order = 1
    for v in range(n):
        cand = [w for w in range(n) if col[w] == col[v]]
        if len(cand) == 1: continue
        big = max(col) + 1
        lev = []                                   # generateurs trouves a ce niveau (fixent les points de base precedents)
        orb = {v}
        for w in cand:
            if w in orb: continue
            d1 = list(col); d2 = list(col); d1[v] = big; d2[w] = big
            r = search(d1, d2)
            if r is not None:
                lev.append(tuple(r)); orb = orbit(v, lev)
        order *= len(orb); gens_all.extend(lev)
        d1 = list(col); d1[v] = big
        col, _ = refine_pair(d1, d1)
    return gens_all, order
