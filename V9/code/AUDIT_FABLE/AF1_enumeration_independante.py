#!/usr/bin/env python3
"""AUDIT FABLE -- A1, registre F1-F2 : enumeration independante des FHT, ecrite a partir de la seule Definition 2.1 (A1)-(A5) de l'article A1.
Aucune ligne de FHT_Enumeration.py n'est reutilisee.

Methode 1 (force brute au niveau des axiomes, v <= 4) : on parcourt TOUTES les familles non vides de parties de taille >= 2 de V
  (2^11 - 1 familles pour v = 4) et on garde celles qui sont la realisation d'une FHT ; on compte a isomorphisme pres.
  Test de realisation (deduit des axiomes, voir le rapport AUDIT_FABLE_A1.md, section 2) :
    - les aretes du sommet sont les aretes maximales pour l'inclusion ; elles doivent former un clutter connexe couvrant V ;
    - les enfants de la racine sont les composantes connexes de la famille des aretes non maximales (les sommets non couverts
      sont des feuilles) ; chaque arete maximale qui rencontre un enfant doit le contenir strictement (A4) ;
    - chaque enfant porte recursivement une FHT.
Methode 2 (decomposition recursive, v <= 5) : partition de V en k >= 2 enfants, FHT sur chaque enfant, clutter connexe sur les enfants
  dont chaque arete contient au moins deux enfants ; comptage a isomorphisme pres par forme canonique sous S_v.
Attendu (A1, resume) : 1, 1, 4, 23, 231.
Usage : python3 AF1_enumeration_independante.py            (quelques secondes ; v = 5 : ~1 min)
"""
import itertools, sys, time

def connected(edges, verts):
    verts = set(verts)
    if not verts: return True
    comp = {next(iter(verts))}
    ch = True
    while ch:
        ch = False
        for e in edges:
            if comp & e and not e <= comp:
                comp |= e; ch = True
    return comp == verts

def components(edges):
    """composantes connexes (ensembles de sommets) de l'hypergraphe des aretes donnees"""
    comps = []
    for e in edges:
        e = set(e); merged = [c for c in comps if c & e]
        for c in merged:
            e |= c; comps.remove(c)
        comps.append(e)
    return [frozenset(c) for c in comps]

def is_realisation(E, V):
    """E : ensemble de frozensets (taille >= 2) ; V : frozenset. Vrai ssi E est la realisation d'une FHT sur V (|V| >= 2)."""
    if not E: return False
    top = [e for e in E if not any(e < f for f in E)]
    if frozenset().union(*top) != V or not connected(top, V): return False
    low = [e for e in E if e not in top]
    for C in components(low):
        sub = frozenset(e for e in low if e <= C)
        if not is_realisation(sub, C): return False
        for e in top:
            if e & C and not C < e: return False      # (A4) : inclusion stricte
    return True

def canon(E, v):
    best = None
    for p in itertools.permutations(range(v)):
        key = tuple(sorted(tuple(sorted(p[x] for x in e)) for e in E))
        if best is None or key < best: best = key
    return best

def brute(v):
    V = frozenset(range(v))
    if v == 1: return 1
    subs = [frozenset(s) for r in range(2, v + 1) for s in itertools.combinations(range(v), r)]
    seen = set()
    for mask in range(1, 1 << len(subs)):
        E = frozenset(subs[i] for i in range(len(subs)) if mask >> i & 1)
        if is_realisation(E, V): seen.add(canon(E, v))
    return len(seen)

# ---------------------------------------------------------------- methode 2
def set_partitions(lst):
    if not lst: yield []; return
    a, rest = lst[0], lst[1:]
    for p in set_partitions(rest):
        yield [[a]] + p
        for i in range(len(p)):
            yield p[:i] + [[a] + p[i]] + p[i + 1:]

_cc = {}
def connected_clutters(k):
    """clutters connexes couvrants sur {0..k-1}, aretes de taille >= 2 (k >= 2), par recherche en profondeur sur les antichaines"""
    if k in _cc: return _cc[k]
    subs = [frozenset(s) for r in range(2, k + 1) for s in itertools.combinations(range(k), r)]
    out = []
    def rec(i, cur):
        if i == len(subs):
            if cur and frozenset().union(*cur) == frozenset(range(k)) and connected(cur, range(k)): out.append(list(cur))
            return
        rec(i + 1, cur)
        s = subs[i]
        if not any(s <= c or c <= s for c in cur):
            cur.append(s); rec(i + 1, cur); cur.pop()
    rec(0, [])
    _cc[k] = out; return out

_lab = {}
def labelled(V):
    """toutes les realisations de FHT sur l'ensemble etiquete V (tuple trie)"""
    if V in _lab: return _lab[V]
    if len(V) == 1: _lab[V] = [frozenset()]; return _lab[V]
    res = set()
    for part in set_partitions(list(V)):
        k = len(part)
        if k < 2: continue
        blocks = [tuple(sorted(b)) for b in part]
        subs = [labelled(b) for b in blocks]
        for host in connected_clutters(k):
            top = frozenset(frozenset(x for j in e for x in blocks[j]) for e in host)
            for choice in itertools.product(*subs):
                res.add(top.union(*choice))
    _lab[V] = list(res); return _lab[V]

if __name__ == '__main__':
    attendu = [1, 1, 4, 23, 231]
    t0 = time.time(); ok = True
    print('Methode 1 (force brute sur les axiomes) :')
    for v in range(1, 5):
        n = brute(v); print(f'  v = {v} : {n} FHT a isomorphisme pres (attendu {attendu[v-1]})  {"OK" if n == attendu[v-1] else "ECART"}'); ok &= n == attendu[v-1]
    print('Methode 2 (decomposition recursive) :')
    vmax = 5 if '--v4' not in sys.argv else 4
    for v in range(1, vmax + 1):
        L = labelled(tuple(range(v)))
        n = 1 if v == 1 else len({canon(E, v) for E in L})
        print(f'  v = {v} : {len(L)} FHT etiquetees, {n} a isomorphisme pres (attendu {attendu[v-1]})  {"OK" if n == attendu[v-1] else "ECART"}'); ok &= n == attendu[v-1]
        if 2 <= v <= 4:      # coherence des deux methodes au niveau etiquete : chaque realisation produite passe le test axiomatique
            assert all(is_realisation(E, frozenset(range(v))) for E in L)
    print(f'VERDICT : {"conforme" if ok else "ECART"} ; {time.time()-t0:.1f} s')
