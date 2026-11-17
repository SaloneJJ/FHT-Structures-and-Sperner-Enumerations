"""Enumeration des FHT (A5 par composant) sur v <= 5 sommets, par decomposition recursive
(V0_FHT_DefinitionAndEnumeration, def. recursive decomposition) ; classes d'isomorphisme
par forme canonique brute (permutations). Controle : 1, 1, 4, 23, 231."""
import itertools
from functools import lru_cache

def connected_clutters(n):
    """clutters connexes sur {0..n-1}, aretes de taille >= 2 si n >= 2, couvrant"""
    if n == 1: return [((0,),)]
    subs = [frozenset(s) for k in range(2, n + 1) for s in itertools.combinations(range(n), k)]
    out = []
    def conn(E):
        seen = set(E[0]); ch = True
        while ch:
            ch = False
            for e in E:
                if e & seen and not e <= seen: seen |= e; ch = True
        return seen == set(range(n))
    def bt(i, chosen):
        if chosen and set().union(*chosen) == set(range(n)) and conn(chosen): out.append(tuple(chosen))
        for j in range(i, len(subs)):
            s = subs[j]
            if any(s <= c or c <= s for c in chosen): continue
            bt(j + 1, chosen + [s])
    bt(0, []); return out

def set_partitions(lst):
    if not lst: yield []; return
    first, rest = lst[0], lst[1:]
    for p in set_partitions(rest):
        for i in range(len(p)): yield p[:i] + [[first] + p[i]] + p[i + 1:]
        yield [[first]] + p

CC = {}
def cc(n):
    if n not in CC: CC[n] = connected_clutters(n)
    return CC[n]

@lru_cache(None)
def fhts(V):
    """realisations (frozenset de frozensets) de tous les FHT etiquetes sur le tuple V"""
    V = tuple(V); n = len(V)
    if n == 1: return [frozenset()]
    out = set()
    for P in set_partitions(list(V)):
        if len(P) < 2: continue
        P = [tuple(sorted(b)) for b in P]
        subs = [fhts(b) for b in P]
        for host in cc(len(P)):
            top = [frozenset(x for i in e for x in P[i]) for e in host]
            for choice in itertools.product(*subs):
                out.add(frozenset(top).union(*choice))
    return list(out)

def canon(n, E):
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(sorted(tuple(sorted(p[x] for x in e)) for e in E))
        if best is None or key < best: best = key
    return best

if __name__ == '__main__':
    for n in range(1, 6):
        L = fhts(tuple(range(n)))
        classes = {canon(n, E) for E in L}
        print(n, len(classes))
