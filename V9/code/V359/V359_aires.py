"""V359 : contrôle des aires de V350 §2 (triangles des règles R1-R3 : rempli = aire 1, sinon aire 3 par parité et cône).
Usage : python3 V359_aires.py <dossier JSON>  (numpy non requis)."""
import json, sys, itertools as it
D = sys.argv[1] if len(sys.argv) > 1 else '.'
for nm, f in (('V37', 'V37_FHT_K26volmin.json'), ('V36', 'V36_FHT_K26vol.json')):
    F = [frozenset(x) for x in json.load(open(f"{D}/{f}"))]
    T = set()
    for x in F:
        T.update(frozenset(c) for c in it.combinations(sorted(x), 3))
    E = set()
    for x in F:
        E.update(frozenset(c) for c in it.combinations(sorted(x), 2))
    V = sorted(set().union(*F))
    def area(c):
        c = frozenset(c)
        if not all(frozenset(e) in E for e in it.combinations(c, 2)): return None
        if c in T: return 1
        cones = [v for v in V if v not in c and all(frozenset((v,)+e) in T for e in it.combinations(sorted(c), 2))]
        return 3 if cones else '>3'
    G = ['g0', 'g1', 'g2', 'g3']
    for side, lab in (('R', 'I'), ('L', 'II')):
        x = [f'{side}_x{i}' for i in range(4)]; z = f'{side}_z'
        r = {'face xxx': {area(c) for c in it.combinations(x, 3)},
             'apex zxx': {area((z,)+c) for c in it.combinations(x, 2)},
             'marche E xxg': {area(c+(g,)) for c in it.combinations(x, 2) for g in G},
             'marche F xgg': {area((xi,)+c) for xi in x for c in it.combinations(G, 2)}}
        print(f"{nm} côté {side} ({lab}) : " + " ; ".join(f"{k} : {sorted(v, key=str)}" for k, v in r.items()))
    print(f"{nm} : générations ggg : {sorted({area(c) for c in it.combinations(G, 3)}, key=str)}")
