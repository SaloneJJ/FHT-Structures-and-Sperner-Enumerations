"""V65 : origine des classes de K^H (cascade T1) : Betti des pieces A, B, G4, K4 et des parties I, II, I^II sous H. JJ SALONE, FHT 2026"""
import sys, json, itertools
sys.path.insert(0, '.')
from fht_tools import betti_z2, maximal, actside
d = json.load(open('V25_K3_pieces.json'))
A = [frozenset(a) for a in d['A']]; B = [frozenset(b) for b in d['B'] if not any(v.startswith('y') for v in b)]
G4 = d['generation']['feuilles']; E4 = [frozenset(e) for e in d['generation']['K4']]
S4 = list(itertools.permutations(range(4)))
SUB = {'S4': S4, 'S3': [h for h in S4 if h[0] == 0], 'D4': [h for h in S4 if {frozenset({h[0], h[1]}), frozenset({h[2], h[3]})} == {frozenset({0, 1}), frozenset({2, 3})}],
       'S2': [h for h in S4 if h[0] == 0 and h[1] == 1], '1': [tuple(range(4))]}
def orb(verts, act, H):
    O = {}
    for v in verts:
        o = frozenset(act(h, v) for h in H); O[o] = o
    return list(O)
def KH(F, verts, act, H):
    O = orb(verts, act, H)
    return maximal([s for s in (frozenset(o for o in O if o <= f) for f in F) if s]), O
side = lambda h, v: actside(h, v)
gact = lambda h, v: G4[h[G4.index(v)]] if v in G4 else v
VA = sorted(set(v for f in A + B for v in f))
print('A :', betti_z2(A), ' B :', betti_z2(B), ' G4 :', betti_z2([frozenset([g]) for g in G4]), ' K4 :', betti_z2(E4))
for nm, H in SUB.items():
    a, _ = KH(A, VA, side, H); b, _ = KH(B, VA, side, H)
    g, _ = KH([frozenset([x]) for x in G4], G4, gact, H); k, _ = KH(E4, G4, gact, H)
    print(f'H={nm:2s}  A_H {betti_z2(a)}  B_H {betti_z2(b)}  G4_H {betti_z2(g)}  K4_H {betti_z2(k)}', flush=True)
# parties completes
def parts():
    I = [frozenset(['L_' + v for v in a] + ['R_' + v for v in b] + [g]) for a in A for b in B for g in G4]
    II = [frozenset(['L_' + v for v in b] + ['R_' + v for v in a] + list(e)) for b in B for a in A for e in E4]
    return maximal(I), maximal(II)
I, II = parts(); V = sorted(set(v for f in I + II for v in f))
def faces(F):
    s = set()
    for f in F:
        for k in range(1, len(f) + 1): s.update(frozenset(c) for c in itertools.combinations(sorted(f), k))
    return s
IIx = maximal(list(faces(I) & faces(II)))
print('I', betti_z2(I), ' II', betti_z2(II), ' I^II', betti_z2(IIx), flush=True)
def act3(hL, hR, hG):
    def f(v):
        if v.startswith('L_'): return 'L_' + actside(hL, v[2:])
        if v.startswith('R_'): return 'R_' + actside(hR, v[2:])
        return gact(hG, v)
    return f
for (a, b, c) in [('S3', 'S4', 'S4'), ('S4', 'S3', '1'), ('S4', '1', '1'), ('S3', 'S3', 'S4')]:
    H = [act3(x, y, z) for x in SUB[a] for y in SUB[b] for z in SUB[c]]
    out = []
    for nm, F in [('I', I), ('II', II), ('I^II', IIx)]:
        O = {}
        for v in V: o = frozenset(h(v) for h in H); O[o] = o
        O = list(O)
        out.append(f"{nm}_H {betti_z2(maximal([s for s in (frozenset(o for o in O if o <= f) for f in F) if s]))}")
    print(f'H={a}x{b}x{c} : ' + '  '.join(out), flush=True)
