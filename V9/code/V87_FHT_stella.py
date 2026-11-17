"""V87 : cube 'stella octangula' (directions = 3 appariements, L/R = les deux tetraedres). JJ SALONE, FHT 2026"""
import itertools
from V86_FHT_cube import analyse, rules_cube, glue, FK, xn
from V86_FHT_k26 import betti_from_faces
def H(k): return {0, k}
def rules_stella():
    R = {}
    for j, k in enumerate((1, 2, 3)):
        plus = [('L', c) for c in range(4) if c in H(k)] + [('R', c) for c in range(4) if c not in H(k)]
        R[j] = [(xn(S, c), xn('R' if S == 'L' else 'L', c ^ k)) for S, c in plus]
    return R
def faces_of(rules_dirs):
    fam = []
    for d, pairs in rules_dirs.items():
        fam.append(frozenset(v for v, w in pairs)); fam.append(frozenset(w for v, w in pairs))
    return fam
# stabilisateurs dans S4^3
PAIR = [((0,1),(2,3)), ((0,2),(1,3)), ((0,3),(1,2))]
S4 = list(itertools.permutations(range(4)))
def s3img(g):
    out = []
    for P in PAIR:
        img = frozenset(frozenset(g[x] for x in b) for b in P)
        out.append([frozenset(frozenset(b) for b in Q) for Q in PAIR].index(img))
    return out
def act(gL, gR, v):
    S, name = v.split('_')
    g = gL if S == 'L' else gR
    if name[0] == 'x': return f'{S}_x{g[int(name[1])]}'
    return v
def stab(fam):
    F = set(fam); cnt = 0
    for gL in S4:
        for gR in S4:
            if set(frozenset(act(gL, gR, v) for v in f) for f in F) == F: cnt += 1
    return cnt
if __name__ == '__main__':
    RA = rules_cube(); RA0 = dict(RA); RA0[0] = [p for p in RA[0] if 'z' not in p[0]]
    RB = rules_stella()
    print('faces (B) :', [sorted(f) for f in faces_of(RB)])
    print('stabilisateur dans S4_L x S4_R (x S4_gen, facteur 24 non compte) : (A) =', stab(faces_of(RA0)), ' (B) =', stab(faces_of(RB)), flush=True)
    analyse([(0,0,0),(1,0,0)], RB, '(B) chaine, 2 tuiles')
    analyse([(0,0,0),(1,0,0),(0,1,0),(1,1,0)], RB, '(B) plaquette')
    box = list(itertools.product(range(2), repeat=3))
    analyse(box, RB, '(B) bloc 2x2x2')
    f = glue(list(itertools.product(range(3), repeat=3)), RB)
    fs = set()
    for t in itertools.product(range(3), repeat=3):
        for s in FK: fs.add(tuple(sorted(str(f((t, v))) for v in s)))
    print('(B) bloc 3x3x3 :', betti_from_faces(fs), 'attendu 324, 972, 972, 324', flush=True)
