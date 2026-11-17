"""H1.1 : K26 sous la hierarchie arborescente, lecture R1.
(a) K26 admet-il une decomposition R1 non triviale ? test du critere de module sur des candidats structures ;
(b) version arborescente 'coeur duplique' : enfants A'^L, q^L, B^L, A'^R, q^R, B^R (coeurs disjoints), feuilles g_i."""
import sys, json, itertools, time
sys.path.insert(0, '.')
from fht_tools import *
from fht_hier import *
d = json.load(open('V25_K3_pieces.json'))
K, V = build_K3('V25_K3_pieces.json', with_pendants=False)
Kset = set(K)
def is_module(C, K):
    C = frozenset(C)
    ctx = {}
    for s in K:
        t = s & C
        if not t: continue
        ctx.setdefault(t, set()).add(s - C)
    T = list(ctx)
    if any(a < b for a in T for b in T): return False
    cs = list(ctx.values())
    return all(c == cs[0] for c in cs)
# (a) candidats : unions d'orbites de feuilles ; toutes les paires et triples
grp = group_perms_K3(V)
orbs = []
seen = set()
for x in V:
    if x in seen: continue
    o = frozenset(g[x] for g in grp); seen |= o; orbs.append(o)
print("orbites de feuilles :", len(orbs), [len(o) for o in orbs])
t0 = time.time(); found = []
for r in range(1, len(orbs)):
    for combo in itertools.combinations(orbs, r):
        C = frozenset().union(*combo)
        if 2 <= len(C) < len(V) and is_module(C, K): found.append(sorted(C))
print("modules parmi les unions d'orbites :", len(found), found[:5], f"({time.time()-t0:.0f} s)")
small = [c for k in (2, 3) for c in itertools.combinations(V, k) if is_module(c, K)]
print("modules a 2 ou 3 feuilles :", len(small), small[:5])
# (b) version arborescente, coeurs dupliques
A_ = [frozenset(f) for f in maximal(d['A'])]
Ap = [f for f in A_ if f != frozenset(['q'])]
B_ = [frozenset(f) for f in maximal([b for b in d['B'] if not any(v.startswith('y') for v in b)])]
def tag(F, s): return [frozenset(s + v for v in f) for f in F]
AL, AR = tag(Ap, 'LA_'), tag(Ap, 'RA_'); BL, BR = tag(B_, 'LB_'), tag(B_, 'RB_')
qL, qR = frozenset(['L_q']), frozenset(['R_q'])
G4 = ['g0', 'g1', 'g2', 'g3']
lv = lambda F: frozenset().union(*F)
top = [lv(AL) | lv(BR) | {g} for g in G4] + [qL | lv(BR) | {g} for g in G4]
top += [lv(BL) | lv(AR) | {a, b} for a, b in itertools.combinations(G4, 2)] + [lv(BL) | qR | {a, b} for a, b in itertools.combinations(G4, 2)]
E = AL + AR + BL + BR + top
Vt = sorted(lv(E))
t = decompose(set(Vt), E)
print("arbre : profondeur", depth(t), "; enfants de la racine :", len(t[3]), "; feuilles :", len(Vt))
R = R1(t)
print("R1 : facettes", len(R), "; Betti Z/2 :", betti_z2(R))
print("|Aut(F)| (realisation) :", aut_order(Vt, E), " ; |Aut(ombre plate R1)| :", aut_order(Vt, R))
kown = len(E) - 0
print("kappa propre (aretes de tous les composants) :", len(E), " ; kappa plat (facettes R1) :", len(R))
print("K26 : feuilles 26, facettes 716, Betti", betti_z2(K))
