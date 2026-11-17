#!/usr/bin/env python3
"""D16 -- assemblage de la piece compressee A' (motif A0, |Aut| = 24) avec la piece couleur B de K26 ; test AF9 complet avec C2 (caracteres de K26).
Criteres : CRITERES_D14_D16.md (D16').  Usage : python3 D16_assemblage.py | tee D16.log"""
import sys, os, json, itertools, time
H_ = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(H_, 'D16p_compression.py')).read().split("seen = set(); orbs = []")[0]; exec(src)
# piece retenue : isoler {l0, t6} puis cone sur le sous-complexe induit sur les orbites de l0, l2, t0, v2, v5
orbit = lambda v: sorted({g.get(v, v) for g in VPd})
Iset = set(orbit('l0')) | set(orbit(ta))
Fi = [f for f in C0 if not (set(f) & Iset)] + [[v] for v in sorted(Iset)]
W = set().union(*[set(orbit(v)) for v in ('l0', 'l2', 't0', 'v2', 'v5')])
Y = [list(s) for s in maximal([frozenset(s) for s in faces_of(Fi) if set(s) <= W])]
Ap = Fi + [['N'] + y for y in Y]
VA = sorted({x for f in Ap for x in f})
print(f"A' : sommets {len(VA)}, Betti {betti_z2(Ap)}, |Aut| = {aut_order(VA, Ap)}, homologie sous S4 diagonal {sorted(eh(Ap).items())}", flush=True)
json.dump(Ap, open(os.path.join(H_, 'A_compresse.json'), 'w'))
d = json.load(open(os.path.join(H_, 'V25_K3_pieces.json'))); Bk = d['B']; A0 = d['A']
print('controle A0 (K26) : homologie sous S4 (sur les x) :', 'Betti', betti_z2(A0))
G4, E4 = d['generation']['feuilles'], d['generation']['K4']
def build(A, B):
    I = [frozenset(['L_' + v for v in a] + ['R_' + v for v in b] + [g]) for a in A for b in B for g in G4]
    II = [frozenset(['L_' + v for v in b] + ['R_' + v for v in a] + e) for b in B for a in A for e in E4]
    K = maximal(I + II); return [sorted(f) for f in K]
Xc = build(Ap, Bk); Vc = sorted({x for f in Xc for x in f}); t = time.time()
print(f"X_comp : sommets {len(Vc)}, facettes {len(Xc)}, Betti {betti_z2(Xc)}, |Aut| = {aut_order(Vc, Xc)} (prevu 13824) ({time.time()-t:.0f} s)", flush=True)
json.dump(Xc, open(os.path.join(H_, 'X_compresse.json'), 'w'))
sys.path.insert(0, os.path.join(H_, '..', 'AF9_test_triplet'))
import AF9_porteurs_triplet as AF9
from AF21b_aut_raffine import aut_sgs_raffine
AF9.aut_sgs = aut_sgs_raffine
K26 = build(A0, Bk)
t = time.time(); oref, ref = AF9.test('K26', K26, want_chars=True)
print('K26 (reference) :', {k: oref.get(k) for k in ['v', 'aut', 'betti', 'P1', 'P2', 'P3', 'P4', 'dimC1', 'C1', 'C1_shift']}, f'({time.time()-t:.0f} s)', flush=True)
t = time.time(); o, _ = AF9.test('X_comp', Xc, ref_chars=ref)
print('X_comp :', {k: o.get(k) for k in ['v', 'aut', 'betti', 'P1', 'P2', 'P3', 'P4', 'dimC1', 'C1', 'C1_shift', 'C2', 'C2_labels', 'C2_incoherent']}, f'({time.time()-t:.0f} s)')
print('verdict AF9 :', AF9.verdict(o))
