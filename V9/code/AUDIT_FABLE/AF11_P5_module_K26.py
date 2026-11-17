#!/usr/bin/env python3
"""PROPOSITION P5 (05_audit/AUDIT_FABLE_PROPOSITIONS.md, criteres ecrits le 13/11/2026 a 20 h 30) : structure de module de H~(K26) sous
S4_L x S4_R x S4_gen, degre par degre, par un code independant ; comparaison exacte avec l'enonce du theoreme « Homology » de l'article 2 :
  degre 3 : (1,1,3) + (1,3,3) ; degre 4 : (2,1,3) + (2,3,3) + (1,1,3') + (3,1,3') ; degre 5 : (s,1,3) + (s,3,3) + (1,2,3') + (3,2,3') ;
  degre 6 : (1,s,3') + (3,s,3')      [(L, R, gen) ; 1, s = signe, 2, 3 = standard, 3' = 3 x s]
Methode : homologie sur F_65521 avec representants (AF9), caractere de chaque degre sur un element par classe de conjugaison de S4^3
(classe = types de cycle sur L_x0..3, R_x0..3, g0..3), multiplicites par orthogonalite des caracteres.  Critere : egalite exacte.
Usage : python3 AF11_P5_module_K26.py RACINE_FHT26          (numpy ; n'ecrit aucun fichier)
"""
import os, sys, json, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import AF9_porteurs_triplet as T
from AF3_K26_homologie_automorphismes import K26_definition

CT = [(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)]
SIZE = {(1, 1, 1, 1): 1, (2, 1, 1): 6, (2, 2): 3, (3, 1): 8, (4,): 6}
CHI = {'1': [1, 1, 1, 1, 1], 's': [1, -1, 1, 1, -1], '2': [2, 0, 2, -1, 0], '3': [3, 1, -1, 0, -1], "3'": [3, -1, -1, 0, 1]}
CHI = {k: dict(zip(CT, v)) for k, v in CHI.items()}
EXPECT = {3: {('1', '1', '3'), ('1', '3', '3')},
          4: {('2', '1', '3'), ('2', '3', '3'), ('1', '1', "3'"), ('3', '1', "3'")},
          5: {('s', '1', '3'), ('s', '3', '3'), ('1', '2', "3'"), ('3', '2', "3'")},
          6: {('1', 's', "3'"), ('3', 's', "3'")}}

def run(name, facets):
    C = T.Cx(name, facets); gens, order = T.aut_sgs(C); nv = len(C.V)
    assert order == 13824, order
    G = T.enumerate_group(gens, nv); H = T.Homology(C)
    Lx = [C.ix[f'L_x{i}'] for i in range(4)]; Rx = [C.ix[f'R_x{i}'] for i in range(4)]; Gx = [C.ix[f'g{i}'] for i in range(4)]
    reps = {}
    for g in G:
        lab = (T.ctype(g, Lx), T.ctype(g, Rx), T.ctype(g, Gx))
        if None in lab: raise SystemExit('un automorphisme ne preserve pas L_x, R_x ou G0')
        reps.setdefault(lab, g)
    assert len(reps) == 125, len(reps)
    degs = [d for d, b in sorted(H.betti.items()) if b]
    chars = {d: {lab: T.signed(np.trace(H.mat(d, g))) for lab, g in reps.items()} for d in degs}
    ok = True; out = {}
    for d in degs:
        mult = {}
        for a, b, c in itertools.product(CHI, repeat=3):
            s = sum(SIZE[l] * SIZE[r] * SIZE[q] * chars[d][(l, r, q)] * CHI[a][l] * CHI[b][r] * CHI[c][q] for (l, r, q) in reps)
            assert s % 13824 == 0, (d, a, b, c, s)
            if s: mult[(a, b, c)] = s // 13824
        found = {k for k, v in mult.items()}
        good = (found == EXPECT.get(d, set())) and all(v == 1 for v in mult.values())
        ok &= good; out[d] = mult
        print(f'  {name} degre {d} (Betti {H.betti[d]}) : ' + ' + '.join(f"{v}x({a},{b},{c})" for (a, b, c), v in sorted(mult.items()))
              + ('   == enonce' if good else f'   ECART (attendu {sorted(EXPECT.get(d, set()))})'), flush=True)
    print(f'  {name} : ' + ('CONFORME au theoreme « Homology »' if ok and set(degs) == {3, 4, 5, 6} else 'ECART'), flush=True)
    return ok

if __name__ == '__main__':
    root = sys.argv[1]
    print('=== P5 : structure de module de H~(K26) sous S4_L x S4_R x S4_gen (F_65521)')
    ok = run('K^(4)', K26_definition())
    d = json.load(open(os.path.join(root, '04_objets', 'V37_FHT_K26volmin.json'))); fac = d if isinstance(d, list) else d.get('facets', d.get('facettes'))
    ok &= run('K26oo (V37)', fac)
    print('VERDICT P5 :', 'conforme (critere : egalite exacte des caracteres en chaque degre 3-6)' if ok else 'ECART : erreur dans T4')
