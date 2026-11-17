#!/usr/bin/env python3
"""AUDIT 14/11/2026 -- les cinq etoiles exclues de la classe K26 (K_27_559, K_27_1021, K_27_1084, K_28_767, K_28_789) : que passent-elles,
et en quoi different-elles de K26 ?  Pour chacune et pour K^(4), K_27_546 (reference) :
  * groupe d'automorphismes : ordre, groupes induits sur L_x0..3, R_x0..3, G0 = {g0..g3} (S4_L, S4_R, S4_gen intacts ou non) ;
  * tests du registre V133 : N0 (algebre tordue de S4_L = A_F), gamma (motif de Betti), N1 (Betti), J necessaire (sigma Aut sigma = Aut) ;
  * homologie : caracteres des elements qui ne bougent que les sommets L (resp. R) compares a K^(4) (le module « gauche » et « droit » du bimodule) ;
  * difference avec K^(4) : facettes ajoutees / retirees, sommets ajoutes.
numpy ; reutilise AF9, AF3.  Duree : ~1 min.  N'ecrit aucun fichier.   Usage : python3 AF14_cinq_etoiles.py RACINE_FHT26
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import AF9_porteurs_triplet as T
from AF3_K26_homologie_automorphismes import K26_definition

def side_chars(C, H, G, side):
    """caracteres (pair, impair) des elements qui ne bougent que les sommets du cote donne, classes par type de cycle sur side_x0..3"""
    pts = [C.ix[f'{side}_x{i}'] for i in range(4) if f'{side}_x{i}' in C.ix]
    others = [C.ix[v] for v in C.V if not v.startswith(side + '_')]
    out = {}
    for g in G:
        if any(g[o] != o for o in others): continue
        ct = T.ctype(g, pts)
        if ct is None or ct in out: continue
        ME, _ = H.block(0, g); MO, _ = H.block(1, g)
        out[ct] = (T.signed(np.trace(ME)), T.signed(np.trace(MO)))
    return out

def describe(name, facets, ref=None):
    C = T.Cx(name, facets); gens, order = T.aut_sgs(C); nv = len(C.V); G = T.enumerate_group(gens, nv)
    img = lambda labs: len({tuple(g[C.ix[l]] for l in labs if l in C.ix) for g in G})
    Lx = [f'L_x{i}' for i in range(4)]; Rx = [f'R_x{i}' for i in range(4)]; G0 = [f'g{i}' for i in range(4)]
    H = T.Homology(C)
    betti = {d: b for d, b in H.betti.items() if b}
    P1 = all(T.sw(v) in C.ix for v in C.V); sig = [C.ix[T.sw(v)] for v in C.V] if P1 else None
    bad = sum(1 for g in gens if not (P1 and C.is_aut([sig[g[sig[i]]] for i in range(nv)])))
    cL, cR = side_chars(C, H, G, 'L'), side_chars(C, H, G, 'R')
    d = dict(name=name, v=nv, aretes=len(C.F), simplexes=C.nsimp, aut=order, SL=img(Lx), SR=img(Rx), Sgen=img(G0), betti=betti,
             J_nec_defauts=f'{bad}/{len(gens)}', cL=cL, cR=cR, extra=sorted(set(C.V) - set(ref['V'])) if ref else [])
    if ref:
        F = {frozenset(f) for f in facets}; Fr = ref['F']
        d['ajoutees'] = sorted([sorted(f) for f in F - Fr]); d['retirees'] = sorted([sorted(f) for f in Fr - F])
        d['memes_L'] = cL == ref['cL']; d['memes_R'] = cR == ref['cR']
    return d, C

if __name__ == '__main__':
    root = sys.argv[1]
    et = json.load(open(os.path.join(root, '03_code', 'catalogues', 'V145_etoiles_facettes.json')))['etoiles']
    K4 = K26_definition()
    d0, C0 = describe('K^(4)', K4)
    ref = dict(V=C0.V, F={frozenset(f) for f in K4}, cL=d0['cL'], cR=d0['cR'])
    print(f"K^(4) : v={d0['v']} aretes={d0['aretes']} |Aut|={d0['aut']} induits L/R/gen = {d0['SL']}/{d0['SR']}/{d0['Sgen']} Betti={d0['betti']}")
    print(f"   caracteres L (type -> pair, impair) : {d0['cL']}\n   caracteres R : {d0['cR']}")
    for nm in ('K_27_546', 'K_27_559', 'K_27_1021', 'K_27_1084', 'K_28_767', 'K_28_789'):
        d, _ = describe(nm, et[nm]['facets'], ref)
        print(f"\n{nm} : v={d['v']} aretes={d['aretes']} simplexes={d['simplexes']} |Aut|={d['aut']} ; groupes induits L/R/gen = {d['SL']}/{d['SR']}/{d['Sgen']} ; Betti={d['betti']}")
        print(f"   J necessaire (sigma Aut sigma = Aut) : generateurs en defaut {d['J_nec_defauts']} ; module gauche = K^(4) : {d['memes_L']} ; module droit = K^(4) : {d['memes_R']}")
        if not d['memes_R']: print(f"   caracteres R : {d['cR']}")
        if not d['memes_L']: print(f"   caracteres L : {d['cL']}")
        print(f"   sommets ajoutes : {d['extra']} ; aretes ajoutees : {d['ajoutees']} ; retirees : {d['retirees']}")
