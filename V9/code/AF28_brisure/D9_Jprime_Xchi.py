#!/usr/bin/env python3
"""D9 -- J' exact sur X_chi (test AF9 avec Aut par raffinement, comme AF21 v2).  Criteres : 08_descente_Xgrand/CRITERES_D9_Jprime_Xchi.md.
Usage : python3 D9_Jprime_Xchi.py | tee D9.log   (quelques minutes)"""
import json, time, os, sys
H_ = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H_, '..', 'AF9_test_triplet'))
import AF9_porteurs_triplet as AF9
from AF21b_aut_raffine import aut_sgs_raffine
AF9.aut_sgs = aut_sgs_raffine
for p in ('X_chi.json', os.path.join('..', '..', '06_objets', 'X_chi.json')):
    if os.path.exists(os.path.join(H_, p)): X = json.load(open(os.path.join(H_, p))); break
t = time.time(); o, _ = AF9.test('X_chi', X)
print({k: o.get(k) for k in ['v', 'facets', 'simplices', 'aut', 'ngens', 'betti', 'P1', 'P2', 'P3', 'P4', 'dims', 'dimC1', 'C1', 'C1_shift', 'check_eq']}, f'{time.time()-t:.0f} s')
