"""Test complet AF9 (sigma, Aut, P1-P4, J' exact C1) sur X_grand (56 sommets, assemblage de la grande symetrie).
Duree estimee > 10 min (homologie sur F_p avec action d'un groupe d'ordre 796 262 400 par generateurs). Usage : python3 AF21_grand_Jprime.py"""
import json, time, os
from AF9_porteurs_triplet import test
X = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '06_objets', 'X_grand.json')))
t = time.time(); o, _ = test('X_grand', X)
print({k: o.get(k) for k in ['v', 'facets', 'simplices', 'aut', 'betti', 'P1', 'P2', 'P3', 'P4', 'dims', 'dimC1', 'C1', 'C1_shift', 'check_eq']}, f'{time.time()-t:.0f} s')
