#!/usr/bin/env python3
"""GS1, G-f : |Aut(X_GS)|, X_GS = X_grand * D8. Attendu |Aut(X_grand)| x 6 = 4 777 574 400. Usage : python3 gs1_aut.py | tee GS1_aut.log"""
import sys, os, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..', 'FHT_28', '06_code', 'AF20_algebres_tordues')))
from fht_tools import aut_order
X = json.load(open(os.path.join(HERE, '..', '..', '..', 'FHT_28', '06_code', 'AF20_algebres_tordues', 'X_grand.json')))
X = X['facettes'] if isinstance(X, dict) and 'facettes' in X else (X['facets'] if isinstance(X, dict) else X)
T8 = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)]
D = [[f'D{i}_{a}', f'D{i}_{b}'] for i in range(3) for a, b in T8]
J = [list(f) + d for f in X for d in D]
V = sorted({x for f in J for x in f}, key=str)
print('X_grand : facettes', len(X), '; X_GS : sommets', len(V), 'facettes', len(J), flush=True)
t = time.time(); print('|Aut(X_GS)| =', aut_order(V, J), f'({time.time()-t:.0f} s)', flush=True)
