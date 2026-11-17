"""AF21 v2 (14/11/2026 soir) : meme test que AF21, mais Aut calcule par raffinement de couleurs (AF21b_aut_raffine).
Diagnostic : AF21 ne terminait pas a cause de aut_sgs (AF9), retour arriere sans raffinement, explosif pour |Aut| = 796 262 400.
Duree : ~2 min.  Resultat : voir AF21_grand_Jprime_v2.log.  Usage : python3 AF21_grand_Jprime_v2.py"""
import json, time, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import AF9_porteurs_triplet as AF9
from AF21b_aut_raffine import aut_sgs_raffine
AF9.aut_sgs = aut_sgs_raffine
for p in ('AF20_algebres_tordues/X_grand.json', '../../06_objets/X_grand.json', 'X_grand.json'):
    if os.path.exists(os.path.join(HERE, p)): X = json.load(open(os.path.join(HERE, p))); break
t = time.time(); o, _ = AF9.test('X_grand', X)
print({k: o.get(k) for k in ['v', 'facets', 'simplices', 'aut', 'ngens', 'betti', 'P1', 'P2', 'P3', 'P4', 'dims', 'dimC1', 'C1', 'C1_shift', 'check_eq']}, f'{time.time()-t:.0f} s')
print('NB : le verdict AF9 "BIMODULE DIFFERENT" vient de ce que C2 (comparaison aux caracteres de K26) n est pas demande ici.')
