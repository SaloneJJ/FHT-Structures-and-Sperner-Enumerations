import json, sys, time
import fht_tools
D = json.load(open('V145_etoiles_facettes.json'))
for nom in sys.argv[1:]:
    F = [frozenset(f) for f in D['etoiles'][nom]['facets']]; V = sorted(set().union(*F), key=str)
    t = time.time(); o = fht_tools.aut_order(V, F)
    print(nom, '|Aut| (individualisation-raffinement) =', o, '; table :', D['meta'][nom]['aut'], f'({time.time()-t:.0f} s)', flush=True)
