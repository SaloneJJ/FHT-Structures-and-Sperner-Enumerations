"""V38 / B1 : N3 complet. Lemme : K_H ne depend que de la partition en H-orbites, et pour H <= S4_L x S4_R x S4_g
cette partition est celle du produit des projections pi_L(H) x pi_R(H) x pi_g(H). Donc les 11^3 = 1331 produits
(a conjugaison pres dans chaque facteur) epuisent les donnees de points fixes de TOUS les sous-groupes de G.
Calcule les Betti Z/2 des K_H pour K26 et ses representants, et compare."""
import sys, json, itertools, collections
sys.path.insert(0, '.')
from fht_tools import *
from chain_inv import SUB4, factors_of, orbits, obj_json
K0, V = build_K3('V25_K3_pieces.json', with_pendants=False)
facs = factors_of(V)
def gens_of(fac, H): return [{v: fac[2](h, v) for v in fac[1]} for h in H]
PARTS = []
for i, HL in enumerate(SUB4):
    for j, HR in enumerate(SUB4):
        for k, HG in enumerate(SUB4):
            PARTS.append(((i, j, k), orbits(V, gens_of(facs[0], HL) + gens_of(facs[1], HR) + gens_of(facs[2], HG))))
def KH(F, O): return maximal([s for s in (frozenset(i for i, o in enumerate(O) if o <= f) for f in F) if s])
objs = collections.OrderedDict([('K26', K0), ('K26*', obj_json('V32_FHT_K26etoile.json')[0]),
        ('K26vol', obj_json('e4_K26_final.json')[0]), ('K26volmin', obj_json(sys.argv[1] if len(sys.argv) > 1 else 'a1_e4_K26_final_s3_best.json')[0])])
res = {}
for nm, K in objs.items():
    res[nm] = [tuple(betti_z2(KH(K, O))) for _, O in PARTS]
    print(nm, 'acycliques', sum(1 for b in res[nm] if not any(b)), 'total Betti', sum(sum(b) for b in res[nm]), flush=True)
ref = res['K26']
for nm in objs:
    d = [PARTS[t][0] for t in range(len(PARTS)) if res[nm][t] != ref[t]]
    print(nm, 'differences avec K26 :', len(d), d[:10])
json.dump({nm: [list(b) for b in r] for nm, r in res.items()}, open('b1_N3_complete.json', 'w'))
