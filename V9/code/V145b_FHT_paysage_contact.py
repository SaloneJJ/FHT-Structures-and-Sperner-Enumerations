"""V145b : tri du paysage (V144 §6, criteres L1-L5) pour les etoiles vertes de V52_paysage_v5_f30.
L1 structure de contact stella (8 x, faces (B) de V87 simplexes de la realisation, Q1-Q4 de V85) ;
L2 stabilisateur de la famille des faces dans Aut(F), action sur les directions et les coins ;
L3 4-cubes de contact (16 coins, 8 faces simplexes de 8 sommets), recherche exhaustive si v <= 28 ;
L4 dimension de H_F (Betti) ; L5 aretes emboitees (profondeur 2) : quels sommets touchent-elles ?
Donnees : V145_etoiles_facettes.json (extrait en lecture seule de V133_FHT_facettes.jsonl et de la table, 11/11 13 h 10).
JJ SALONE, FHT 2026"""
import json, itertools, time, sys
from collections import Counter, deque
from V86_FHT_k26 import build_K
from V87_FHT_stella import rules_stella, faces_of
import fht_tools

D = json.load(open('V145_etoiles_facettes.json'))
ETO = D['etoiles']; META = D['meta']
FB = [frozenset(f) for f in faces_of(rules_stella())]          # 6 faces (B), par paires (F+, F-)
XS = frozenset().union(*FB)
TMAX4 = float(sys.argv[1]) if len(sys.argv) > 1 else 900.0      # budget (s) de la recherche de 4-cubes par objet

def groupe(V, gens):
    """enumeration du groupe engendre (permutations sur V, en tuples d'indices)"""
    ix = {v: i for i, v in enumerate(V)}
    G0 = [tuple(ix[g[v]] for v in V) for g in gens]
    e = tuple(range(len(V))); seen = {e}; q = deque([e])
    while q:
        p = q.popleft()
        for g in G0:
            r = tuple(g[i] for i in p)
            if r not in seen: seen.add(r); q.append(r)
    return seen, ix

def analyse(nom, F, gens, v_att=None):
    t0 = time.time()
    F = [frozenset(f) for f in F]
    V = sorted(set().union(*F), key=str); n = len(V)
    maxi = fht_tools.maximal(F)
    emb = [f for f in F if any(f < g for g in F)]
    def simplexe(s): return any(s <= f for f in maxi)
    res = {'v': n, 'facettes': len(F), 'maximales': len(maxi), 'emboitees': len(emb)}
    # L1
    res['x_presents'] = XS <= set(V)
    res['faces_B_simplexes'] = [simplexe(f) for f in FB] if res['x_presents'] else None
    res['Q1_opposees_disjointes'] = all(not (FB[2*k] & FB[2*k+1]) for k in range(3))
    res['tout_x_dans_un_simplexe'] = simplexe(XS) if res['x_presents'] else None
    # L2
    if gens:
        G, ix = groupe(V, gens)
        res['Aut'] = len(G)
        FBi = [frozenset(ix[v] for v in f) for f in FB] if res['x_presents'] else []
        fam = set(FBi); stab = []
        for p in G:
            if fam and set(frozenset(p[i] for i in f) for f in FBi) == fam: stab.append(p)
        res['stab_faces_B'] = len(stab)
        if stab:
            dirs = [frozenset((FBi[2*k], FBi[2*k+1])) for k in range(3)]
            orb_dir = set()
            for p in stab:
                img = frozenset(frozenset(frozenset(p[i] for i in f) for f in d) for d in [dirs[0]])
                orb_dir.add(next(k for k in range(3) if img == frozenset([dirs[k]])))
            res['stab_transitif_directions'] = (len(orb_dir) == 3)
            xi = sorted(ix[x] for x in XS)
            orb = {}
            for x in xi:
                orb[x] = frozenset(p[x] for p in stab)
            res['stab_orbites_coins'] = sorted(Counter(len(o) for o in set(orb.values())).items())
        # Aut preserve-t-il l'ensemble des x ?
        if res['x_presents']:
            Xi = frozenset(ix[x] for x in XS)
            res['Aut_preserve_x'] = all(frozenset(p[i] for i in Xi) == Xi for p in G)
    # L5 : aretes emboitees
    res['emboitees_detail'] = [sorted(f, key=str) for f in emb][:6]
    res['emboitees_touchent_x'] = any(f & XS for f in emb)
    # L3 : 4-cubes
    idx = {v: i for i, v in enumerate(V)}
    FM = [sum(1 << idx[v] for v in f) for f in maxi]
    def is_face(m): return any(m & Fm == m for Fm in FM)
    nonadj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and not is_face((1 << i) | (1 << j)): nonadj[i] |= 1 << j
    res['paires_non_adjacentes'] = sum(bin(x).count('1') for x in nonadj) // 2
    # recherche par retour arriere : 16 sommets, chacun non adjacent a au plus un autre du lot
    cand = []; t1 = time.time(); fini = True
    order = sorted(range(n), key=lambda i: -bin(nonadj[i]).count('1'))
    def bt(start, chosen_mask, chosen, cnt):
        nonlocal fini
        if time.time() - t1 > TMAX4: fini = False; return
        if len(chosen) == 16: cand.append(chosen_mask); return
        for pos in range(start, n):
            if n - pos < 16 - len(chosen): return
            i = order[pos]
            c = bin(nonadj[i] & chosen_mask).count('1')
            if c > 1: continue
            # chaque deja choisi touche par i ne doit pas depasser 1
            if any(cnt[j] + 1 > 1 for j in chosen if nonadj[i] >> j & 1): continue
            for j in chosen:
                if nonadj[i] >> j & 1: cnt[j] += 1
            cnt[i] = c
            bt(pos + 1, chosen_mask | (1 << i), chosen + [i], cnt)
            for j in chosen:
                if nonadj[i] >> j & 1: cnt[j] -= 1
            cnt[i] = 0
            if not fini: return
    bt(0, 0, [], [0] * n)
    cubes = 0
    for m in cand:
        bits = [i for i in range(n) if m >> i & 1]
        # 8-faces contenues dans m : sous-ensembles des intersections (facette & m) de taille >= 8
        inter = {Fm & m for Fm in FM}
        inter = [x for x in inter if bin(x).count('1') >= 8]
        f8 = set()
        for x in inter:
            xb = [i for i in range(n) if x >> i & 1]
            for S in itertools.combinations(xb, 8): f8.add(sum(1 << i for i in S))
        halves = list({min(s, m ^ s) for s in f8 if (m ^ s) in f8})
        for quad in itertools.combinations(halves, 4):
            if all(bin(a & c).count('1') == 4 for a, c in itertools.combinations(quad, 2)):
                sig = set()
                ok = True
                for i in bits:
                    key = tuple((q >> i) & 1 for q in quad)
                    if key in sig: ok = False; break
                    sig.add(key)
                if ok: cubes += 1
    res['4cubes'] = {'16-ensembles_compatibles': len(cand), 'recherche_complete': fini, '4-cubes': cubes}
    res['duree_s'] = round(time.time() - t0, 1)
    return res

K, I, II = build_K()
objets = [('K26', [sorted(f) for f in K], None)]
for nom in ['K_27_546', 'K_27_559', 'K_27_1021', 'K_27_1084', 'K_28_767', 'K_28_789', 'K_36_44', 'K_34_1', 'K_36_39']:
    objets.append((nom, ETO[nom]['facets'], ETO[nom]['gens']))
# generateurs de Aut(K26) : S4_L x S4_R x S4_gen (actions de fht_tools), pour le stabilisateur
def gens_k26(V):
    out = []
    S4g = [(1, 0, 2, 3), (1, 2, 3, 0)]
    idn = (0, 1, 2, 3)
    for gL, gR, gg in [(g, idn, idn) for g in S4g] + [(idn, g, idn) for g in S4g] + [(idn, idn, g) for g in S4g]:
        out.append({v: fht_tools.act3(gL, gR, gg, v) for v in V})
    return out
for nom, F, gens in objets:
    Vn = sorted(set().union(*map(set, F)), key=str)
    if nom == 'K26':
        gens = gens_k26(Vn)
        # verification : ce sont des automorphismes
        Fs = set(frozenset(f) for f in F)
        assert all(set(frozenset(g[v] for v in f) for f in Fs) == Fs for g in gens), 'generateurs K26 invalides'
    r = analyse(nom, F, gens)
    m = META.get(nom, {})
    r['H_F'] = sum((m.get('betti') or {}).values()) if m else None
    r['profondeur'] = m.get('depth'); r['mouvement'] = m.get('move'); r['parent'] = m.get('parent')
    print(nom, json.dumps(r, ensure_ascii=False), flush=True)
