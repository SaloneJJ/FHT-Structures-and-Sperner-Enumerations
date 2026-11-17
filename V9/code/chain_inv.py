"""V34 / E1 : invariants N2-N6 de la chaine de classes pour les objets du programme.
Objets : K26, K3, K26*, K^(2), K^(2)*, K^(6), K^(6)*, K25, K33, K35.
N3 : Betti Z/2 de K_H (complexe des H-orbites, lemme de V34) pour 264 sous-groupes H."""
import sys, json, itertools, math, collections
sys.path.insert(0, '.')
from fht_tools import *
from g_family import Kn

# ---------------- sous-groupes de S4 a conjugaison pres (11 classes)
def mul(g, h): return tuple(g[h[i]] for i in range(len(g)))
def inv(g):
    r = [0] * len(g)
    for i, x in enumerate(g): r[x] = i
    return tuple(r)
def closure(gens, n):
    e = tuple(range(n)); S = {e}; fr = [e]
    while fr:
        nf = []
        for a in fr:
            for g in gens:
                b = mul(g, a)
                if b not in S: S.add(b); nf.append(b)
        fr = nf
    return frozenset(S)
def subgroup_classes(n):
    G = list(itertools.permutations(range(n)))
    subs = set()
    for a in G:
        for b in G: subs.add(closure([a, b], n))
    cls = []
    for H in sorted(subs, key=len):
        conj = {frozenset(mul(mul(g, h), inv(g)) for h in H) for g in G}
        if not any(c in conj for c in cls): cls.append(H)
    return cls
SUB4 = subgroup_classes(4)
assert len(SUB4) == 11
NAMES4 = {1: '1', 2: 'C2', 3: 'C3', 4: 'V4|C4', 6: 'S3', 8: 'D4', 12: 'A4', 24: 'S4'}

# ---------------- objets : (facettes, feuilles, facteurs) ; facteur = (nom, feuilles, action(g, v), groupe complet)
def side_factor(name, pre, V):
    vs = [v for v in V if v.startswith(pre)]
    return (name, vs, lambda g, v, pre=pre: pre + actside(g, v[len(pre):]), list(itertools.permutations(range(4))))
def gen_factor(V):
    gs = sorted(v for v in V if v.startswith('g'))
    n = len(gs)
    return ('gen', gs, lambda g, v: 'g%d' % g[int(v[1:])], list(itertools.permutations(range(n))))
def obj_K3(pend):
    K, V = build_K3('V25_K3_pieces.json', with_pendants=pend); return K, V
def obj_json(fn):
    d = json.load(open(fn)); F = d['facettes'] if isinstance(d, dict) else d
    K = [frozenset(f) for f in F]; return K, sorted(set(v for f in K for v in f))
def obj_gadget(name):
    if name == 'K35':
        d = json.load(open('V22_FHT_K35_biS4.json')); V = d['feuilles']
        K = [frozenset(V[i] for i in e) for e in d['aretes']]
    else:
        d = json.load(open('V23_K33_pieces.json')); A_, B_ = d['A'], d['B']
        if name == 'K25': B_ = [b for b in B_ if not any(v.startswith('y') for v in b)]
        K = [frozenset(['L_' + v for v in a] + [s] + ['R_' + v for v in b]) for a in A_ for b in B_ for s in ('s+', 's-')]
        K += [frozenset(['L_' + v for v in b] + ['R_' + v for v in a]) for b in B_ for a in A_] + [frozenset(['s+', 't'])]
        K = maximal(K)
    return K, sorted(set(v for f in K for v in f))
OBJ = collections.OrderedDict()
OBJ['K26'] = obj_K3(False); OBJ['K3'] = obj_K3(True); OBJ['K26*'] = obj_json('V32_FHT_K26etoile.json')
OBJ['K(2)'] = Kn(2); OBJ['K(2)*'] = obj_json('K2_etoile.json'); OBJ['K(6)'] = Kn(6); OBJ['K(6)*'] = obj_json('K6_etoile.json')
for nm in ('K25', 'K33', 'K35'): OBJ[nm] = obj_gadget(nm)

def factors_of(V):
    f = [side_factor('L', 'L_', V), side_factor('R', 'R_', V)]
    if any(v.startswith('g') for v in V): f.append(gen_factor(V))
    return f

def orbits(V, gens):
    """orbites des feuilles sous le groupe engendre par une liste de permutations (dicts partiels)"""
    par = {v: v for v in V}
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for g in gens:
        for v, w in g.items(): par[find(v)] = find(w)
    out = collections.defaultdict(set)
    for v in V: out[find(v)].add(v)
    return [frozenset(o) for o in out.values()]

def fixed_complex(K, V, gens):
    O = orbits(V, gens)
    idx = {}
    for i, o in enumerate(O):
        for v in o: idx[v] = i
    F = []
    for f in K:
        s = frozenset(i for i, o in enumerate(O) if o <= f)
        if s: F.append(s)
    return maximal(F) if F else []

def subgroup_list(facs):
    """(etiquette, generateurs) pour H_L x H_R x H_g et diagonales"""
    L, R = facs[0], facs[1]; Gf = facs[2] if len(facs) > 2 else None
    def gens_of(fac, H): return [{v: fac[2](h, v) for v in fac[1]} for h in H]
    gH = [('1', [])] if Gf is None else [('1', []), ('Sgen', gens_of(Gf, Gf[3]))]
    out = []
    for i, HL in enumerate(SUB4):
        for j, HR in enumerate(SUB4):
            for lg, gg in gH:
                out.append((f"L{i}xR{j}x{lg}", gens_of(L, HL) + gens_of(R, HR) + gg))
    for i, H in enumerate(SUB4):
        for lg, gg in gH:
            dg = []
            for h in H:
                m = {v: L[2](h, v) for v in L[1]}; m.update({v: R[2](h, v) for v in R[1]}); dg.append(m)
            out.append((f"D{i}x{lg}", dg + gg))
    return out

def N4(V, facs):
    rows = []
    for name, vs, act, Gfull in facs:
        seen = set()
        for v in vs:
            if v in seen: continue
            o = {act(g, v) for g in Gfull}; seen |= o; rows.append((name, len(o)))
    other = [v for v in V if not any(v in f[1] for f in facs)]
    rows += [('fixe', 1)] * len(other)
    orb = collections.Counter(s for _, s in rows)
    fixed = sum(1 for _, s in rows if s == 1)
    hV = sum(s * math.log(s) for _, s in rows) / len(V)
    return rows, fixed, hV

def N5(K, facs):
    def canon(f):
        return tuple(min(tuple(sorted(act(g, v) for v in f if v in set(vs))) for g in Gfull) for _, vs, act, Gfull in facs) + (tuple(sorted(v for v in f if not any(v in set(x[1]) for x in facs))),)
    cnt = collections.Counter(canon(f) for f in K)
    sizes = sorted(cnt.values())
    Gord = math.prod(len(f[3]) for f in facs)
    return sizes, (sum(sizes) / len(sizes)) / Gord

def nsimplices(K):
    S = set()
    for f in K:
        f = tuple(sorted(f))
        for k in range(1, len(f) + 1): S.update(itertools.combinations(f, k))
    return len(S)

if __name__ == '__main__':
    res = {}
    for nm, (K, V) in OBJ.items():
        facs = factors_of(V)
        Gord = math.prod(len(f[3]) for f in facs)
        b = betti_z2(K)
        rows, fixed, hV = N4(V, facs)
        fsz, phi = N5(K, facs)
        ns = nsimplices(K)
        a = aut_order(V, K)
        N3 = {}
        for lab, gens in subgroup_list(facs):
            KH = fixed_complex(K, V, gens)
            N3[lab] = betti_z2(KH) if KH else 'vide'
        tot3 = sum(sum(x.values()) for x in N3.values() if x != 'vide')
        res[nm] = dict(v=len(V), kappa=len(K), simplexes=ns, G=Gord, aut=a, betti=b, fixes=fixed, hV=hV,
                       orbites_feuilles=sorted(s for _, s in rows), orbites_facettes=len(fsz), phi=phi, N3=N3, N3_total=tot3,
                       N3_vides=sum(1 for x in N3.values() if x == 'vide'))
        r = res[nm]
        print(f"{nm:6s} v={r['v']:2d} kappa={r['kappa']:4d} simplexes={ns:6d} |G|={Gord:6d} |Aut|={a:6d} fixes={fixed} hV={hV:.3f} "
              f"orb.fac={len(fsz)} phi={phi:.4f} N3: total={tot3} vides={r['N3_vides']} Betti={b}", flush=True)
    json.dump(res, open('chain_inv.json', 'w'), default=str)
