#!/usr/bin/env python3
"""AUDIT -- Test des porteurs du triplet fini de K26 (structure reelle J' exacte, theoreme J de l'article 2).
Criteres ecrits AVANT calcul dans 05_audit/TEST_PORTEURS_TRIPLET.md (14/11/2026, 7 h 05, Guadeloupe).

Pour chaque complexe :
  * Aut : ensemble fort de generateurs (chaine de stabilisateurs, retour arriere elague par types de sommets et co-degres) ;
  * homologie reduite sur F_p (p = 65 521) avec representants de cycles (reduction de colonnes avec effacement et suivi des
    combinaisons) et action geometrique rho de Aut sur l'homologie (chaines orientees) ;
  * prealables : P1 sigma (echange L_ <-> R_) defini sur les sommets ; P2 sigma Aut sigma = Aut ; P3 Aut preserve {g0..g3} ;
    P4 dim H_pair = dim H_impair ;
  * C1 : existe-t-il A : H_impair -> H_pair inversible avec A s(g) rho_impair(sigma g sigma) = rho_pair(g) A pour tous les
    generateurs g (s(g) = signe de g sur {g0..g3}) ?  Espace des solutions calcule exactement sur F_p ; inversibilite testee sur
    3 combinaisons aleatoires.  Variante : solutions ne reliant que des degres voisins (d -> d +- 1) ;
  * C2 : caracteres de H_pair et H_impair sur les elements qui fixent tous les sommets hors L_/R_, classes par
    (type de cycle sur L_x0..x3, type de cycle sur R_x0..x3), compares a ceux de K26.
Controles : K26 (Definition), V37, V36, V32 doivent passer ; K26 avec s = 1 doit echouer a C1.
Python 3 + numpy.  N'ecrit aucun fichier.  Usage : python3 AF9_porteurs_triplet.py RACINE_FHT26 [noms ...]
"""
import itertools, json, os, sys, time, random
from collections import Counter
import numpy as np

P = 65521
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from AF3_K26_homologie_automorphismes import K26_definition   # construction d'apres la Definition de A2 (script d'audit AF3)

random.seed(20261114)

def sw(v):                                   # sigma : echange des etiquettes L_ <-> R_
    return v.replace('L_', '\x00').replace('R_', 'L_').replace('\x00', 'R_')

# ------------------------------------------------------------------ complexe
class Cx:
    def __init__(self, name, facets):
        self.name = name
        self.V = sorted(set().union(*[set(f) for f in facets]))
        self.ix = {v: i for i, v in enumerate(self.V)}
        self.F = [frozenset(self.ix[v] for v in f) for f in facets]
        self.Fset = set(self.F)
        S = set()
        for f in self.F:
            f = sorted(f)
            for r in range(1, len(f) + 1):
                S.update(itertools.combinations(f, r))
        self.S = {}
        for s in S: self.S.setdefault(len(s) - 1, []).append(s)
        for d in self.S: self.S[d].sort()
        self.idx = {d: {s: i for i, s in enumerate(L)} for d, L in self.S.items()}
        self.dmax = max(self.S)
        self.nsimp = sum(len(L) for L in self.S.values())

    def is_aut(self, perm):
        return all(frozenset(perm[v] for v in f) in self.Fset for f in self.F)

# ------------------------------------------------------------------ automorphismes
def aut_sgs(C):
    nv = len(C.V); F = C.F
    vtype = [tuple(sorted(Counter(len(f) for f in F if v in f).items())) for v in range(nv)]
    cod = [[tuple(sorted(Counter(len(f) for f in F if u in f and v in f).items())) for v in range(nv)] for u in range(nv)]
    def extend(img, order, pos):
        if pos == len(order):
            return dict(img) if all(frozenset(img[v] for v in f) in C.Fset for f in F) else None
        u = order[pos]; used = set(img.values())
        for y in range(nv):
            if y in used or vtype[y] != vtype[u]: continue
            if all(cod[u][w] == cod[y][img[w]] for w in img):
                img[u] = y
                r = extend(img, order, pos + 1)
                if r is not None: return r
                del img[u]
        return None
    gens = []; order = 1; fixed = {}
    base = sorted(range(nv), key=lambda v: (vtype[v], v))
    for b in base:
        rest = [v for v in base if v not in fixed and v != b]
        orbit = 0
        for y in range(nv):
            if y in fixed.values() or vtype[y] != vtype[b]: continue
            if not all(cod[b][w] == cod[y][fixed[w]] for w in fixed): continue
            img = dict(fixed); img[b] = y
            r = extend(img, rest, 0)
            if r is not None:
                orbit += 1
                if y != b: gens.append(tuple(r[i] for i in range(nv)))
        order *= orbit; fixed[b] = b
    return gens, order

def enumerate_group(gens, nv, cap=60000):
    e = tuple(range(nv)); seen = {e}; frontier = [e]
    while frontier:
        nf = []
        for p in frontier:
            for g in gens:
                q = tuple(g[p[i]] for i in range(nv))
                if q not in seen:
                    seen.add(q); nf.append(q)
                    if len(seen) > cap: return None
        frontier = nf
    return seen

# ------------------------------------------------------------------ algebre lineaire sur F_p (colonnes creuses)
def bcol(C, d, s):
    col = {}
    if d == 0: return {0: 1}                                     # bord augmente vers C_{-1}
    for i in range(len(s)):
        col[C.idx[d - 1][s[:i] + s[i + 1:]]] = 1 if i % 2 == 0 else P - 1
    return col

def reduce_cols(cols, clear=(), track=False):
    piv = {}; zeros = []
    for j, c in enumerate(cols):
        if j in clear: continue
        c = dict(c); v = {j: 1} if track else None
        while c:
            l = max(c)
            if l not in piv: break
            pc, inv, pv = piv[l]; fac = c[l] * inv % P
            for k, val in pc.items():
                nv_ = (c.get(k, 0) - fac * val) % P
                if nv_: c[k] = nv_
                else: c.pop(k, None)
            if track:
                for k, val in pv.items():
                    nv_ = (v.get(k, 0) - fac * val) % P
                    if nv_: v[k] = nv_
                    else: v.pop(k, None)
        if c:
            l = max(c); piv[l] = (c, pow(c[l], P - 2, P), v)
        else:
            zeros.append((j, v))
    return piv, zeros

class Homology:
    """H_d reduit sur F_p avec representants et coordonnees."""
    def __init__(self, C):
        self.C = C
        cols = {d: [bcol(C, d, s) for s in C.S[d]] for d in range(0, C.dmax + 1)}
        self.pivB = {}; rk = {}
        for d in range(0, C.dmax + 1):
            piv, _ = reduce_cols(cols[d]); rk[d] = len(piv)
            self.pivB[d - 1] = piv                              # pivots de bord dans C_{d-1}
        rk[C.dmax + 1] = 0; self.pivB.setdefault(C.dmax, {})
        n = {d: len(C.S[d]) for d in C.S}
        self.betti = {d: n[d] - rk[d] - rk[d + 1] for d in range(0, C.dmax + 1)}
        self.reps = {}; self.piv = {}
        for d, b in self.betti.items():
            if b <= 0: continue
            clear = set(self.pivB[d].keys())
            _, zeros = reduce_cols(cols[d], clear=clear, track=True)
            assert len(zeros) == b, (d, len(zeros), b)
            reps = [z[1] for z in zeros]                       # cycles (dict indice -> coeff)
            piv = {l: (c, inv, None) for l, (c, inv, _) in self.pivB[d].items()}
            for i, z in enumerate(reps):
                c = dict(z); co = np.zeros(b, dtype=np.int64); co[i] = 1
                while c:
                    l = max(c)
                    if l not in piv: break
                    pc, inv, pco = piv[l]; fac = c[l] * inv % P
                    for k, val in pc.items():
                        nv_ = (c.get(k, 0) - fac * val) % P
                        if nv_: c[k] = nv_
                        else: c.pop(k, None)
                    if pco is not None: co = (co - fac * pco) % P
                assert c, 'representant dans les bords'
                l = max(c); piv[l] = (c, pow(c[l], P - 2, P), co)
            self.reps[d] = reps; self.piv[d] = piv

    def coords(self, d, x):
        piv = self.piv[d]; b = self.betti[d]; co = np.zeros(b, dtype=np.int64); c = dict(x)
        while c:
            l = max(c)
            pc, inv, pco = piv[l]                               # KeyError = x n'est pas un cycle
            fac = c[l] * inv % P
            for k, val in pc.items():
                nv_ = (c.get(k, 0) - fac * val) % P
                if nv_: c[k] = nv_
                else: c.pop(k, None)
            if pco is not None: co = (co + fac * pco) % P
        return co

    def act(self, d, g, chain):
        C = self.C; out = {}
        for k, coef in chain.items():
            t = [g[i] for i in C.S[d][k]]
            inv = sum(1 for a in range(len(t)) for b_ in range(a + 1, len(t)) if t[a] > t[b_])
            j = C.idx[d][tuple(sorted(t))]
            val = (coef if inv % 2 == 0 else -coef) % P
            nv_ = (out.get(j, 0) + val) % P
            if nv_: out[j] = nv_
            else: out.pop(j, None)
        return out

    def mat(self, d, g):
        return np.array([self.coords(d, self.act(d, g, z)) for z in self.reps[d]], dtype=np.int64).T % P

    def degrees(self, parity):
        return [d for d, b in sorted(self.betti.items()) if b > 0 and d % 2 == parity]

    def block(self, parity, g):
        ds = self.degrees(parity); bs = [self.betti[d] for d in ds]; n = sum(bs)
        M = np.zeros((n, n), dtype=np.int64); o = 0
        for d, b in zip(ds, bs):
            M[o:o + b, o:o + b] = self.mat(d, g); o += b
        lab = [d for d, b in zip(ds, bs) for _ in range(b)]
        return M, lab

def nullspace_mod(M):
    """base (lignes) de {x : M x = 0} sur F_p"""
    M = (M % P).astype(np.int64).copy(); m, k = M.shape; pc = []; row = 0
    for col in range(k):
        if row >= m: break
        nz = np.nonzero(M[row:, col])[0]
        if len(nz) == 0: continue
        r = row + nz[0]
        if r != row: M[[row, r]] = M[[r, row]]
        M[row] = (M[row] * pow(int(M[row, col]), P - 2, P)) % P
        oth = np.nonzero(M[:, col])[0]; oth = oth[oth != row]
        if len(oth): M[oth] = (M[oth] - np.outer(M[oth, col], M[row]) % P) % P
        pc.append(col); row += 1
    pcs = set(pc); B = []
    for f in range(k):
        if f in pcs: continue
        x = np.zeros(k, dtype=np.int64); x[f] = 1
        for i, c in enumerate(pc): x[c] = (-M[i, f]) % P
        B.append(x)
    return np.array(B, dtype=np.int64).reshape(len(B), k)

def rank_mod(M):
    return M.shape[1] - nullspace_mod(M).shape[0]

def intertwiners(pairs, labE, labO, shift_only=False):
    """pairs : liste de (M_g, N_g) ; cherche A (nE x nO) avec A N_g = M_g A.  Renvoie la base (k, nE, nO)."""
    nE, nO = len(labE), len(labO)
    allowed = [(r, c) for r in range(nE) for c in range(nO) if (not shift_only or abs(labE[r] - labO[c]) == 1)]
    k = len(allowed)
    Bas = np.zeros((k, nE * nO), dtype=np.int64)
    for j, (r, c) in enumerate(allowed): Bas[j, r * nO + c] = 1
    for M, N in pairs:
        if Bas.shape[0] == 0: break
        L = np.zeros((Bas.shape[0], nE * nO), dtype=np.int64)
        for j in range(Bas.shape[0]):
            A = Bas[j].reshape(nE, nO)
            L[j] = ((A @ N) % P - (M @ A) % P).reshape(-1) % P
        Cn = nullspace_mod(L.T)                                 # combinaisons c avec sum c_j L_j = 0
        Bas = (Cn @ Bas) % P if Cn.shape[0] else np.zeros((0, nE * nO), dtype=np.int64)
    return Bas.reshape(Bas.shape[0], nE, nO)

def has_invertible(Bas, n, tries=3):
    if Bas.shape[0] == 0 or Bas.shape[1] != Bas.shape[2]: return False
    for _ in range(tries):
        c = np.array([random.randrange(P) for _ in range(Bas.shape[0])], dtype=np.int64)
        A = np.tensordot(c, Bas, axes=1) % P
        if rank_mod(A) == n: return True
    return False

def signed(x):
    x = int(x) % P
    return x - P if x > P // 2 else x

def ctype(perm, pts):
    if any(perm[p] not in pts for p in pts): return None
    seen = set(); out = []
    for p in pts:
        if p in seen: continue
        l = 0; q = p
        while q not in seen: seen.add(q); q = perm[q]; l += 1
        out.append(l)
    return tuple(sorted(out, reverse=True))

# ------------------------------------------------------------------ test d'un complexe
def test(name, facets, ref_chars=None, s_trivial=False, want_chars=False):
    t0 = time.time(); C = Cx(name, facets); out = dict(name=name, v=len(C.V), facets=len(C.F), simplices=C.nsimp)
    gens, order = aut_sgs(C); out['aut'] = order; out['ngens'] = len(gens)
    nv = len(C.V); V = C.V
    gidx = [C.ix[v] for v in ('g0', 'g1', 'g2', 'g3') if v in C.ix]
    # P1 -- P3
    P1 = all(sw(v) in C.ix for v in V); out['P1'] = P1
    sig = [C.ix[sw(v)] for v in V] if P1 else None
    P2 = P1 and all(C.is_aut([sig[g[sig[i]]] for i in range(nv)]) for g in gens); out['P2'] = P2
    P3 = len(gidx) > 0 and all(set(g[i] for i in gidx) == set(gidx) for g in gens); out['P3'] = P3
    # action sur {g0..g3} : ordre du groupe image
    H = Homology(C); out['betti'] = {d: b for d, b in H.betti.items() if b}
    nE = sum(b for d, b in H.betti.items() if d % 2 == 0); nO = sum(b for d, b in H.betti.items() if d % 2 == 1)
    P4 = (nE == nO and nE > 0); out['P4'] = P4; out['dims'] = (nE, nO)
    out['C1'] = out['C1_shift'] = None; out['dimC1'] = None
    if P1 and P2 and P3 and P4:
        def sgn_g(g):
            pos = {p: k for k, p in enumerate(gidx)}; perm = [pos[g[p]] for p in gidx]
            inv = sum(1 for a in range(len(perm)) for b in range(a + 1, len(perm)) if perm[a] > perm[b])
            return -1 if inv % 2 else 1
        pairs = []
        for g in gens:
            M, labE = H.block(0, g)
            h = tuple(sig[g[sig[i]]] for i in range(nv))
            N, labO = H.block(1, h)
            s = 1 if s_trivial else sgn_g(g)
            pairs.append((M, (s * N) % P))
        Bas = intertwiners(pairs, labE, labO)
        out['dimC1'] = Bas.shape[0]; out['C1'] = has_invertible(Bas, nE)
        Bs = intertwiners(pairs, labE, labO, shift_only=True)
        out['C1_shift'] = has_invertible(Bs, nE)
        # controle de J' construit : U^2 = 1 et equivariance sur les generateurs (par construction ; verifie)
        if out['C1']:
            c = np.array([random.randrange(P) for _ in range(Bas.shape[0])], dtype=np.int64)
            A = np.tensordot(c, Bas, axes=1) % P
            ok = all(np.all((A @ N - M @ A) % P == 0) for M, N in pairs); out['check_eq'] = ok
    # C2 : caracteres sur les elements L/R purs
    chars = None
    if want_chars or ref_chars is not None:
        G = enumerate_group(gens, nv)
        others = [C.ix[v] for v in V if not (v.startswith('L_') or v.startswith('R_'))]
        Lx = [C.ix[f'L_x{i}'] for i in range(4) if f'L_x{i}' in C.ix]; Rx = [C.ix[f'R_x{i}'] for i in range(4) if f'R_x{i}' in C.ix]
        bylab = {}
        if G is not None:
            for g in G:
                if any(g[o] != o for o in others): continue
                lab = (ctype(g, Lx), ctype(g, Rx))
                if None in lab: continue
                bylab.setdefault(lab, []).append(g)
        chars = {}; incons = False
        for lab, els in sorted(bylab.items()):
            vals = set()
            for g in els[:3]:
                ME, _ = H.block(0, g); MO, _ = H.block(1, g)
                vals.add((signed(np.trace(ME)), signed(np.trace(MO))))
            if len(vals) > 1: incons = True
            chars[lab] = sorted(vals)[0]
        out['C2_incoherent'] = incons
        if ref_chars is not None:
            out['C2'] = (set(chars) == set(ref_chars)) and all(chars[k] == ref_chars[k] for k in ref_chars) and not incons
            out['C2_labels'] = (len(chars), len(ref_chars))
    out['t'] = round(time.time() - t0, 1)
    return out, chars

def verdict(o):
    if not o['P1']: return 'NON PORTEUR (P1 : sigma non defini)'
    if not o['P2']: return 'NON PORTEUR (P2 : sigma Aut sigma != Aut)'
    if not o['P3']: return 'NON PORTEUR (P3 : {g0..g3} non preserve)'
    if not o['P4']: return 'NON PORTEUR (P4 : dimensions paire/impaire %s)' % (o['dims'],)
    if not o['C1']: return 'NON PORTEUR (C1 : pas de J\' exact)'
    if o.get('C2'): return 'PORTEUR DU TRIPLET DE K26' + ('' if o['aut'] >= 13824 else ' (groupe plus petit que S4^3)')
    return 'STRUCTURE REELLE, BIMODULE DIFFERENT'

def show(o, v=None):
    print(f"{o['name']:12s} v={o['v']:2d} facettes={o['facets']:4d} simplexes={o['simplices']:6d} |Aut|={o['aut']:6d} "
          f"Betti={o['betti']} P1..P4={int(o['P1'])}{int(o['P2'])}{int(o['P3'])}{int(o['P4'])} dimC1={o['dimC1']} C1={o['C1']} "
          f"C1(d+-1)={o['C1_shift']} C2={o.get('C2')} ({o['t']} s)" + (f"  -> {v}" if v else ''), flush=True)

if __name__ == '__main__':
    root = sys.argv[1]
    only = set(sys.argv[2:])
    cat = os.path.join(root, '03_code', 'catalogues'); obj = os.path.join(root, '04_objets')
    print('=== Controles positifs (classe de K26) et controle negatif')
    o, ref = test('K26', K26_definition(), want_chars=True); o['C2'] = True; show(o, verdict(o))
    print('    caracteres de reference (type L_x, type R_x) -> (tr pair, tr impair) :', ref)
    okctl = verdict(o).startswith('PORTEUR')
    on, _ = test('K26 s=1', K26_definition(), s_trivial=True); show(on)
    neg_ok = (on['C1'] is False)
    print('    controle negatif (doit echouer a C1) :', 'echoue, conforme' if neg_ok else 'PASSE : TEST INVALIDE')
    for fn, nm in (('V37_FHT_K26volmin.json', 'V37'), ('V36_FHT_K26vol.json', 'V36'), ('V32_FHT_K26etoile.json', 'V32')):
        d = json.load(open(os.path.join(obj, fn))); fac = d if isinstance(d, list) else d.get('facets', d.get('facettes'))
        o, _ = test(nm, fac, ref_chars=ref); v = verdict(o); show(o, v); okctl &= v.startswith('PORTEUR')
    print('    controles positifs :', 'tous passent' if okctl else 'ECHEC : TEST INVALIDE')
    print('=== Echantillon')
    et = json.load(open(os.path.join(cat, 'V145_etoiles_facettes.json')))['etoiles']
    sample = [(k, et[k]['facets']) for k in ('K_27_546', 'K_27_559', 'K_27_1021', 'K_27_1084', 'K_28_767', 'K_28_789', 'K_36_44')]
    v129 = {json.loads(l)['id']: json.loads(l)['facets'] for l in open(os.path.join(cat, 'V129_candidats_facettes.jsonl'))}
    sample += [(k, v129[k]) for k in ('V129-011', 'V129-017', 'V129-018')]
    for nm, fac in sample:
        if only and nm not in only: continue
        o, _ = test(nm, fac, ref_chars=ref); show(o, verdict(o))
    print('VALIDITE DU TEST :', 'valide' if (okctl and neg_ok) else 'INVALIDE')
