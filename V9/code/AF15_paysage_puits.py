"""AF15 -- Paysage revu (05_audit/PAYSAGE_REVU_14-11.md, criteres ecrits le 14/11/2026 a 14 h 10 avant calcul).
Fonctionnelle E' lexicographique, sans parametre :
   e_T (tests du triplet ; = 0 dans la classe K26)  >  e_R = rayon du graphe primal PG  >  e_V = nombre de sommets
   >  e_C = (rayon, diametre, excentricite moyenne) du graphe des facettes LG, puis du graphe d'incidence IG  >  e_S = (volume, kappa).
Mouvements : orbites (sous Aut = S4^3, 13 824) de suppression, reduction, agrandissement elementaire (<= 2 facettes absorbees), comme e5/a1 (A2).
Appartenance a la classe pour un mouvement qui ne monte pas : recouvrement connexe, Betti Z/2 egaux, cone de l'inclusion acyclique sur F2, |Aut| = 13 824
(A2, Prop. 'Orbital moves keeping the class' => bimodule, gamma, J' transportes ; Aut = G => sigma normalise Aut).
Usage :
   python3 AF15_paysage_puits.py cle                 # E' des representants connus
   python3 AF15_paysage_puits.py puits K4|V37|V32|V36|fichier.json      # Q1 : le depart est-il un puits ?
   python3 AF15_paysage_puits.py descente K4|... [pas_max]              # Q2 : descente de plus grande pente jusqu'a un puits
"""
import sys, os, json, time, itertools, collections
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'articleV_FHT21'); sys.path.insert(0, SRC)
from fht_tools import betti_z2, aut_order, build_K3, group_perms_K3
OBJ = os.path.join(HERE, '..', '..', '..', '04_objets')
K0, V = build_K3(os.path.join(SRC, 'V25_K3_pieces.json'), with_pendants=False)
n = len(V); VI = {v: i for i, v in enumerate(V)}
GRP = [[VI[g[v]] for v in V] for g in group_perms_K3(V)]          # 13 824 permutations d'indices
G_ORDER = len(GRP)
def tomask(f): return sum(1 << VI[v] for v in f)
def bits(m): return [i for i in range(n) if m >> i & 1]
def names(m): return [V[i] for i in bits(m)]
def img(p, m):
    r = 0
    for i in bits(m): r |= 1 << p[i]
    return r
def load(src):
    if src == 'K4': return [tomask(f) for f in K0]
    fn = {'V37': 'V37_FHT_K26volmin.json', 'V32': 'V32_FHT_K26etoile.json', 'V36': 'V36_FHT_K26vol.json'}.get(src)
    d = json.load(open(os.path.join(OBJ, fn) if fn else src))
    d = d['facettes'] if isinstance(d, dict) else d
    return [tomask(f) for f in d]
# ------------------------------------------------------------------ invariants
def incid(F): return np.array([[m >> i & 1 for i in range(n)] for m in F], dtype=np.float32)
def ecc(A):
    """excentricites exactes d'un graphe (matrice d'adjacence booleenne)"""
    N = A.shape[0]; Af = A.astype(np.float32); R = A | np.eye(N, dtype=bool); e = np.full(N, -1); k = 0
    e[R.all(axis=1)] = 1 if N > 1 else 0
    while (e < 0).any():
        k += 1
        if k > N: return None
        R = R | ((R.astype(np.float32) @ Af) > 0)
        e[(e < 0) & R.all(axis=1)] = k + 1
    return e
def triple(e): return (int(e.min()), int(e.max()), Fraction(int(e.sum()), len(e)))
def rad_PG(B):
    C = (B.T @ B) > 0; np.fill_diagonal(C, False); return int(ecc(C).min())
def LG(B):
    A = (B @ B.T) > 0; np.fill_diagonal(A, False)
    univ = (A | np.eye(A.shape[0], dtype=bool)).all(axis=1); u = int(univ.sum())
    m = A.shape[0]
    if u: return (1, 1 if u == m else 2, Fraction(u * 1 + (m - u) * 2, m)) if u < m else (1, 1, Fraction(1))
    return triple(ecc(A))
def IG(B):
    m = B.shape[0]; N = n + m; A = np.zeros((N, N), dtype=bool)
    A[:n, n:] = B.T > 0; A[n:, :n] = B > 0
    return triple(ecc(A))
def volume(F):
    S = set()
    for f in F:
        b = bits(f); sub = 0
        # sous-masques non vides
        s = f
        while s:
            S.add(s); s = (s - 1) & f
    return len(S)
class Key:
    """E' paresseuse : niveaux calcules a la demande"""
    def __init__(s, F): s.F = F; s.B = incid(F); s.c = {}
    def lvl(s, k):
        if k not in s.c:
            s.c[k] = {0: lambda: rad_PG(s.B), 1: lambda: n, 2: lambda: LG(s.B), 3: lambda: IG(s.B),
                      4: lambda: volume(s.F), 5: lambda: len(s.F)}[k]()
        return s.c[k]
    def full(s): return tuple(s.lvl(k) for k in range(6))
def cmp(a, b):
    for k in range(6):
        x, y = a.lvl(k), b.lvl(k)
        if x != y: return (-1 if x < y else 1), k
    return 0, 6
LEVEL = ['e_R (rayon PG)', 'e_V (sommets)', 'e_C (LG)', 'e_C (IG)', 'e_S (volume)', 'e_S (kappa)', 'egal']
# ------------------------------------------------------------------ mouvements d'orbite
def orbit_reps(F):
    seen, reps = set(), []
    for f in F:
        if f in seen: continue
        o = {img(p, f) for p in GRP}; seen |= o; reps.append(f)
    return reps
def moves(F):
    out, seen = [], set()
    for e in orbit_reps(F):
        cand = [('suppr', (e,))]
        if len(bits(e)) >= 3: cand += [('reduc', (e, 1 << i)) for i in bits(e)]
        cand += [('agrand', (e, 1 << i)) for i in range(n) if not e >> i & 1]
        for kind, data in cand:
            o = frozenset(tuple(img(p, x) for x in data) for p in GRP)
            if (kind, o) in seen: continue
            seen.add((kind, o)); out.append((kind, o))
    return out
def apply(F, kind, o):
    Fs = set(F)
    if kind == 'suppr':
        rem = {d[0] for d in o}; return [f for f in F if f not in rem], True
    if kind == 'reduc':
        rem = {d[0] for d in o}; keep = [f for f in F if f not in rem]
        new = {d[0] & ~d[1] for d in o}
        new = [g for g in new if not any(g & h == g for h in keep)]
        new = [g for g in new if not any(g != h and g & h == g for h in new)]
        return keep + new, True
    new = {d[0] | d[1] for d in o}
    for g in new:
        if sum(1 for f in F if f & g == f) > 2: return None, False
    keep = [f for f in F if not any(f & g == f for g in new)]
    return keep + sorted(new), False
def in_class(cur, nE, small_is_new, b0):
    if set().union(*[set(bits(f)) for f in nE]) != set(range(n)): return 'non couvrant'
    NN = [names(f) for f in nE]
    if betti_z2(NN) != b0: return 'Betti'
    small, big = (nE, cur) if small_is_new else (cur, nE)
    if betti_z2([names(f) for f in big] + [names(f) + ['*c*'] for f in small]): return 'cone'
    a = aut_order(V, NN)
    if a != G_ORDER: return 'Aut %d' % a
    return 'classe'
def fmt(k):
    if isinstance(k, tuple): return tuple(fmt(x) for x in k)
    return '%s (%.4f)' % (k, float(k)) if isinstance(k, Fraction) else k
# ------------------------------------------------------------------ modes
def P(*a, log=None):
    print(*a, flush=True)
    if log: print(*a, file=log, flush=True)
def puits(src, F=None, log=None, quiet=False):
    F = F or load(src); K = Key(F); t0 = time.time()
    b0 = betti_z2([names(f) for f in F])
    P(f"[{src}] facettes {len(F)} ; E' = {fmt(K.full())} ; Betti {b0} ; |Aut| {aut_order(V, [names(f) for f in F])}", log=log)
    M = moves(F); P(f"   {len(M)} classes de mouvements d'orbite ({time.time()-t0:.0f} s)", log=log)
    tally = collections.Counter(); down, flat = [], []
    for kind, o in M:
        nE, small = apply(F, kind, o)
        if nE is None: tally['agrand non elementaire'] += 1; continue
        Kn = Key(nE); c, lev = cmp(Kn, K)
        if c > 0: tally['monte (%s)' % LEVEL[lev]] += 1; continue
        verdict = in_class(F, nE, small, b0)
        if verdict != 'classe': tally['ne monte pas mais sort de la classe (%s) => e_T > 0, monte' % verdict] += 1; continue
        (down if c < 0 else flat).append((Kn, kind, o, nE, lev)); tally[('DESCEND (%s)' if c < 0 else 'EGAL') % LEVEL[lev]] += 1
    for k, v in sorted(tally.items()): P(f"   {k:70s} {v}", log=log)
    P(f"   duree {time.time()-t0:.0f} s", log=log)
    if not down and not flat: P(f"   VERDICT Q1 [{src}] : PUITS STRICT de E' pour les mouvements d'orbite", log=log)
    elif not down: P(f"   VERDICT Q1 [{src}] : PUITS PLAT ({len(flat)} mouvements a egalite)", log=log)
    else:
        best = min(down, key=lambda d: d[0].full())
        P(f"   VERDICT Q1 [{src}] : PAS UN PUITS ({len(down)} mouvements descendants ; meilleur : {best[1]}, E' = {fmt(best[0].full())})", log=log)
    return down, flat
def descente(src, maxsteps):
    F = load(src); tag = os.path.basename(src).replace('.json', '')
    log = open(os.path.join(HERE, f'AF15_descente_{tag}.log'), 'a')
    P(f"=== descente depuis {src} ({time.strftime('%d/%m %H:%M')})", log=log)
    for step in range(maxsteps):
        down, flat = puits(f'{tag} pas {step}', F, log=log)
        if not down:
            out = os.path.join(HERE, f'AF15_fond_{tag}.json'); json.dump([names(f) for f in F], open(out, 'w'))
            P(f"=== FOND atteint en {step} pas ; facettes ecrites dans {out}", log=log); return
        F = min(down, key=lambda d: d[0].full())[3]
        json.dump([names(f) for f in F], open(os.path.join(HERE, f'AF15_descente_{tag}_pas{step+1}.json'), 'w'))
    P("=== nombre maximal de pas atteint", log=log)
if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'cle'
    if mode == 'cle':
        for s in ['K4', 'V36', 'V37', 'V32']:
            F = load(s); K = Key(F); t = time.time()
            print(f"{s:4s} facettes {len(F):4d}  E' (e_R, e_V, LG, IG, volume, kappa) = {fmt(K.full())}  ({time.time()-t:.1f} s)", flush=True)
    elif mode == 'puits':
        log = open(os.path.join(HERE, f'AF15_puits_{os.path.basename(sys.argv[2]).replace(".json","")}.log'), 'w')
        puits(sys.argv[2], log=log)
    elif mode == 'descente': descente(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 200)
