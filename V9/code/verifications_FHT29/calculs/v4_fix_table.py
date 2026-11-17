# Verification (Fable) : table Fix(lambda) de l'annexe de A1 (clutters connexes invariants)
import re, itertools, math, sys
from fractions import Fraction
tex = open(sys.argv[1], encoding='utf-8').read().split('\n')
start = next(i for i,l in enumerate(tex) if 'label{tab:fix}' in l); rows = {}; m = None
for l in tex[start:start+90]:
    mm = re.match(r'\s*(\d*)\s*&\s*\$\(([\d,]+)\)\$\s*&\s*([\d\\, ]+)\\\\', l)
    if mm:
        if mm.group(1): m = int(mm.group(1))
        lam = tuple(int(x) for x in mm.group(2).split(','))
        rows[(m, lam)] = int(re.sub(r'[^\d]', '', mm.group(3)))
def z(lam):
    r = 1
    for k in set(lam): c = lam.count(k); r *= k**c * math.factorial(c)
    return r
print("lignes lues :", len(rows))
for m in sorted({k[0] for k in rows}):
    lams = [k[1] for k in rows if k[0] == m]
    S = sum(Fraction(rows[(m,l)], z(l)) for l in lams)
    print(f"m={m}: {len(lams)} partitions (attendu {sum(1 for _ in __import__('sympy').utilities.iterables.partitions(m))}), Burnside = {S}")
# force brute m <= 5
def clutters(m):
    subs = [s for s in range(1, 1<<m)]
    subs.sort(key=lambda s: bin(s).count('1'))
    out = []
    def rec(i, chosen):
        if i == len(subs):
            out.append(tuple(chosen)); return
        s = subs[i]
        rec(i+1, chosen)
        if all((s & c) != c and (s & c) != s for c in chosen): rec(i+1, chosen+[s])
    rec(0, [])
    return out
def connected_cover(cl, m):
    if not cl: return False
    u = 0
    for e in cl: u |= e
    if u != (1<<m)-1: return False
    comp = cl[0]; rest = list(cl[1:]); ch = True
    while ch:
        ch = False
        for e in rest[:]:
            if e & comp: comp |= e; rest.remove(e); ch = True
    return not rest
def cyc_perm(lam):
    p = []; a = 0
    for k in lam: p += [a+((i+1)%k) for i in range(k)]; a += k
    return p
for m in (2,3,4,5):
    cc = [frozenset(c) for c in clutters(m) if connected_cover(c, m)]
    S = set(cc)
    for (mm, lam), val in sorted(rows.items()):
        if mm != m: continue
        p = cyc_perm(lam)
        img = lambda e: sum(1<<p[i] for i in range(m) if e>>i & 1)
        fix = sum(1 for c in cc if frozenset(img(e) for e in c) == c)
        print(f"  m={m} lambda={lam}: table {val}, force brute {fix}", "OK" if fix == val else "ECART")

# --- suite des FHT par l'equation d'especes F(x) = x + sum_m Z_m(F(x), F(x^2), ...), a partir de la table
def read_all(path):
    tex = open(path, encoding='utf-8').read().split('\n')
    start = next(i for i,l in enumerate(tex) if 'label{tab:fix}' in l); rows = {}; m = None
    for l in tex[start:start+95]:
        if 'end{longtable}' in l: break
        mm = re.match(r'\s*(\d*)\s*&\s*\$\(([\d,]+)\)\$\s*&\s*(.*)$', l)
        if mm:
            if mm.group(1): m = int(mm.group(1))
            lam = tuple(int(x) for x in mm.group(2).split(','))
            v = re.sub(r'\$\^\\d?dagger\$', '', mm.group(3)); v = re.sub(r'\\begin\{tabular\}\[t\]\{@\{\}r@\{\}\}|\\end\{tabular\}', '', v)
            rows[(m, lam)] = int(re.sub(r'[^\d]', '', v))
    return rows
R = read_all(sys.argv[1]); N = 8
print("Burnside m=7, 8 :", [sum(Fraction(R[k], z(k[1])) for k in R if k[0] == m) for m in (7, 8)])
def mul(a, b):
    c = [Fraction(0)]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if i+j <= N and y: c[i+j] += x*y
    return c
F = [Fraction(0)]*(N+1); F[1] = Fraction(1)
for it in range(N):
    G = [Fraction(0)]*(N+1); G[1] = Fraction(1)
    for (m, lam), val in R.items():
        term = [Fraction(0)]*(N+1); term[0] = Fraction(val, z(lam))
        for k in lam:
            Fk = [Fraction(0)]*(N+1)
            for i,x in enumerate(F):
                if i*k <= N: Fk[i*k] = x
            term = mul(term, Fk)
        G = [g+t for g,t in zip(G, term)]
    F = G
print("suite des FHT, n = 1..8 :", [int(x) for x in F[1:]], "| somme n<=5 :", int(sum(F[1:6])))
