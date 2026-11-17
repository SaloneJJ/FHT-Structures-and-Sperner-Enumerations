# (c1) recount of K^(4) from the definition (independent rewrite; no file is written)
import itertools
from collections import Counter
def A(s):
    v = lambda n: s + n
    Z = frozenset([v('z')] + [v('x%d' % i) for i in range(4)])
    F = [Z, frozenset([v('p0'), v('p1'), v('p2')]), frozenset([v('q')])]
    F += [frozenset([v('f'), v('p%d' % i), v('p%d' % j)]) for i, j in ((0, 1), (0, 2), (1, 2))]
    pair = {0: ((0, 1), (2, 3)), 1: ((0, 2), (1, 3)), 2: ((0, 3), (1, 2))}
    for k, halves in pair.items():
        for a, b in halves:
            F.append(frozenset([v('p%d' % k), v('z'), v('x%d' % a), v('x%d' % b)]))
    return F
def B(s):
    v = lambda n: s + n
    Z = [v('z')] + [v('x%d' % i) for i in range(4)]
    return [frozenset(Z)] + [frozenset([v('w'), z]) for z in Z]
G0 = [frozenset(['g%d' % i]) for i in range(4)]
K4 = [frozenset(['g%d' % i, 'g%d' % j]) for i, j in itertools.combinations(range(4), 2)]
def join(*fams):
    return [frozenset().union(*c) for c in itertools.product(*fams)]
I = join(A('L'), B('R'), G0); II = join(B('L'), A('R'), K4)
print('facets of A, B:', len(A('L')), len(B('L')), '| generating faces of I, II:', len(I), len(II))
allf = set(I) | set(II)
maxf = [f for f in allf if not any(f < g for g in allf)]
simp = set()
for f in maxf:
    fs = sorted(f)
    for r in range(1, len(fs) + 1):
        simp.update(itertools.combinations(fs, r))
fv = Counter(len(s) for s in simp)
V = set().union(*maxf)
print('vertices', len(V), '| edges', fv[2], '| maximal faces', len(maxf), '| simplices', len(simp))
print('f-vector', [fv[k] for k in sorted(fv)])
print('reduced Euler characteristic', -1 + sum((-1) ** (k - 1) * fv[k] for k in fv), ' expected -12+36-36+12 = 0')
absorbed = [f for f in I if f not in maxf]
print('faces of I absorbed by II:', len(absorbed), sorted(map(len, absorbed)))
# same count for X_grand maximal faces (definition of the article): 32*6*4 + 6*32*6
print('X_grand facets by formula:', 32 * 6 * 4 + 6 * 32 * 6, '| X_GS:', (32 * 6 * 4 + 6 * 32 * 6) * 21, '| vertices', 2 * (21 + 6 - 1) + 4, 56 + 24)
print('orders:', (48 * 120) ** 2 * 24, 6 * (48 * 120) ** 2 * 24, (16 * 120) ** 2 * 24, (16 * 24) ** 2 * 24, (12 * 24) ** 2 * 24, 24 ** 3)
print('dims:', 16 + 32, 64 + 128, 4 + 4 + 32, 4 + 4 + 2 + 18, 2 + 4 + 2 + 18, 2 + 4 + 18)
