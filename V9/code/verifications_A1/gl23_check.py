"""A1-16 / 1.5 : R^det[GL(2,3)] ; det = sign of the permutation of the 9 points of F_3^2 (article).
Also checks the counts used in the proof of prop:spinsector (ii)."""
import itertools
from twisted_decomp import decompose, sign, total_dim
pts = [(a, b) for a in range(3) for b in range(3)]
mats = [m for m in itertools.product(range(3), repeat=4) if (m[0]*m[3]-m[1]*m[2]) % 3]
def perm(m):
    a, b, c, d = m
    return tuple(pts.index(((a*x+b*y) % 3, (c*x+d*y) % 3)) for x, y in pts)
G = sorted(perm(m) for m in mats)
det = {perm(m): (1 if (m[0]*m[3]-m[1]*m[2]) % 3 == 1 else -1) for m in mats}
assert all(det[g] == sign(g) for g in G), "det != sign on 9 points"
print("|GL(2,3)| =", len(G), "; det = sign of the permutation of the 9 points : True")
zc, parts = decompose(G, lambda g: det[g])
print("R^det[GL(2,3)] =", " + ".join(parts), "; dim", total_dim(parts))
def mul(m, n):
    a,b,c,d=m; e,f,g,h=n
    return ((a*e+b*g)%3,(a*f+b*h)%3,(c*e+d*g)%3,(c*f+d*h)%3)
I=(1,0,0,1); mI=(2,0,0,2)
sqI=[m for m in mats if mul(m,m)==I]; sqmI=[m for m in mats if mul(m,m)==mI]
d=lambda m: 1 if (m[0]*m[3]-m[1]*m[2])%3==1 else -1
print("#{g^2=I} =", len(sqI), " f(I) =", sum(map(d,sqI)), "; #{g^2=-I} =", len(sqmI), " f(-I) =", sum(map(d,sqmI)))
# does GL(2,3) contain a subgroup isomorphic to S4 ? (order-24 subgroups)
from itertools import combinations
