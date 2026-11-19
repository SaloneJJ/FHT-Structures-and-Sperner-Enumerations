"""Character norms used in A1 (prop:spinsector): (1) Lambda^2 of the standard representation W of S5
is irreducible; (2) the restriction to 2O (binary octahedral, 48 unit quaternions) of the spin-3/2
representation Sym^3(C^2) is irreducible; the faithful irreducible representations of GL(2,3) and 2O
have degrees 2,2,4 (8 classes, quotient S4)."""
import itertools, numpy as np
# (1) S5
S5 = list(itertools.permutations(range(5)))
def fix(p): return sum(1 for i in range(5) if p[i] == i)
def chiW(p): return fix(p) - 1
def sq(p): return tuple(p[p[i]] for i in range(5))
# chi_{Lambda^2 W}(g) = (chi(g)^2 - chi(g^2))/2
n1 = sum(((chiW(p)**2 - chiW(sq(p)))/2)**2 for p in S5)/120
print("S5: <chi_Lambda2W, chi_Lambda2W> =", n1, "; dim =", (16-4)/2 if False else (chiW(tuple(range(5)))**2-chiW(tuple(range(5))))//2)
# (2) 2O
s = 1/np.sqrt(2); units = []
base = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
for b in base:
    for sg in (1,-1): units.append(tuple(sg*x for x in b))
for sgn in itertools.product((1,-1), repeat=4): units.append(tuple(0.5*x for x in sgn))
for i,j in itertools.combinations(range(4),2):
    for a,b in itertools.product((1,-1), repeat=2):
        q=[0,0,0,0]; q[i]=a*s; q[j]=b*s; units.append(tuple(q))
assert len(units) == 48
def chi32(q):  # spin j=3/2: chi = sin(4 theta)/sin(theta), q = cos(theta) + ...
    th = np.arccos(np.clip(q[0], -1, 1))
    return 4*np.cos(th)**3*1 if abs(np.sin(th))<1e-12 and False else (4.0*np.sign(np.cos(th))**3 if abs(np.sin(th))<1e-12 else np.sin(4*th)/np.sin(th))
n2 = sum(chi32(q)**2 for q in units)/48
print("2O: <chi_3/2, chi_3/2> =", round(n2, 12))
print("degrees: 24 = a^2+b^2+c^2 with 1<=a<=b<=c :", [t for t in itertools.combinations_with_replacement(range(1,5),3) if sum(x*x for x in t)==24])
