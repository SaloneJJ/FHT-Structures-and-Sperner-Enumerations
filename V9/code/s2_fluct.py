"""V40 / S2 : fluctuations internes du triplet fini (1 generation, modele unitairement equivalent a H(K26))
(a) algebre A_F (gauche) ; (b) algebre de tous les automorphismes cote L et R : A_F (x) A_F° (gauche et droite independantes),
sans condition d'ordre un. Support des fluctuations par blocs (composantes) : Higgs, Majorana (sigma), autres."""
import contextlib, io, sys
import numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('/tmp/claude-0/-home-claude/94269512-dfd2-59fb-9fc2-b4a8e01f5b00/scratchpad/dspace.py').read().split("# V4 : commute with left rho(S4)")[0])
rng = np.random.default_rng(7)
X = left(1j * (central('1') - central('s')) + gprod(central('2'), u((1, 0, 2, 3))))
M1c = lin_constraints(lambda D: np.concatenate([order1(D), (D @ X - X @ D).ravel()]))
Bc, s_ = complex_basis(M1c)
def herm_J(Bs):
    sol = real_dim(Bs); return sol[0] if isinstance(sol, tuple) else sol
sols = herm_J(Bc)
print("famille photon (1 gen) dim_R =", len(sols))
D0 = sum(rng.normal() * S for S in sols)
names = [c[3] for c in comps]
def support(Ms, tol=1e-8):
    out = {}
    for c1 in range(12):
        for c2 in range(12):
            n = max(np.abs(M[blk(c1), blk(c2)]).max() for M in Ms)
            if n > tol: out[(names[c1], names[c2])] = n
    return out
print("support de D0 :", sorted(support([D0]).keys()))
def vecR(M): return np.concatenate([M.real.ravel(), M.imag.ravel()])
def span_rank(Ms):
    F = np.array([vecR(M) for M in Ms]); return np.linalg.matrix_rank(F, tol=1e-7)
AFb = AF_basis
def rnd_af(): return sum(rng.normal() * c for c in AFb)
# (a) A_F
forms = []
for a in AFb:
    for b in AFb:
        La, Lb = left(a), left(b); w = La @ (D0 @ Lb - Lb @ D0); forms.append(w + w.conj().T)
fl = [A + Jop(A) for A in forms]
print("(a) A_F : dim_R Omega^1 (herm) =", span_rank(forms), "; A + JAJ^-1 :", span_rank(fl))
print("    support des fluctuations :", sorted(support(fl).keys()))
# (b) A_F (x) A_F°
def rnd_C(): return left(rnd_af()) @ right(rnd_af())
formsC = []
for _ in range(1500):
    a, b = rnd_C(), rnd_C(); w = a @ (D0 @ b - b @ D0); formsC.append(w + w.conj().T)
rC = span_rank(formsC)
print("(b) A_F x A_F° : dim_R Omega^1 (herm) =", rC)
sup = support(formsC)
print("    blocs touches :", len(sup))
for k in sorted(sup): print("     ", k)
# compare with D-space without order one (odd, self-adjoint, J) : which blocks are allowed at all
V1 = herm_J([Dfrom(np.eye(nU)[k]) for k in range(nU)])
print("D impairs autoadjoints J-symetriques sans ordre un : dim_R =", len(V1))

# (c) D sans ordre un mais invariant photon : espace et fluctuations par A_F (x) J A_F J^-1 (CCvS)
M2c = lin_constraints(lambda D: (D @ X - X @ D).ravel())
Bc2, _ = complex_basis(M2c)
sol2 = herm_J(Bc2)
print("(c) D impairs, J, photon, SANS ordre un : dim_R =", len(sol2))
extra = support(sol2); base = support(sols)
print("    blocs permis en plus de la famille d'ordre un :", sorted(set(extra) - set(base)))
D1 = sum(rng.normal() * S for S in sol2)
formsD1 = []
for _ in range(1500):
    a, b = rnd_C(), rnd_C(); w = a @ (D1 @ b - b @ D1); formsD1.append(w + w.conj().T)
print("    fluctuations (A_F x A_F°) d'un D generique sans ordre un : dim_R =", span_rank(formsD1))
s1 = support(formsD1)
print("    blocs touches en plus des Yukawa :", sorted(set(s1) - set(sup)))
print("    bloc Majorana (nuR,nubarR) touche :", ('nuR', 'nubarR') in s1)
