"""V353 : ce qui pourrait fixer l'angle φ. Partie algèbre (L1–L3, b1, c1). Critères : 01_synthese/V353 §0–1 (écrits avant calcul).
numpy seul. Les relevés de O_h sont ceux de V331 (rvA_common.py, à placer dans le même dossier).   JJ SALONE, FHT 2026"""
import numpy as np, itertools as it
import rvA_common as R
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
X = [k3(sx, I2, I2), k3(sz, sx, I2), k3(sz, sz, sx)]; Y = [k3(sy, I2, I2), k3(sz, sy, I2), k3(sz, sz, sy)]
G5 = sum(X)/np.sqrt(3); A = (X[0]-X[1])/np.sqrt(2); B = (X[0]+X[1]-2*X[2])/np.sqrt(6)
Y123 = Y[0]@Y[1]@Y[2]; G4 = G5@Y123; F = 1j*A@B; FG = F@G5
I8 = np.eye(8)
def O(phi): return (np.cos(phi)*A + np.sin(phi)*B)@Y123
def comm(a, b): return np.abs(a@b - b@a).max()
def acomm(a, b): return np.abs(a@b + b@a).max()
from scipy.linalg import expm
print("== L1 : F = iAB commute avec Γ5, Y_d, Γ4, FΓ5 ; anticommute avec A, B, O_φ")
print("  [F,Γ5] %.1e  max[F,Y_d] %.1e  [F,Γ4] %.1e  [F,FΓ5] %.1e" % (comm(F, G5), max(comm(F, y) for y in Y), comm(F, G4), comm(F, FG)))
print("  {F,A} %.1e  {F,B} %.1e  {F,O_0.7} %.1e   F²=1 : %s" % (acomm(F, A), acomm(F, B), acomm(F, O(0.7)), np.allclose(F@F, I8)))
for th in (0.3, 1.1):
    U = expm(1j*th*F); d = [np.abs(U@O(p)@U.conj().T - O(p + s*2*th)).max() for p in (0.0, 0.9) for s in (+1, -1)]
    print(f"  θ={th}: |U O_φ U† − O_(φ+2θ)| = {d[0]:.1e}, {d[2]:.1e} ; |… − O_(φ−2θ)| = {d[1]:.1e}, {d[3]:.1e}")
print("== L2 : O_φ = O_0 (cos φ − i sin φ F)")
print("  écart max sur φ ∈ [0,2π) :", "%.1e" % max(np.abs(O(p) - O(0)@(np.cos(p)*I8 - 1j*np.sin(p)*F)).max() for p in np.linspace(0, 6.28, 13)))
w, V = np.linalg.eigh(F); Pp = V[:, w > 0]; Pm = V[:, w < 0]
for p in (0.0, 0.5, 1.2):
    blk = Pm.conj().T@O(p)@Pp; ref = Pm.conj().T@O(0)@Pp
    ph = np.angle(np.trace(ref.conj().T@blk)); print(f"  φ = {p:.2f} : bloc ⟨f=−1|O_φ|f=+1⟩ = e^(i·{ph:+.4f}) × bloc(φ=0) ; écart {np.abs(blk - np.exp(1j*ph)*ref).max():.1e}")
print("== c1 : Γ4 et U_I (inversion) vis-à-vis de O_φ et F")
print("  {Γ4,O_φ} %.1e  [Γ4,F] %.1e  {FΓ5,O_φ} %.1e" % (acomm(G4, O(0.7)), comm(G4, F), acomm(FG, O(0.7))))
print("== L3 / b1 : relevés de O_h (V331) : action sur (A·Y123, B·Y123), sur FΓ5, sur Γ4")
rows = []
for g in R.ohgroup():
    U, ok = R.lift(g)
    if U is None: continue
    det, keepT = R.classify(g)
    M = np.array([[np.real(np.trace((U@O(a)@U.conj().T)@O(b).conj().T))/8 for a in (0, np.pi/2)] for b in (0, np.pi/2)])
    sFG = np.real(np.trace(U@FG@U.conj().T@FG))/8; sG4 = np.real(np.trace(U@G4@U.conj().T@G4))/8; sG5 = np.real(np.trace(U@G5@U.conj().T@G5))/8
    rot = det == 1
    rows.append((det, keepT, round(np.trace(g)), np.round(M, 4), round(sFG, 4), round(sG4, 4), round(sG5, 4)))
print("  nombre de relevés :", len(rows))
cls = {}
for det, kT, tr, M, sFG, sG4, sG5 in rows:
    key = (det, kT, tr, round(np.linalg.det(M), 3), round(np.trace(M), 3), sFG, sG4, sG5); cls[key] = cls.get(key, 0) + 1
print("  (det g, g∈T?, tr g, det M_E, tr M_E, FΓ5→±, Γ4→±, Γ5→±) : multiplicité")
for k, v in sorted(cls.items()): print("   ", k, ":", v)
# groupe du modèle (FΓ5 conservé) : image dans GL(E)
keepFG = [r for r in rows if r[4] > 0.99]
print("  éléments qui conservent FΓ5 :", len(keepFG), "; images dans GL(E) : det M_E =", sorted(set(round(np.linalg.det(r[3]), 3) for r in keepFG)),
      "; tr M_E =", sorted(set(round(np.trace(r[3]), 3) for r in keepFG)))
print("  → une réflexion de E (det M_E = −1) existe dans le groupe du modèle :", any(np.linalg.det(r[3]) < -0.5 for r in keepFG))
# ---- b1 : caractère de S4 sur H~1 de la pièce A (V22) et comparaison aux axes du cube
import json, os
pa = '/mnt/user-data/uploads/FHT_22/02_objets/V22_FHT_biS4_pieces.json'
if os.path.exists(pa):
    Ap = json.load(open(pa))['A']
else:
    Ap = [["z","x0","x1","x2","x3"],["z","p0","p1"],["f2","p0","p1"],["z","p0","p2"],["f2","p0","p2"],["z","p1","p2"],["f2","p1","p2"],
          ["p0","z","x0","x1"],["p0","z","x2","x3"],["p1","z","x0","x2"],["p1","z","x1","x3"],["p2","z","x0","x3"],["p2","z","x1","x2"],["u","p0"],["u","p1"],["u","p2"],["q"]]
simp = set()
for f in Ap:
    for r in range(1, len(f) + 1):
        for c in it.combinations(sorted(f), r): simp.add(c)
byd = {}
for s in simp: byd.setdefault(len(s) - 1, []).append(s)
for d in byd: byd[d].sort()
idx = {d: {s: i for i, s in enumerate(byd[d])} for d in byd}
def bd(d):
    M = np.zeros((len(byd[d-1]), len(byd[d])))
    for j, s in enumerate(byd[d]):
        for i in range(len(s)): M[idx[d-1][s[:i] + s[i+1:]], j] = (-1)**i
    return M
D1, D2 = bd(1), bd(2)
# H1 harmonique : ker D1 ∩ ker D2^T
L1 = D1.T@D1 + D2@D2.T; w1, V1 = np.linalg.eigh(L1); H1 = V1[:, w1 < 1e-9]
print("== b1 : pièce A (V22) : dim H1 =", H1.shape[1])
PAIR = {0: ((0, 1), (2, 3)), 1: ((0, 2), (1, 3)), 2: ((0, 3), (1, 2))}
def act(p, v):
    if v[0] in 'xy' and v[1:].isdigit(): return v[0] + str(p[int(v[1:])])
    if v[0] == 'p':
        a, b = PAIR[int(v[1])]; img = tuple(sorted([tuple(sorted((p[a[0]], p[a[1]]))), tuple(sorted((p[b[0]], p[b[1]])))]))
        return 'p' + str([k for k, pr in PAIR.items() if tuple(sorted(pr)) == img][0])
    return v
def perm_sign(seq):
    s = 1; seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]: s = -s
    return s
def rep1(p):
    M = np.zeros((len(byd[1]), len(byd[1])))
    for j, s in enumerate(byd[1]):
        im = [act(p, v) for v in s]; srt = tuple(sorted(im)); M[idx[1][srt], j] = perm_sign([srt.index(v) for v in im])
    return M
def cyc(p):
    seen, out = set(), []
    for i in range(4):
        if i in seen: continue
        l = 0; j = i
        while j not in seen: seen.add(j); j = p[j]; l += 1
        out.append(l)
    return tuple(sorted(out))
chi = {}
for p in it.permutations(range(4)):
    chi.setdefault(cyc(p), round(np.trace(H1.T@rep1(p)@H1), 6))
print("  caractère de S4 sur H1(A) par type de cycle :", chi, " (E de S4 : (1,1,1,1):2, (1,1,2):0, (2,2):2, (1,3):−1, (4):0)")
# axes du cube = appariements des quatre grandes diagonales
diag = [np.array(v) for v in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
print("  axes ↔ appariements des diagonales (coins du tétraèdre L : (+++), (+−−), (−+−), (−−+)) :")
for k, (a, b) in PAIR.items():
    ax = [d for d in range(3) if (diag[a[0]] + diag[a[1]])[d] != 0][0]
    print(f"    p{k} = {{{a[0]}{a[1]}|{b[0]}{b[1]}}} ↔ axe {'xyz'[ax]} (milieu de l'arête x{a[0]}x{a[1]} = {(diag[a[0]] + diag[a[1]])//2})")
# stabilisateur d'un coin (S3 = Stab(x0)) : image dans S4/V4
V4 = [(0, 1, 2, 3), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0)]
def compose(p, q): return tuple(p[q[i]] for i in range(4))
stab0 = [p for p in it.permutations(range(4)) if p[0] == 0]
cos_ = set(frozenset(compose(p, v) for v in V4) for p in stab0)
print("  Stab(x0) : ordre", len(stab0), "; classes modulo V4 atteintes :", len(cos_), "(6 ⇒ Stab(x0) → S4/V4 est un isomorphisme)")
stab01 = [p for p in it.permutations(range(4)) if p[0] == 0 and p[1] == 1]
print("  Stab(x0,x1) : ordre", len(stab01), "; caractère E sur ses éléments :", [chi[cyc(p)] for p in stab01], "⇒ droite fixe de E :", any(chi[cyc(p)] == 0 for p in stab01))
# ---- (ajout après le résultat de d1/d2, 14 h 30) : origine de la symétrie φ → −φ observée (non prévue en §1)
print("== symétrie antiunitaire Θ = W_g K (K : conjugaison dans la base géométrique ; g ∈ O \\ T)")
isreal = lambda M: np.abs(M.imag).max() < 1e-12; isimag = lambda M: np.abs(M.real).max() < 1e-12
print("  X_d réels :", all(isreal(x) for x in X), "; Y_d imaginaires :", all(isimag(y) for y in Y), "; Γ4 imaginaire :", isimag(G4),
      "; FΓ5 imaginaire :", isimag(FG), "; O_φ imaginaire :", isimag(O(0.7)))
print("  ⇒ K : termes sin (i·antisym ⊗ Y) et cos invariants ; vertical Γ4 invariant ; μFΓ5 → −μFΓ5 ; O_φ → −O_φ")
for g in R.ohgroup():
    U, ok = R.lift(g)
    if U is None: continue
    det, kT = R.classify(g)
    if det == 1 and not kT:
        Ug = U; break
Kc = lambda M: M.conj()
FGt = Ug@Kc(FG)@Ug.conj().T; G4t = Ug@Kc(G4)@Ug.conj().T
print("  g ∈ O\\T choisi : W_g K FΓ5 = %+.0f FΓ5 ; W_g K Γ4 = %+.0f Γ4" % (np.real(np.trace(FGt@FG))/8, np.real(np.trace(G4t@G4))/8))
for p in (0.3, 1.0):
    Ot = Ug@Kc(O(p))@Ug.conj().T
    a, b = np.real(np.trace(Ot@O(0)))/8, np.real(np.trace(Ot@O(np.pi/2)))/8
    print(f"  φ = {np.rad2deg(p):.1f}° → Θ O_φ Θ⁻¹ = O_ψ, ψ = {np.rad2deg(np.arctan2(b, a)):.1f}° (|O_ψ| = {np.hypot(a, b):.3f})")
