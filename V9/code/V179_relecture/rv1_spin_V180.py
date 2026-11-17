# Relecture V180 : classe projective des relevés (octet, K26, tétraèdre), sauts γ_i.
import numpy as np, itertools as it
from rv_common import *

def surface(kind):
    if kind == "cube":
        pts = [np.array(p, float) for p in it.product([-1, 1], repeat=3)]
        edges = [(a, b) for a, b in it.combinations(range(8), 2) if np.isclose(np.linalg.norm(pts[a] - pts[b]), 2)]
        faces = [orient_out(pts, [i for i in range(8) if pts[i][ax] == s]) for ax in range(3) for s in (-1, 1)]
    elif kind == "K26":
        pts = [np.array(p, float) for p in it.product([-1, 0, 1], repeat=3) if any(p)]
        edges = [(a, b) for a, b in it.combinations(range(26), 2) if np.isclose(np.linalg.norm(pts[a] - pts[b]), 1)]
        faces = []
        # carrés unité contenus dans les faces du cube [-1,1]^3 : 4 sommets de la coquille formant un carré de côté 1
        for quad in it.combinations(range(26), 4):
            P = np.array([pts[i] for i in quad])
            for ax in range(3):
                if np.all(P[:, ax] == P[0, ax]) and abs(P[0, ax]) == 1:
                    o = [k for k in range(3) if k != ax]
                    if np.ptp(P[:, o[0]]) == 1 and np.ptp(P[:, o[1]]) == 1 and len({tuple(p) for p in P[:, o]}) == 4:
                        faces.append(orient_out(pts, list(quad)))
    elif kind == "tet":
        pts = [np.array(p, float) for p in it.product([-1, 1], repeat=3) if np.prod(p) > 0]
        edges = list(it.combinations(range(4), 2))
        faces = [orient_out(pts, list(f)) for f in it.combinations(range(4), 3)]
    return pts, edges, faces

def ham(kind, phi_per_face, total_m):
    pts, edges, faces = surface(kind)
    nf = len(faces); flux = np.full(nf, phi_per_face); flux[-1] -= 2 * np.pi * total_m   # somme nulle exacte
    A = peierls(edges, faces, flux, len(edges))
    return ham_from_A(len(pts), edges, A), pts, edges, faces

def klass(H, pts, group=None):
    Ux, _ = lift_diag(H, pts, HX); Uy, _ = lift_diag(H, pts, HY)
    lam, sc = commutator_class(Ux, Uy)
    return lam, sc

print("=== Octet : flux uniforme sortant 2πm/6, m = 0..5 ===")
for m in range(6):
    H, pts, E, F = ham("cube", 2 * np.pi * m / 6, m)
    lam, sc = klass(H, pts)
    Us = [lift_diag(H, pts, R)[0] for R in ROTS]
    print(" m=%d : [Ux,Uy] = %+.4f%+.4fi scalaire=%s ; relevés 24/24 : %s ; niveaux %s" % (m, lam.real, lam.imag, sc, all(U is not None for U in Us), level_structure(H, Us)))

# Contrôle de la remarque V182 §0 : ℤ/2 avec un relèvement NON uniforme (+π sur 5 faces, −π sur une) : même H à une jauge près
H3, pts, E, F = ham("cube", np.pi, 3)
flux = np.full(6, np.pi); flux[0] = -np.pi; flux[-1] -= 2*np.pi*2      # somme = 4π - 4π = 0 ; « m = 2 »
H2b = ham_from_A(8, E, peierls(E, F, flux, len(E)))
print(" ℤ/2 relevé non uniforme (« m = 2 ») : [Ux,Uy] = %+.3f ; spectre identique à m=3 : %s"
      % (klass(H2b, pts)[0].real, np.allclose(np.linalg.eigvalsh(H2b), np.linalg.eigvalsh(H3))))

print("=== Q3 : sauts par direction (flux π) ===")
G = []
for ax in range(3):
    g = np.zeros((8, 8), complex)
    for a, b in E:
        if abs(pts[a][ax] - pts[b][ax]) > 1: g[a, b] = H3[a, b]; g[b, a] = H3[b, a]
    G.append(g)
print(" γ_i² = 1 :", [bool(np.allclose(g @ g, np.eye(8))) for g in G],
      "; anticommutation :", [bool(np.allclose(G[i] @ G[j] + G[j] @ G[i], 0)) for i, j in [(0, 1), (0, 2), (1, 2)]],
      "; Σγ = H :", bool(np.allclose(sum(G), H3)))
ok = True; signs = []
for R in ROTS:
    U, _ = lift_diag(H3, pts, R)
    for i in range(3):
        j = int(np.argmax(np.abs(R[:, i])))           # R e_i = ± e_j
        X = U @ G[i] @ U.conj().T
        s = np.trace(G[j].conj().T @ X) / 8
        ok &= bool(np.isclose(abs(s), 1) and np.allclose(X, s * G[j], atol=1e-9)); signs.append(np.round(s.real, 6))
print(" covariance U γ_i U† = s γ_{R(i)}, 24 rotations :", ok, "; signes s rencontrés :", sorted(set(signs)))

print("=== Q4 : coquille K26 ===")
pts26, E26, F26 = surface("K26")
print(" sites %d, arêtes %d, carrés %d, Euler %d" % (len(pts26), len(E26), len(F26), len(pts26) - len(E26) + len(F26)))
for lab, phi, m in [("ℤ/2 π (m=12)", np.pi, 12), ("U(1) m=1", 2*np.pi/24, 1), ("U(1) m=2", 2*np.pi*2/24, 2), ("U(1) m=3", 2*np.pi*3/24, 3), ("U(1) m=5", 2*np.pi*5/24, 5)]:
    H, p, _, _ = ham("K26", phi, m); lam, sc = klass(H, p)
    extra = ""
    if m == 1:
        Us = [lift_diag(H, p, R)[0] for R in ROTS]; extra = " niveaux bas : %s" % level_structure(H, Us)[:3]
    print(" %-14s [Ux,Uy] = %+.4f (scalaire %s)%s" % (lab, lam.real, sc, extra))

print("=== Q5 : tétraèdre R, groupe T ===")
pts4 = surface("tet")[0]
T = [R for R in ROTS if site_perm(pts4, R) is not None]
for lab, phi, m in [("π/2 (m=1)", np.pi/2, 1), ("0", 0.0, 0), ("π (m=2)", np.pi, 2), ("3π/2 (m=3)", 3*np.pi/2, 3)]:
    H, p, _, _ = ham("tet", phi, m); lam, sc = klass(H, p)
    Us = [lift_diag(H, p, R)[0] for R in T]
    print(" %-11s |T|=%d [Ux,Uy] = %+.4f ; niveaux %s" % (lab, len(T), lam.real, level_structure(H, Us)))
