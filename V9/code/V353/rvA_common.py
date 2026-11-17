# V331 relecture -- outils indépendants (aucune fonction des scripts relus n'est importée).
# BBH écrit sous la forme séparable de l'article VII (Thm bbh) : H = h⊗1⊗1 + Z⊗h⊗1 + Z⊗Z⊗h.
# Base : s = (s1,s2,s3) ∈ {0,1}^3, indice 4*s1+2*s2+s3 ; position du site dans la cellule : s - 1/2.
import numpy as np, itertools as it

TF, TC = 0.4, 1.0
sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
SITES = list(it.product([0, 1], repeat=3))
POS = np.array([np.array(s) - 0.5 for s in SITES])

def h1(k, tf=TF, tc=TC):
    return np.array([[0, tf + tc*np.exp(-1j*k)], [tf + tc*np.exp(1j*k), 0]])

def kron3(a, b, c): return np.kron(np.kron(a, b), c)

def H(k, tf=TF, tc=TC):
    k = np.asarray(k, float)
    return kron3(h1(k[0], tf, tc), I2, I2) + kron3(sz, h1(k[1], tf, tc), I2) + kron3(sz, sz, h1(k[2], tf, tc))

KPT = np.array([np.pi]*3)
S = np.diag([(-1)**sum(s) for s in SITES]).astype(complex)       # signe de sous-réseau

def comm(A, B): return A@B - B@A
def acomm(A, B): return A@B + B@A

def nullspace(M, tol=1e-9):
    u, s, vh = np.linalg.svd(M)
    r = int(np.sum(s > tol*max(1., s[0])))
    return vh[r:].conj()            # conjugaison : M v = 0 pour v = vh[i].conj()

def commutant(ops, n=8):
    M = np.vstack([np.kron(A, np.eye(n)) - np.kron(np.eye(n), A.T) for A in ops])   # vec ligne : A X - X A
    return [v.reshape(n, n) for v in nullspace(M)]

# --- Clifford au point K : composantes paires (cos) et impaires (sin) en q
def clifford():
    G = H(KPT)                         # = m Γ5 (seul terme en K)
    G5 = G/np.sqrt(np.real(np.trace(G@G))/8)
    eps = 1e-4; al = []; cs = []
    for d in range(3):
        e = np.zeros(3); e[d] = eps
        al.append((H(KPT+e) - H(KPT-e))/(2*eps))
        cs.append((H(KPT+e) + H(KPT-e) - 2*H(KPT))/eps**2)   # = t_c Γ'_{c,d} (dérivée seconde)
    al = [a/np.sqrt(np.real(np.trace(a@a))/8) for a in al]
    cs = [c/np.sqrt(np.real(np.trace(c@c))/8) for c in cs]
    # Γ5 doit être ∝ somme des Γ'_c ; A, B : composantes orthogonales (rep E de S3)
    A = (cs[0] - cs[1])/np.sqrt(2); B = (cs[0] + cs[1] - 2*cs[2])/np.sqrt(6)
    A /= np.sqrt(np.real(np.trace(A@A))/8); B /= np.sqrt(np.real(np.trace(B@B))/8)
    G4 = G5@al[0]@al[1]@al[2]
    return dict(G5=G5, al=al, cs=cs, A=A, B=B, G4=G4)

# --- relevés des éléments de O_h : U = D P avec U H(k) U† = H(g k) pour tout k (D diagonale, k-indépendante)
def ohgroup():
    out = []
    for p in it.permutations(range(3)):
        P = np.eye(3)[list(p)]
        for sg in it.product([1, -1], repeat=3):
            out.append(np.diag(sg) @ P)
    return out

def site_perm(g):
    img = POS @ g.T
    return [int(np.argmin(np.linalg.norm(POS - q, axis=1))) for q in img]

def lift(g, nk=6, seed=0):
    perm = site_perm(g); n = 8
    P = np.zeros((n, n)); P[perm, np.arange(n)] = 1
    rng = np.random.default_rng(seed); rows = []
    for _ in range(nk):
        k = rng.uniform(-np.pi, np.pi, 3)
        Hk, Hgk = H(k), H(g@k)
        M = np.zeros((n*n, n), complex)
        for j in range(n):
            E = np.zeros((n, n), complex); E[j, j] = 1
            M[:, j] = (Hgk@E@P - E@P@Hk).ravel()
        rows.append(M)
    ns = nullspace(np.vstack(rows))
    if len(ns) != 1: return None, len(ns)
    d = ns[0]/ns[0][np.argmax(np.abs(ns[0]))]
    if not np.allclose(np.abs(d), 1, atol=1e-8): return None, -1
    U = np.diag(d)@P
    for k in rng.uniform(-np.pi, np.pi, (4, 3)):
        assert np.allclose(U@H(k)@U.conj().T, H(g@k), atol=1e-9)
    return U, 1

def classify(g):
    """T : rotations préservant chaque tétraèdre ; etc."""
    det = round(np.linalg.det(g)); perm = site_perm(g)
    even = {i for i, s in enumerate(SITES) if sum(s) % 2 == 0}
    keeps = {perm[i] for i in even} == even
    return det, keeps
