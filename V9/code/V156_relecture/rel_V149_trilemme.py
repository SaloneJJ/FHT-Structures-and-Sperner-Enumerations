"""Relecture de V149 §6 : (a) KO-1 <=> s impaire et w paire pour gamma_F (x) (w + sigma.s) ; s(k*) = 0 aux 8 TRIM (N pair).
(b) contre-exemple au pas (i) : terme 1_F (x) u(k).sigma avec u PAIRE : compatible KO-1, invariant de jauge, local, invariant par translation,
ne cree pas de surface nodale et reduit les noeuds de 8 a 2 par espece (doublage minimal, chiralites opposees : NN respecte).
Utilise le F reconstruit (rel_F.py), J = (J_F gamma_F) (x) (i sigma_y o conj), k -> -k."""
import itertools, numpy as np
from rel_F import *
rng = np.random.default_rng(7); DF = build_DF(rng)
sx, sy, sz = pauli; I2 = np.eye(2)
U = np.kron(UF @ GF, 1j * sy)                       # partie unitaire de J (agit sur F (x) spin), J = U o conj o (k -> -k)
def D(k, DFm, u_on=True, r=1.0):
    s = np.sin(k); u = np.array([0, 0, r * np.sum(1 - np.cos(k))]) if u_on else np.zeros(3)
    return np.kron(DFm, I2) + np.kron(GF, s[0] * sx + s[1] * sy + s[2] * sz) + np.kron(np.eye(dF), u[0] * sx + u[1] * sy + u[2] * sz)
Nl = 8; ks = [2 * np.pi * np.array(n) / Nl for n in itertools.product(range(Nl), repeat=3)]
err = max(np.abs(U @ D(k, DF).conj() @ U.conj().T + D(-k, DF)).max() for k in ks[::7])
print(f'KO-1 avec le terme 1 (x) u.sigma (u paire, u_z = r sum(1 - cos k)) : max |J D(k) J^-1 + D(-k)| = {err:.1e} ; J^2 = +1 : {np.allclose(U @ U.conj(), np.eye(2*dF))}')
bad = max(np.abs(U @ np.kron(np.eye(dF), sz).conj() @ U.conj().T - np.kron(np.eye(dF), sz)).max(), 0)
print(f'(controle : le meme terme avec u IMPAIRE violerait KO-1 ; J (1 (x) sigma_z) J^-1 = -(1 (x) sigma_z) : {np.allclose(U @ np.kron(np.eye(dF), sz).conj() @ U.conj().T, -np.kron(np.eye(dF), sz))})')
basis = real_basis(); Jb = [Jconj(UF, b) for b in basis]
k0 = ks[37]; o1 = max(np.abs((D(k0, DF) @ np.kron(a, I2) - np.kron(a, I2) @ D(k0, DF)) @ np.kron(b, I2) - np.kron(b, I2) @ (D(k0, DF) @ np.kron(a, I2) - np.kron(a, I2) @ D(k0, DF))).max() for a in basis for b in Jb)
print(f'ordre un pour A_F : {o1:.1e} ; le terme commute avec pi(A_F) et J pi(A_F) J^-1 (identite sur F)')
# noeuds (D_F = 0) : nombre de k du reseau N = 8 ou D(k) a une valeur propre nulle, par signe de gamma_F
for u_on in (False, True):
    nodes = {+1: [], -1: []}
    for k in ks:
        s = np.sin(k); u = np.array([0, 0, np.sum(1 - np.cos(k))]) if u_on else np.zeros(3)
        for g in (+1, -1):
            if np.linalg.norm(g * s + u) < 1e-12: nodes[g].append(tuple(np.round(k / np.pi, 3)))
    print(f'{"avec" if u_on else "sans"} terme u : noeuds (k/pi) gamma_F=+1 : {nodes[1]} ; gamma_F=-1 : {nodes[-1]}')
# chiralite des noeuds : signe du jacobien de f = g s + u
def chir(k, g):
    h = 1e-6; J = np.zeros((3, 3))
    f = lambda q: g * np.sin(q) + np.array([0, 0, np.sum(1 - np.cos(q))])
    for j in range(3):
        e = np.zeros(3); e[j] = h; J[:, j] = (f(k + e) - f(k - e)) / (2 * h)
    return int(np.sign(np.linalg.det(J)))
print('chiralites (gamma_F = +1) : k = 0 ->', chir(np.zeros(3), 1), '; k = (0,0,-pi/2) ->', chir(np.array([0, 0, -np.pi / 2]), 1))
# surface nodale ? minimum de |f| hors noeuds sur un reseau fin
Nf = 64; mins = []
for n in itertools.product(range(Nf), repeat=3):
    k = 2 * np.pi * np.array(n) / Nf; f = np.sin(k) + np.array([0, 0, np.sum(1 - np.cos(k))]); mins.append(np.linalg.norm(f))
mins = np.array(mins); print(f'reseau fin {Nf}^3 : nombre de k avec |f| < 0.05 : {np.sum(mins < 0.05)} (deux cones ponctuels => O(1), pas une surface)')

# (c) meme terme avec la derivee SLAC (mode de Nyquist -> 0, KO-1 exact) : un seul noeud par espece sur N pair
def s_slac(k):
    q = (k + np.pi) % (2 * np.pi) - np.pi
    return np.where(np.isclose(np.abs(q), np.pi), 0.0, q)
for Nl in (8, 16):
    ks = [2 * np.pi * np.array(n) / Nl for n in itertools.product(range(Nl), repeat=3)]
    nodes = {+1: 0, -1: 0}; minf = {+1: 9, -1: 9}
    for k in ks:
        u = np.array([0, 0, np.sum(1 - np.cos(k))])
        for g in (1, -1):
            f = g * s_slac(k) + u; nf = np.linalg.norm(f)
            if nf < 1e-12: nodes[g] += 1
            elif nf < minf[g]: minf[g] = nf
    # KO-1 : s_slac(-k) = -s_slac(k) sur le reseau ?
    ko = max(np.abs(s_slac(-k) + s_slac(k)).max() for k in ks)
    print(f'SLAC + u, N = {Nl} : noeuds gamma_F=+1 : {nodes[1]}, gamma_F=-1 : {nodes[-1]} ; plus petit |f| hors noeud : {minf[1]:.3f} ; ecart KO-1 de s : {ko:.1e}')
