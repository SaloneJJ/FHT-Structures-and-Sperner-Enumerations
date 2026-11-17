"""V155 : anomalies du contenu chiral du bord (piste B). Criteres et predictions : V155 (ecrits avant calcul).
Modele F : celui de V146_FHT_produit_bord.py (partie finie seulement, extraite du script).
Helicite au bord : h = gamma_F pour toutes les especes (V152 par. 1.5, V154 par. 1).
Generateurs de jauge : G = pi(X) - J_F pi(X) J_F^{-1} pour X dans u(A_F) (action Ad(u) = pi(u) J pi(u) J^{-1}).
Usage : python3 V155_FHT_anomalies_bord.py      JJ SALONE, FHT 2026"""
import itertools
from fractions import Fraction
import numpy as np

src = open('V146_FHT_produit_bord.py').read().split('# ------------------------------------------------------------------ espace Sigma')[0]
exec(src)                                    # definit dI, dF, gamma_F, GF, UF, rep, iI, DF, ...
np.set_printoptions(precision=4, suppress=True, linewidth=160)
print('=== V155 : anomalies du bord, modele F de V146 (dim', dF, ') ===')

def gen(lam, q, m):
    """generateur hermitien G sur H_F pour u = exp(i eps X), X = (lam, q, m) hermitien ; pi est R-lineaire"""
    Gp = rep(1j * lam, 1j * q, 1j * m) / 1j          # pi(iX)/i
    return Gp - UF @ Gp.conj() @ UF.T               # J X J^{-1} = U conj(X) U^T, et (i)* = -i

sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
GM = []                                             # matrices de Gell-Mann
for (i, j) in itertools.combinations(range(3), 2):
    A = np.zeros((3, 3), complex); A[i, j] = A[j, i] = 1; GM.append(A)
    A = np.zeros((3, 3), complex); A[i, j] = -1j; A[j, i] = 1j; GM.append(A)
GM.append(np.diag([1, -1, 0]).astype(complex)); GM.append(np.diag([1, 1, -2]).astype(complex) / np.sqrt(3))
Z2, Z3 = np.zeros((2, 2), complex), np.zeros((3, 3), complex)
gens = {'theta': gen(1.0, Z2, Z3), 'phi': gen(0.0, Z2, np.eye(3, dtype=complex))}
for a, s in enumerate((sx, sy, sz)): gens[f'T{a+1}'] = gen(0.0, s / 2, Z3)
for a, l in enumerate(GM): gens[f'L{a+1}'] = gen(0.0, Z2, l / 2)
names = list(gens)
h = np.diag(gamma_F).astype(complex)

# A-P1 : hermiticite, commutation avec gamma_F, forme bloc (I, II)
e_h = max(np.abs(G - G.conj().T).max() for G in gens.values())
e_g = max(np.abs(G @ GF - GF @ G).max() for G in gens.values())
e_b = max(max(np.abs(G[:dI, dI:]).max(), np.abs(G[dI:, :dI]).max()) for G in gens.values())
print(f'A-P1 : hermitiens {e_h:.1e} ; [G, gamma_F] {e_g:.1e} ; blocs hors diagonale I/II {e_b:.1e}')
# algebre de Lie : su(2) et su(3) fermees, u(1) centraux
e_c = max(np.abs(gens['T1'] @ gens['T2'] - gens['T2'] @ gens['T1'] - 1j * gens['T3']).max(),
          max(np.abs(gens[u] @ gens[n] - gens[n] @ gens[u]).max() for u in ('theta', 'phi') for n in names))
print(f'       [T1, T2] = i T3 et u(1) centraux : ecart {e_c:.1e}')

# A-P0 (partiel, lecture) : h = gamma_F uniforme par espece -> charges par etat de I
def frac(x, den=6): return Fraction(int(round(x * den)), den)
print('\nCharges sur I (diagonales), h = gamma_F :')
Qt, Qp = np.real(np.diag(gens['theta'])), np.real(np.diag(gens['phi']))
assert np.allclose(gens['theta'], np.diag(np.diag(gens['theta']))) and np.allclose(gens['phi'], np.diag(np.diag(gens['phi'])))
nomsA = ['s_a', 'd1', 'd2', 's_b']; nomsB = ['lepton', 'c1', 'c2', 'c3']
vus = set()
for a, b in itertools.product(range(4), range(4)):
    i = iI(a, b, 0); cle = (nomsA[a], 'lepton' if b == 0 else 'quark')
    if cle in vus: continue
    vus.add(cle)
    print(f'   {cle[1]:6s} {cle[0]:4s} h = {int(gamma_F[i]):+d} : Q = {Qt[i]:+.0f} theta {Qp[i]:+.0f} phi')

def T(P, *Gs, S=slice(None)):
    """Tr_S [ P . sym(G1 ... Gk) ] (symetrise)"""
    perms = list(itertools.permutations(range(len(Gs))))
    tot = 0
    for p in perms:
        M = P
        for k in p: M = M @ Gs[k]
        tot += np.trace(M[S, S])
    return np.real(tot / len(perms))

for nom, S, fac in (('I (48 champs de Weyl)', slice(0, dI), 1), ('H_F entier (96)', slice(None), 2)):
    print(f'\n=== Coefficients sur {nom} ===')
    print(f'  gravitationnel Tr[h Q] : theta {T(h, gens["theta"], S=S):+.6f} ; phi {T(h, gens["phi"], S=S):+.6f}')
    c = {}
    for k in range(4):                             # forme cubique : coefficients de theta^(3-k) phi^k
        Gs = [gens['theta']] * (3 - k) + [gens['phi']] * k
        c[k] = T(h, *Gs, S=S) * [1, 3, 3, 1][k]
    print(f'  U(1)^3 : C(theta, phi) = {c[0]:+.4f} th^3 {c[1]:+.4f} th^2 ph {c[2]:+.4f} th ph^2 {c[3]:+.4f} ph^3')
    s2 = {u: sum(T(h, gens[f'T{a}'], gens[f'T{a}'], gens[u], S=S) for a in (1, 2, 3)) / 3 for u in ('theta', 'phi')}
    s3 = {u: sum(T(h, gens[f'L{a}'], gens[f'L{a}'], gens[u], S=S) for a in range(1, 9)) / 8 for u in ('theta', 'phi')}
    print(f'  SU(2)^2 U(1) : theta {s2["theta"]:+.4f} ; phi {s2["phi"]:+.4f}   (attendu : proportionnel a (1, 3))')
    print(f'  SU(3)^2 U(1) : theta {s3["theta"]:+.4f} ; phi {s3["phi"]:+.4f}')
    d3 = max(abs(T(h, gens[f'L{a}'], gens[f'L{b}'], gens[f'L{cc}'], S=S)) for a, b, cc in itertools.combinations_with_replacement(range(1, 9), 3))
    d2 = max(abs(T(h, gens[f'T{a}'], gens[f'T{b}'], gens[f'T{cc}'], S=S)) for a, b, cc in itertools.combinations_with_replacement(range(1, 4), 3))
    mix = max(abs(T(h, gens[f'T{a}'], gens[f'T{b}'], gens[f'L{cc}'], S=S)) for a in (1, 2, 3) for b in (1, 2, 3) for cc in range(1, 9))
    mix2 = max(abs(T(h, gens[f'L{a}'], gens[f'L{b}'], gens[f'T{cc}'], S=S)) for a in range(1, 9) for b in range(1, 9) for cc in (1, 2, 3))
    mixu = max(abs(T(h, gens[u], gens[f'T{a}'], gens[f'L{b}'], S=S)) for u in ('theta', 'phi') for a in (1, 2, 3) for b in range(1, 9))
    mixu2 = max(abs(T(h, gens[u], gens[v], gens[x], S=S)) for u in ('theta', 'phi') for v in ('theta', 'phi')
                for x in [f'T{a}' for a in (1, 2, 3)] + [f'L{b}' for b in range(1, 9)])
    print(f'  SU(3)^3 max |d| : {d3:.1e} ; SU(2)^3 : {d2:.1e} ; SU(2)^2 SU(3) : {mix:.1e} ; SU(3)^2 SU(2) : {mix2:.1e} ;'
          f' U(1) SU(2) SU(3) : {mixu:.1e} ; U(1)^2 (SU(2) ou SU(3)) : {mixu2:.1e}')
    if fac == 1: ref = (c, s2, s3)
    else:
        c1, s21, s31 = ref
        e2 = max([abs(c[k] - 2 * c1[k]) for k in range(4)] + [abs(s2[u] - 2 * s21[u]) + abs(s3[u] - 2 * s31[u]) for u in ('theta', 'phi')])
        print(f'  A-P2 / E3 : ecart au facteur 2 exact (H_F = 2 x I) : {e2:.1e}')

# A-P3 : directions sans anomalie dans u(1)^2 : Q = a theta + b phi
c, s2, s3 = ref
print('\n=== A-P3 : directions de u(1)^2 sans anomalie (Q = a.theta + b.phi), sur I ===')
print(f'  forme cubique / (-6 x 3 generations) : th^3 {c[0]/-18:+.4f}, th^2 ph {c[1]/-18:+.4f}, th ph^2 {c[2]/-18:+.4f}, ph^3 {c[3]/-18:+.4f}'
      f'   (attendu : th^2 (th + 3 ph) -> 1, 3, 0, 0)')
for (a, b), lab in (((1, 0), 'theta seul'), ((0, 1), 'phi seul (nombre baryonique, au signe pres)'), ((3, -1), 'unimodulaire (theta + 3 phi = 0)'),
                    ((1, 1), 'controle (1, 1)')):
    cub = c[0] * a**3 + c[1] * a**2 * b + c[2] * a * b**2 + c[3] * b**3
    su2 = s2['theta'] * a + s2['phi'] * b; su3 = s3['theta'] * a + s3['phi'] * b
    print(f'  {lab:45s} : U(1)^3 {cub:+9.4f} ; SU(2)^2 U(1) {su2:+8.4f} ; SU(3)^2 U(1) {su3:+8.4f}')
# droite unimodulaire : charges = hypercharge ?
Y = 3 * Qt - Qp                                     # a = 3, b = -1, normalise : unites de theta (b = -a/3)
Y = Y / 3
print('  charges sur la droite unimodulaire (unites de theta) :', sorted(set(str(frac(x)) for x in Y[:dI])))
# gravitationnel le long de la droite
print(f'  gravitationnel le long de la droite : {T(h, (3 * gens["theta"] - gens["phi"]) / 3, S=slice(0, dI)):+.1e}')

# A-P4 : anomalies globales
C2 = sum(gens[f'T{a}'] @ gens[f'T{a}'] for a in (1, 2, 3))
cas = np.real(np.diag(C2))[:dI]; hI = gamma_F[:dI]
ndoub = int(round(np.sum((np.abs(cas - 0.75) < 1e-9) & (hI < 0)))) // 2
ndoub_p = int(round(np.sum((np.abs(cas - 0.75) < 1e-9) & (hI > 0)))) // 2
print('\n=== A-P4 : anomalies globales (sur I) ===')
print(f'  doublets de SU(2) : {ndoub} d\'helicite -, {ndoub_p} d\'helicite + ; par generation : {ndoub // 3} (Witten : pair requis)')
print(f'  nombre de champs de Weyl : {dI} ; modulo 16 : {dI % 16}')
