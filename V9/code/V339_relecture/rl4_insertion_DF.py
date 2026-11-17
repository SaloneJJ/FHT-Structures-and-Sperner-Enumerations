"""V339 relecture, partie 4 : insertion de D_F sur le survivant (point 5, non testé dans V334-V338).
Modèle réduit de Kaplan (chaîne verticale ouverte de Nz sites, Wilson r = 1 habillé par γ_F, Γ4 nu), octet natif au cône (t_f = t_c, q latéral fixé),
fibre interne ℂ² : γ_F = σz, D_F = σx. H = H_impair ⊗ 1 + H_pair ⊗ γ_F + μ FΓ5 ⊗ 1 + ε O ⊗ D_F, avec O = Γ5 (insertion de V317/V327) ou d'autres O.
Question : les modes nuls de bord (un Weyl par secteur γ_F, goût f = −γ_F) prennent-ils une masse d'ordre ε ?"""
import numpy as np, sys, itertools as it
sys.path.insert(0, '.')
from rl_common import *
H0, M = octet_mats(1, 1)
Xd = [M[d] + M[d].T for d in range(3)]; Yd = [1j * (M[d] - M[d].T) for d in range(3)]
G5 = sum(Xd) / np.sqrt(3); G4 = G5 @ Yd[0] @ Yd[1] @ Yd[2]
A = (Xd[0] - Xd[1]) / np.sqrt(2); B = (Xd[0] + Xd[1] - 2 * Xd[2]) / np.sqrt(6); F = 1j * A @ B; FG = F @ G5
gz = np.diag([1., -1.]); DF = np.array([[0, 1.], [1, 0]])
def Hchain(q, mu, eps, O, Nz=24, mtot=0.0):
    K = np.full(3, np.pi); Hq = Hk(K + q); Hm = Hk(K - q)
    ev, od = (Hq + Hm) / 2, (Hq - Hm) / 2
    T = np.diag(np.ones(Nz - 1), -1)
    odd = np.kron(np.eye(Nz), od) + np.kron(-0.5j * (T - T.T), G4)
    even = np.kron(np.eye(Nz), ev + mtot * G5) + np.kron(0.5 * (2 * np.eye(Nz) - T - T.T), G5)
    return (np.kron(odd, np.eye(2)) + np.kron(even, gz) + mu * np.kron(np.kron(np.eye(Nz), FG), np.eye(2))
            + eps * np.kron(np.kron(np.eye(Nz), O), DF))
def low(Hm, k=8): return np.sort(np.abs(np.linalg.eigvalsh(Hm)))[:k]
q0 = np.zeros(3); q1 = 0.05 * np.array([1, .6, .3])
print('Contrôle sans D_F : μ = 0, m_tot = −1 (deux copies) : 8 plus petits |E| à q = 0 :', np.round(low(Hchain(q0, 0, 0, G5, mtot=-1)), 5))
print('                    μ = 1 (un survivant par secteur) :', np.round(low(Hchain(q0, 1, 0, G5)), 5))
for eps in (0.05, 0.1, 0.2):
    print(f'ε = {eps} :  O = Γ5 : μ = 0, m_tot = −1 → {np.round(low(Hchain(q0, 0, eps, G5, mtot=-1), 4), 5)} ; μ = 1 → {np.round(low(Hchain(q0, 1, eps, G5), 4), 5)} (q = 0) ; '
          f'{np.round(low(Hchain(q1, 1, eps, G5), 4), 5)} (|q| = 0,05)')
# recherche des O (monômes de Clifford sur l'octet) : {O, α_d} = {O, Γ4} = 0, [O, Γ5] = 0, {O, F} = 0
gens = {'G5': G5, 'A': A, 'B': B, 'Y1': Yd[0], 'Y2': Yd[1], 'Y3': Yd[2]}
print('\nMonômes O (produits des 6 générateurs) qui anticommutent avec α_d et Γ4, commutent avec Γ5, et anticommutent avec F :')
found = []
for r in range(1, 7):
    for cmb in it.combinations(gens, r):
        O = np.eye(8, dtype=complex)
        for g in cmb: O = O @ gens[g]
        if not np.allclose(O, O.conj().T): O = 1j * O
        ok = all(np.allclose(O @ Z, -Z @ O) for Z in Yd + [G4]) and np.allclose(O @ G5, G5 @ O) and np.allclose(O @ F, -F @ O)
        if ok: found.append(('·'.join(cmb), O))
for nm, O in found:
    print(f'  O = {nm} : ε = 0,1, μ = 1, q = 0 → {np.round(low(Hchain(q0, 1, 0.1, O), 4), 5)}')
print('Lecture : A et B portent la représentation E de O (et de T, où E se scinde en deux représentations complexes conjuguées) :'
      '\n  aucun O de ce type n\'est invariant sous T ; un couplage de Yukawa qui donne une masse au survivant doit briser T (ou passer par un autre canal).')

U = G5 @ A @ B     # J_lat = U K (convention « cos », V331)
for nm, O in [('Γ5', G5)] + found:
    s_ = '+' if np.allclose(U @ O.conj() @ np.linalg.inv(U), O) else ('−' if np.allclose(U @ O.conj() @ np.linalg.inv(U), -O) else '?')
    print(f'  J_lat O J_lat⁻¹ = {s_}O pour O = {nm}')
print('=> le O qui donne la masse au survivant a la MÊME parité sous J_lat que Γ5 si les signes coïncident (compatibilité KO-1 à vérifier dans ce cas).')
# contrôle : terme de goût HABILLÉ (μ FΓ5 ⊗ γ_F, exclu par V331 pour la structure réelle) avec O = Γ5
def Hchain_dressed(q, mu, eps, O, Nz=24):
    Hm0 = Hchain(q, 0.0, eps, O, Nz)
    return Hm0 + mu * np.kron(np.kron(np.eye(Nz), FG), gz)
print('\nContrôle : μ FΓ5 habillé par γ_F, O = Γ5, ε = 0,1, q = 0 →', np.round(low(Hchain_dressed(q0, 1, 0.1, G5), 4), 5),
      '(les survivants ont alors le même goût dans les deux secteurs et prennent la masse ε ; mais ce terme détruit la structure réelle, V331 §4)')
