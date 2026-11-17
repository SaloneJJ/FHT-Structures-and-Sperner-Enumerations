#!/usr/bin/env python3
"""AUDIT FABLE -- A3, section 10 (registre A3-A8) et section 4 (anomalies, C3) : controles analytiques independants.
 (1) Un niveau plat au point uniforme : H^2 = 4 sum_d sin^2(k_d/2) 1_8  =>  z(t) = 8 [e^{-2t} I_0(2t)]^3 ; n_eff(40) = 8.04, d_s(10) = 3.04.
 (2) Weyl de reseau a deux composantes : n_eff -> 2 (limite continue), correction de reseau a t = 15, 25.
 (3) Traces de jauge d'une generation (16 Weyl) : Tr T^2_SU(3) = Tr T^2_SU(2) = 2, Tr Y^2 = 10/3  => g_3^2 = g_2^2 = 5/3 g_1^2.
 (4) Anomalies du contenu chiral d'une generation avec nu_R (cubiques, mixtes, gravitationnelle, Witten).
 (5) Table des couplages a une boucle (Modele standard) aux quatre echelles de l'article.
 (6) Facteur du temps : Tr e^{-t D_4^2} = beta / sqrt(4 pi t) x Tr e^{-t D_3^2} (somme de Poisson), controle numerique.
Usage : python3 AF7_action_spectrale_controles.py     (numpy)
"""
import numpy as np
from fractions import Fraction as Fr

def i0e(x):                       # e^{-x} I_0(x) par integration numerique (pas de scipy requis)
    th = np.linspace(0, np.pi, 200001); return np.trapezoid(np.exp(x * (np.cos(th) - 1)), th) / np.pi

ok = True
print('(1) Niveau plat : z(t) = 8 [e^{-2t} I_0(2t)]^3')
z = lambda t: 8 * i0e(2 * t) ** 3
neff = lambda t: z(t) * (4 * np.pi * t) ** 1.5
ds = lambda t, h=1e-3: -2 * (np.log(z(t * (1 + h))) - np.log(z(t * (1 - h)))) / (np.log(1 + h) - np.log(1 - h))
print(f'    n_eff(40) = {neff(40):.3f} (article 8.04) ; d_s(10) = {ds(10):.3f} (article 3.04)')
ok &= abs(neff(40) - 8.04) < 6e-3 and abs(ds(10) - 3.04) < 6e-3
# verification directe sur un tore 64^3
k = 2 * np.pi * np.arange(64) / 64; s1 = np.exp(-40 * 4 * np.sin(k / 2) ** 2).sum() / 64
print(f'    tore 64^3, t = 40 : n_eff = {8 * s1**3 * (4 * np.pi * 40) ** 1.5:.3f}')

print('(2) Weyl de reseau (2 composantes, E^2 = 4 sum sin^2(k/2) au voisinage du cone ; ici modele de reference E^2 = sum sin^2 k, 8 cones dont 1 retenu)')
for t in (15, 25, 40):
    print(f'    t = {t} : n_eff d\'une espece a 2 composantes avec dispersion 4 sin^2(k/2) : {2 * i0e(2 * t) ** 3 * (4 * np.pi * t) ** 1.5:.3f}  (limite continue : 2 ; correction ~ 3/(16 t) = {3 / (16 * t):.3f})')

print('(3) Traces de jauge d\'une generation (Weyl gauches : Q, u^c, d^c, L, e^c, nu^c)')
# (dim SU3, dim SU2, Y) avec Q = T3 + Y
gen = [(3, 2, Fr(1, 6)), (3, 1, Fr(-2, 3)), (3, 1, Fr(1, 3)), (1, 2, Fr(-1, 2)), (1, 1, Fr(1)), (1, 1, Fr(0))]
T3 = sum(Fr(1, 2) * d2 for d3, d2, y in gen if d3 == 3); T2 = sum(Fr(1, 2) * d3 for d3, d2, y in gen if d2 == 2); Y2 = sum(d3 * d2 * y * y for d3, d2, y in gen)
print(f'    Tr T^2_SU(3) = {T3}, Tr T^2_SU(2) = {T2}, Tr Y^2 = {Y2}  (article : 2, 2, 10/3) ; Tr Y^2 / Tr T^2 = {Y2 / T2} = 5/3')
ok &= (T3, T2, Y2) == (2, 2, Fr(10, 3)) and sum(d3 * d2 for d3, d2, y in gen) == 16

print('(4) Anomalies')
A = {'U(1)^3': sum(d3 * d2 * y**3 for d3, d2, y in gen), 'grav-U(1)': sum(d3 * d2 * y for d3, d2, y in gen),
     'SU(2)^2 U(1)': sum(d3 * y for d3, d2, y in gen if d2 == 2), 'SU(3)^2 U(1)': sum(d2 * y for d3, d2, y in gen if d3 == 3),
     'SU(3)^3': sum(d2 * (1 if y == Fr(1, 6) else -1) for d3, d2, y in gen if d3 == 3)}
ndoublets = sum(d3 for d3, d2, y in gen if d2 == 2)
print('   ', {k_: str(v) for k_, v in A.items()}, f'; doublets de SU(2) par generation : {ndoublets} (pair : pas d\'anomalie de Witten) ; 16 Weyl par generation (Z/16)')
ok &= all(v == 0 for v in A.values()) and ndoublets % 2 == 0

print('(5) Couplages a une boucle, entrees a M_Z : alpha^-1 = 127.95, sin^2 theta_W = 0.2312, alpha_s = 0.118 ; b = (41/10, -19/6, -7)')
MZ = 91.1876; a1, a2, a3 = (3 / 5) * (1 - 0.2312) * 127.95, 0.2312 * 127.95, 1 / 0.118; b = (41 / 10, -19 / 6, -7)
labels = []
table = {1e16: (37.9, 45.9, 44.5, -17, -3), 1e17: (36.4, 47.0, 47.1, -23, 0), 2.4e18: (34.3, 48.6, 50.6, -29, 4), 1.22e19: (33.3, 49.5, 52.4, -33, 6)}
for L, ref in table.items():
    t = np.log(L / MZ) / (2 * np.pi); v = [a - bb * t for a, bb in zip((a1, a2, a3), b)]
    r1, r3 = (v[1] / v[0] - 1) * 100, (v[1] / v[2] - 1) * 100            # alpha_1/alpha_2 - 1 = alpha_2^-1/alpha_1^-1 - 1 (definition de l'en-tete de la table)
    i1, i3 = (v[0] / v[1] - 1) * 100, (v[2] / v[1] - 1) * 100            # alpha_1^-1/alpha_2^-1 - 1 (rapport des inverses)
    good = all(abs(x - y) < 0.06 for x, y in zip(v, ref[:3])); ok &= good
    lab_ok = abs(r1 - ref[3]) < 0.7 and abs(r3 - ref[4]) < 0.7; inv_ok = abs(i1 - ref[3]) < 0.7 and abs(i3 - ref[4]) < 0.7; labels.append((lab_ok, inv_ok))
    print(f'    Lambda = {L:.3g} GeV : alpha_i^-1 = {v[0]:.1f}, {v[1]:.1f}, {v[2]:.1f} (article {ref[0]}, {ref[1]}, {ref[2]})  {"ok" if good else "ECART"} ; '
          f'alpha_1/alpha_2 - 1 = {r1:+.0f} %, alpha_3/alpha_2 - 1 = {r3:+.0f} % ; rapports des inverses : {i1:+.0f} %, {i3:+.0f} % ; article : {ref[3]:+d} %, {ref[4]:+d} %')
print('    => les deux dernieres colonnes de la table de l\'article sont', 'alpha_1/alpha_2 - 1 et alpha_3/alpha_2 - 1 (en-tete exact)' if all(a for a, _ in labels) else
      ('alpha_1^-1/alpha_2^-1 - 1 et alpha_3^-1/alpha_2^-1 - 1 : EN-TETE A CORRIGER (les valeurs sont justes pour les rapports des inverses)' if all(b_ for _, b_ in labels) else 'ni l\'un ni l\'autre : ECART'))
ok &= all(b_ for _, b_ in labels) or all(a for a, _ in labels)
print('(6) Facteur du temps (frequences de Matsubara fermioniques sur S^1_beta)')
beta, t = 30.0, 2.0; n = np.arange(-4000, 4000); s = np.exp(-t * ((2 * n + 1) * np.pi / beta) ** 2).sum()
print(f'    sum_n e^(-t w_n^2) = {s:.8f} ; beta / sqrt(4 pi t) = {beta / np.sqrt(4 * np.pi * t):.8f}'); ok &= abs(s - beta / np.sqrt(4 * np.pi * t)) < 1e-8
print('VERDICT :', 'conforme' if ok else 'ECART')
