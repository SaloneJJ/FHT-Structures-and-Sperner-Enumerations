#!/usr/bin/env python3
"""Applique les criteres de P1 et P4 (05_audit/AUDIT_FABLE_PROPOSITIONS.md, ecrits le 13/11/2026 a 20 h 30) aux journaux de lancer_P1_P4.sh.
Lit la ligne « flot net par couche : [...] » ; F(N) = flot de la couche fine (dernier nombre).
Usage : python3 evaluer_P1_P4.py DOSSIER_DES_JOURNAUX
"""
import os, re, sys, math
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'P1_P4_logs')

def flot(nom):
    p = os.path.join(D, nom + '.txt')
    if not os.path.exists(p): return None
    m = re.findall(r'flot net par couche : \[([^\]]*)\]', open(p, encoding='utf-8', errors='replace').read())
    return [float(x) for x in m[-1].split()] if m else None

print('=== Bilan P1 / P4 (criteres ecrits avant calcul, 13/11/2026)')
# ---------------------------------------------------------------- P1
F3_eps0 = -0.692                     # N = 3, eps = 0 (V357/log357_eps0.0.txt) ; N = 3, eps = 0,2 : -0.106 (log357_eps0.2.txt)
a, b = flot('P1_N4_eps0.0'), flot('P1_N4_eps0.2')
print('P1 -- branche du jumeau a N = 4 (secteur gamma_F = +1, 64 pas)')
if a is None: print('  eps = 0   : pas de resultat (calcul absent ou interrompu)')
else:
    F = a[-1]
    v = ('CONFIRME' if (-1.00 <= F <= -0.75 and abs(F) > abs(F3_eps0)) else
         'NON ETABLI SUR L\'ARBRE' if abs(F) <= abs(F3_eps0) else 'ni confirme ni infirme (hors [-1,00 ; -0,75] mais |F(4)| > 0,692) : le dire')
    print(f'  eps = 0   : couches {a} ; F(4) = {F:+.3f} (F(3) = {F3_eps0:+.3f}) -> {v}')
if b is None: print('  eps = 0,2 : pas de resultat (calcul absent ou interrompu)')
else:
    F = b[-1]
    v = ('CONFIRME' if (F < 0 and abs(F) >= 0.25) else 'ECHEC' if abs(F) < 0.15 else 'NON CONCLUANT (entre 0,15 et 0,25, ou signe oppose) : le dire')
    print(f'  eps = 0,2 : couches {b} ; F(4) = {F:+.3f} (F(3) = -0.106) -> {v}')
if a is not None and b is not None:
    ok0 = (-1.00 <= a[-1] <= -0.75 and abs(a[-1]) > 0.692); ok2 = (b[-1] < 0 and abs(b[-1]) >= 0.25); fail2 = abs(b[-1]) < 0.15
    if ok0 and ok2: print('  => C14 devient [T] + [E] (N = 3, 4).')
    elif ok0 and fail2: print('  => l\'arbre ne porte le jumeau, aux tailles accessibles, que pour lambda_j <~ 0,3 : a ecrire dans A3.')
    else: print('  => lecture mixte : reporter les deux valeurs telles quelles.')
# ---------------------------------------------------------------- P4
print('P4 -- flot C3\' et pas en theta (R-unif, gamma_F = +1)')
F64, F128 = -0.943, -0.961                      # N = 3 : 64 pas (V336), 128 pas (note V338)
c = flot('P4_N3_256')
if c is None: print('  N = 3, 256 pas : pas de resultat')
else:
    F256 = c[-1]
    inc1, inc2 = F128 - F64, F256 - F128
    mono = abs(F256) > abs(F128) and abs(inc2) < abs(inc1)
    ok = -0.985 <= F256 <= -0.960 and mono
    print(f'  N = 3, 256 pas : couches {c} ; F(3) = {F256:+.3f} (64 pas : {F64}, 128 pas : {F128}) ; croissance en |F| et increments decroissants : {mono} '
          f'-> {"conforme" if ok else "non conforme"}')
    if inc1 != 0 and inc2 != 0 and inc1 * inc2 > 0 and abs(inc2) < abs(inc1):
        p = math.log2(inc1 / inc2); Finf = F256 + inc2 / (2 ** p - 1)
        okR = -1.02 <= Finf <= -0.97
        print(f'  extrapolation de Richardson (64, 128, 256) : ordre p = {p:.2f}, F(pas -> 0) = {Finf:+.3f} -> '
              + ('« poids spectral -> -1 a N fixe » [E]' if okR else 'hors [-1,02 ; -0,97] : ecrire que le flot a N = 3 sature en dessous de 1 et donner la valeur'))
    else:
        print('  extrapolation impossible (increments non monotones) : ecrire que le flot a N = 3 sature, donner la valeur')
d = flot('P4_N4_128')
if d is None: print('  N = 4, 128 pas : pas de resultat')
else:
    F = d[-1]
    print(f'  N = 4, 128 pas : couches {d} ; F(4) = {F:+.3f} (64 pas : -0.904) -> ' + ('conforme (|F(4)| > 0,904)' if abs(F) > 0.904 else 'non conforme'))
