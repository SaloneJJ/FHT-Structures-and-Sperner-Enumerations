"""R3a, version corrigee (14/11/2026, soir).  Correction : un caractere reel d'indicateur de Frobenius-Schur -1 n'est pas le
caractere d'un QG-module ; l'irreductible rationnel correspondant est 2*chi (indice de Schur 2).  La version 1 testait chi
(cas de Q8 : 'False' etait un artefact).  On compare aussi le reseau sans signes (permutations seules, anneau de Burnside)
et le reseau avec signes Lambda_+-(G).  Groupes ajoutes : SL(2,3), Q8 x C3, Dic12 (ou la condition pourrait mordre)."""
from twisted_alg import *
from burnside_signe import subgroups, homs_pm, induced, lattice_basis, in_lattice
import numpy as np
def fs(G, chi_of):
    return round(sum(chi_of(mul(g, g)) for g in G).real / len(G))
def analyse2(G, name):
    C, cl, idx, reps, X = char_table(G)
    k = len(C); e = tuple(range(len(G[0]))); e0 = [l for l in range(k) if reps[l] == e][0]
    cls_of = {}
    for l, c in enumerate(C):
        for g in c: cls_of[g] = l
    gens_pm, gens_p = [], []
    for H in subgroups(G):
        for eps in homs_pm(H):
            v = [int(x) for x in induced(G, sorted(H), eps, reps)]
            gens_pm.append(v)
            if all(s == 1 for s in eps.values()): gens_p.append(v)
    Bpm = lattice_basis(gens_pm); Bp = lattice_basis(gens_p)
    out = []; used = set()
    for a in range(k):
        if a in used: continue
        chi = X[a]
        if np.allclose(chi.imag, 0, atol=1e-6):
            v = chi.real; used.add(a)
            nu = fs(G, lambda g: chi[cls_of[g]]); m = 2 if nu == -1 else 1
        else:
            b = [bb for bb in range(k) if np.allclose(X[bb], np.conj(chi), atol=1e-6)][0]
            v = (chi + X[b]).real; used |= {a, b}; nu = 0; m = 1
        w = np.round(m * v).astype(int)
        out.append((int(w[e0]), nu, m, in_lattice(Bp, w), in_lattice(Bpm, w)))
    print(f'{name} : |G| = {len(G)}, classes {k}, rang Burnside {len(Bp)}, rang Lambda+- {len(Bpm)}')
    for d, nu, m, ip, ipm in out:
        print(f'   dim_Q {d:2d}  (FS {nu:+d}, Schur {m})   permutations seules : {ip!s:5}   avec signes Lambda+- : {ipm}')
    return out
if __name__ == '__main__':
    vecs = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
    def mp(M): return tuple(vecs.index(((M[0][0]*a + M[0][1]*b) % 3, (M[1][0]*a + M[1][1]*b) % 3)) for a, b in vecs)
    analyse2(closure([(1, 0, 2, 3), (1, 2, 3, 0)]), 'S4')
    analyse2(closure([(1, 0, 2, 3, 4), (1, 2, 3, 4, 0)]), 'S5')
    analyse2(closure([mp([[1, 1], [0, 1]]), mp([[0, 1], [2, 0]]), mp([[2, 0], [0, 1]])]), 'GL(2,3)')
    q8 = [mp([[0, 1], [2, 0]]), mp([[1, 1], [1, 2]])]
    analyse2(closure(q8), 'Q8')
    analyse2(closure([mp([[1, 1], [0, 1]]), mp([[0, 1], [2, 0]])]), 'SL(2,3)')
    ext = lambda p: tuple(p) + (8, 9, 10); c3 = tuple(range(8)) + (9, 10, 8)
    analyse2(closure([ext(q8[0]), ext(q8[1]), c3]), 'Q8 x C3')
    # Dic12 = <a, x | a^6 = 1, x^2 = a^3, x a x^-1 = a^-1>, action reguliere
    els = [(i, j) for j in range(2) for i in range(6)]
    def m_(p, q):
        (i, j), (k, l) = p, q
        if j == 0: return ((i + k) % 6, l)
        r = (i - k) % 6
        return ((r + 3) % 6, 0) if l == 1 else (r, 1)
    regl = lambda g: tuple(els.index(m_(g, h)) for h in els)
    analyse2(closure([regl((1, 0)), regl((0, 1))]), 'Dic12 = C3 x| C4')
