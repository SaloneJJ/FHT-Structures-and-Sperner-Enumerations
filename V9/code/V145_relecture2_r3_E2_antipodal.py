# Controle E2 (V145 §2, V144 §3.3) : (1) la "coherence diagonale" est-elle impliquee par la construction ?
# (2) variante antipodale (copies inversees, cadre absolu Sigma = (-1)^profondeur) : memes cellules ? quel test passe ?
import io, contextlib, itertools, numpy as np
from collections import defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    import V145_FHT_contact_bord as B
C, Fp, Fm, tau, sgn = B.C, B.Fp, B.Fm, B.tau, B.sgn
s_of = {c: np.array(sgn[c]) for c in C}; c_of = {tuple(v): c for c, v in s_of.items()}
anti = {c: c_of[tuple(-s_of[c])] for c in C}
def construire(nmax, regle):          # regle : +1 diagonale, -1 antipodale ; cadre absolu d'une cellule de profondeur n : regle**n
    res = {}; cells = [()]; glued = []   # (A, B, k) : A du cote -k, B du cote +k (cadre absolu)
    for n in range(1, nmax + 1):
        Sig = regle ** n; SigM = regle ** (n - 1)
        cells = [w + (c,) for w in cells for c in C]
        new = []
        for w in set(x[:-1] for x in cells):      # (ii) freres voisins (positions absolues SigM * s(c))
            for c in C:
                for k in range(3):
                    if SigM * s_of[c][k] == -1:
                        t = s_of[c].copy(); t[k] = -t[k]; new.append((w + (c,), w + (c_of[tuple(t)],), k))
        for (u, v, k) in glued:                   # (iii) heritage
            for c in C:
                if SigM * s_of[c][k] == +1:       # enfant de u du cote +k
                    t = s_of[c].copy(); t[k] = -t[k]; new.append((u + (c,), v + (c_of[tuple(t)],), k))
        glued = new
        uf = B.UF()
        for (a, b, k) in glued:
            for x in C:
                if Sig * s_of[x][k] == +1:        # coin propre de a sur sa face absolue +k
                    t = s_of[x].copy(); t[k] = -t[k]; y = c_of[tuple(t)]
                    assert (x in Fp[k] and tau[k][x] == y) or (x in Fm[k] and tau[k].get(y) == x), 'recollement hors regle (S -> S-barre)'
                    uf.u((a, x), (b, y))
        res[n] = (cells, uf)
    return res
for regle, nom in ((1, 'diagonale'), (-1, 'antipodale')):
    R = construire(4, regle)
    for n, (cells, uf) in R.items():
        cl = {(w, x): uf.f((w, x)) for w in cells for x in C}
        nsom = len(set(cl.values())); inj = all(len({cl[(w, x)] for x in C}) == 8 for w in cells)
        tests = {}
        if n >= 2:
            cellsP, ufP = R[n - 1]
            for tnom, phi in (('diag', lambda x: x), ('anti', lambda x: anti[x])):
                ok = True; grp = defaultdict(list)
                for w in cellsP:
                    for x in C: grp[ufP.f((w, x))].append((w, x))
                for mem in grp.values():
                    reps = {cl[(w + (x,), phi(x))] for (w, x) in mem}
                    if len(reps) != 1: ok = False; break
                tests[tnom] = ok
        print(f"{nom:10s} n={n} : {len(cells)} cellules, {nsom} classes de coins (attendu {(2**n+1)**3}), injectif {inj}, coherence {tests}")
# (1) preuve que la coherence diagonale est automatique : chaque identification elementaire (u,c)~(v,tau c) au niveau n-1
# provient d'une paire recollee (u,v,k) avec c dans F_k^+ ; l'heritage (iii) recolle (u+(c,), x) ~ (v+(tau c,), tau x) pour
# tout x dans F_k^+, en particulier x = c. La coherence est donc vraie par construction (une ligne), quel que soit K.
