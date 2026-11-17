# Controle E4 (V145 §4) : generateurs = automorphismes ? image du stabilisateur du cube (B) sur les 8 coins :
# groupe ponctuel (ordre, contient-il des elements de det -1 = chiralite de la tuile ?), invariants quadratiques (metrique).
import json, itertools, numpy as np
from collections import deque
from V87_FHT_stella import rules_stella, faces_of
from V86_FHT_k26 import build_K
import fht_tools
D = json.load(open('V145_etoiles_facettes.json'))
FB = [frozenset(f) for f in faces_of(rules_stella())]
Fp = [FB[2*k] for k in range(3)]
C = sorted(frozenset().union(*FB))
sgn = {c: np.array([1 if c in Fp[k] else -1 for k in range(3)]) for c in C}
def objets():
    K, I, II = build_K(); V = sorted(set().union(*K), key=str)
    S4g = [(1, 0, 2, 3), (1, 2, 3, 0)]; idn = (0, 1, 2, 3); gens = []
    for gL, gR, gg in [(g, idn, idn) for g in S4g] + [(idn, g, idn) for g in S4g] + [(idn, idn, g) for g in S4g]:
        gens.append({v: fht_tools.act3(gL, gR, gg, v) for v in V})
    yield 'K26', [sorted(f) for f in K], gens
    for nom, e in D['etoiles'].items(): yield nom, e['facets'], e['gens']
for nom, F, gens in objets():
    Fs = set(frozenset(f) for f in F); V = sorted(set().union(*Fs), key=str)
    ok = all(set(frozenset(g[v] for v in f) for f in Fs) == Fs for g in gens)
    # groupe engendre, stabilisateur de la famille des faces (B), image sur les coins
    ix = {v: i for i, v in enumerate(V)}; G0 = [tuple(ix[g[v]] for v in V) for g in gens]
    e = tuple(range(len(V))); seen = {e}; q = deque([e])
    while q:
        p = q.popleft()
        for g in G0:
            r = tuple(g[i] for i in p)
            if r not in seen: seen.add(r); q.append(r)
    FBi = set(frozenset(ix[v] for v in f) for f in FB)
    stab = [p for p in seen if set(frozenset(p[i] for i in f) for f in FBi) == FBi]
    imgs = set()
    for p in stab: imgs.add(tuple(V[p[ix[c]]] for c in C))
    dets = []; mats = []
    for im in imgs:
        # matrice signee M avec M s(c) = s(im(c))
        S = np.array([sgn[c] for c in C]).T; T = np.array([sgn[x] for x in im]).T
        M = np.round(T @ np.linalg.pinv(S)).astype(int)
        assert np.allclose(M @ S, T), 'pas une isometrie du cube'
        dets.append(int(round(np.linalg.det(M)))); mats.append(M)
    # dimension de l'espace des formes quadratiques invariantes
    A = []
    for M in mats:
        for i, j in itertools.combinations_with_replacement(range(3), 2):
            pass
    basis = []
    for i, j in itertools.combinations_with_replacement(range(3), 2):
        B = np.zeros((3, 3)); B[i, j] = B[j, i] = 1; basis.append(B)
    rows = []
    for M in mats:
        for B in basis:
            rows.append((M.T @ B @ M - B).ravel())
    Rm = np.array(rows).T  # colonnes = images des formes de base
    # resoudre sum_b x_b (M^T B_b M - B_b) = 0 pour tout M
    Mbig = np.vstack([np.array([(M.T @ B @ M - B).ravel() for B in basis]).T for M in mats])
    dim_inv = 6 - np.linalg.matrix_rank(Mbig)
    print(f"{nom:10s} gens automorphismes : {ok} ; |<gens>| = {len(seen)} ; |stab| = {len(stab)} ; image sur les coins : ordre {len(imgs)}, "
          f"det -1 present : {(-1 in dets)} ; formes quadratiques invariantes : dim {dim_inv}")
