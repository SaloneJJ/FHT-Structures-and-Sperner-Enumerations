#!/usr/bin/env python3
"""D2 -- C3' sur une racine fermee : Dirac natif de A3 (V352.Hnat, point R-unif, jauge alternee) sur l'OREILLER (double du cube = deux cubes A, B
recolles face a face, atlas oriente : transitions = demi-tours), compare au TORE revetement (8 cubes racines periodiques, T^3 de cote 2).
Criteres : FHT_28/09_racine_S3/CRITERES_D2_oreiller.md (ecrits le 15/11/2026, 8 h 40, AVANT calcul).
Transports exacts dans Cl(6) : V_f = s_f Y_a Y_b (demi-tour autour de l'axe d de la face, {a, b} = les deux autres axes), s_f = +-1 (signe de la face).
Usage :
  python3 D2_oreiller.py controles                  (C1 tore = V352 Bloch k0 in {0,pi}^3 ; C2 jauge ; C3 hermitien, Gamma4 conserve)
  python3 D2_oreiller.py signes                     (principe de signes : holonomies des 12 aretes, 32 classes)
  python3 D2_oreiller.py mesure N [classe]          (R(Lambda), w_N, Gamma4, rho ; secteurs s = +1, -1 ; mu = 1 et controle C4 mu = 0)
Sorties : D2_N<N>.json + log (| tee)."""
import sys, os, itertools, json, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, '..', '06_code', 'V352')): sys.path.insert(0, p)
import V352_arbre as M
Xg, Yg, G4, G5, FG = M.Xg, M.Yg, M.G4, M.G5, M.FG
I8 = np.eye(8)
BLOCH_K0 = [0.0, 0.0, 0.0]

def h_fwd(d, s, tc=1.0):
    """bloc H[u, v], v = voisin avant de u selon l'axe local d (identique a V352.Hnat en jauge alternee)"""
    return -tc / 2j * Yg[d] - s * tc / 2 * Xg[d]

def V_face(d, sign=1):
    a, b = [c for c in range(3) if c != d]
    return sign * Yg[a] @ Yg[b]

# ------------------------------------------------------------------ geometrie
def cells_of(N, cubes):
    cells = []; index = {}
    for n in range(N + 1):
        m = 2 ** n
        for q in cubes:
            for co in itertools.product(range(m), repeat=3):
                index[(n, q, co)] = len(cells); cells.append((n, q, co))
    return cells, index

def neighbour(mode, n, q, co, d, e, signs):
    """voisin lateral de (n, q, co) dans la direction locale d, sens e ; renvoie (q', co', W) avec W = transport (spineur de v -> repere de u)
    et le sens du lien dans le repere de v (le lien est 'avant' selon d dans les deux reperes : e est le meme)."""
    m = 2 ** n; c2 = list(co); c2[d] += e
    if mode == 'bloch':                        # un cube racine periodique, condition de Bloch k0 (jauge des bords ; equivalente a V352)
        if 0 <= c2[d] < m: return q, tuple(c2), I8
        c2[d] %= m; return q, tuple(c2), np.exp(-1j * e * BLOCH_K0[d]) * I8
    if mode == 'tore':                         # 8 cubes q in {0,1}^3, grille globale periodique de cote 2m
        g = [q[i] * m + co[i] for i in range(3)]; g[d] = (g[d] + e) % (2 * m)
        return tuple(x // m for x in g), tuple(x % m for x in g), I8
    if 0 <= c2[d] < m: return q, tuple(c2), I8
    # oreiller : on sort de q par la face (d, e) ; on entre dans l'autre cube ; reperes : A = coordonnees x, B = y = m-1-x (repere miroir)
    # la face (d, e) dans les coordonnees du cube : x_d = m-1 (e = +1 dans A) ; dans le repere de B la meme face est y_d = 0 (meme sens e)
    q2 = 'B' if q == 'A' else 'A'
    co2 = tuple((0 if e > 0 else m - 1) if i == d else m - 1 - co[i] for i in range(3))
    face = (d, e) if q == 'A' else (d, -e)       # face nommee dans les coordonnees de A
    V = V_face(d, signs[face])                   # V : operateurs du repere A -> repere B
    W = V.conj().T if q == 'A' else V            # u dans A, v dans B : spineur de B -> A = V^+ ; u dans B, v dans A : V
    return q2, co2, W

def build(N, mode='oreiller', s=+1, mu=1.0, signs=None, r=1.0, tf=1.0, tc=1.0):
    """H creuse (8 x nombre de cellules), meme convention que V352.Hnat (jauge alternee t_c -> -t_c)"""
    if signs is None: signs = {(d, e): 1 for d in range(3) for e in (1, -1)}
    cubes = ['A', 'B'] if mode == 'oreiller' else (['R'] if mode == 'bloch' else list(itertools.product((0, 1), repeat=3)))
    cells, index = cells_of(N, cubes); nc = len(cells)
    blocks = {}
    def add(u, v, B): blocks[(u, v)] = blocks.get((u, v), 0) + B
    TCs = -tc
    for (n, q, co) in cells:
        u = index[(n, q, co)]
        add(u, u, s * tf * sum(Xg) + mu * FG + s * r * G5)
        for d in range(3):
            q2, co2, W = neighbour(mode, n, q, co, d, +1, signs)
            v = index[(n, q2, co2)]
            B = (TCs / 2j * Yg[d] + s * TCs / 2 * Xg[d]) @ W     # = h_fwd(d) en jauge alternee
            add(u, v, B); add(v, u, B.conj().T)
    s8 = 1 / np.sqrt(8)
    for (n, q, co) in cells:
        if n == 0: continue
        ch = index[(n, q, co)]; pa = index[(n - 1, q, tuple(c // 2 for c in co))]
        add(ch, pa, s8 / 2j * G4 - s * r / 2 * s8 * G5); add(pa, ch, -s8 / 2j * G4 - s * r / 2 * s8 * G5)
    I, J, Vv = [], [], []
    for (u, v), B in blocks.items():
        rr, cc = np.nonzero(np.abs(B) > 1e-14)
        I += list(8 * u + rr); J += list(8 * v + cc); Vv += list(B[rr, cc])
    return sp.csr_matrix((Vv, (I, J)), shape=(8 * nc, 8 * nc)), cells

# ------------------------------------------------------------------ controles
def controles():
    ok = True
    for N in (1, 2):
        for s in (1, -1):
            H, _ = build(N, 'tore', s=s); E1 = np.sort(np.linalg.eigvalsh(H.toarray()))
            E0 = np.sort(np.concatenate([np.linalg.eigvalsh(M.Hnat(N, k0=np.pi * np.array(k), s=s)) for k in itertools.product((0, 1), repeat=3)]))
            err = np.max(np.abs(E1 - E0)); ok &= err < 1e-10
            print(f'C1 tore N = {N}, s = {s:+d} : ecart a la reunion des 8 secteurs de Bloch de V352.Hnat = {err:.1e}', flush=True)
            Eb = np.sort(np.concatenate([np.linalg.eigvalsh(build_bloch(N, k, s=s).toarray()) for k in itertools.product((0, 1), repeat=3)]))
            err = np.max(np.abs(Eb - E0)); ok &= err < 1e-10
            print(f'C1 mode Bloch creux N = {N}, s = {s:+d} : ecart a V352.Hnat = {err:.1e}', flush=True)
    for N in (1, 2):
        H, _ = build(N); herm = abs(H - H.conj().T).max(); ok &= herm < 1e-12
        sg = {(d, e): -1 for d in range(3) for e in (1, -1)}
        H2, _ = build(N, signs=sg); g = np.max(np.abs(np.sort(np.linalg.eigvalsh(H.toarray())) - np.sort(np.linalg.eigvalsh(H2.toarray()))))
        ok &= g < 1e-10
        print(f'C2/C3 oreiller N = {N} : hermitien {herm:.1e} ; jauge (signe global de B) ecart {g:.1e}', flush=True)
    c3 = all(np.allclose(V_face(d) @ G4 @ V_face(d).conj().T, G4) and np.allclose(V_face(d) @ FG @ V_face(d).conj().T, FG) for d in range(3))
    ok &= c3; print('C3 : V_f Gamma4 V_f^+ = Gamma4 et V_f FGamma5 V_f^+ = FGamma5 pour les faces :', c3)
    print('CONTROLES', 'OK' if ok else 'ECHEC', flush=True); return ok

# ------------------------------------------------------------------ principe de signes
def S_rot(a, b, th):
    """relevement spinoriel de la rotation d'angle th de e_a vers e_b : S Y_a S^+ = cos Y_a + sin Y_b"""
    return np.cos(th / 2) * I8 - np.sin(th / 2) * Yg[a] @ Yg[b]

def edge_holonomies(signs):
    """N = 1 (m = 2) ; pour chaque arete du cube A (axe d, position (pa, pb) dans {0, m-1}^2 sur les deux autres axes), boucle autour de l'arete :
    on part d'une cellule de A au coin de l'arete, on traverse la face a vers B, on revient par la face b dans A, on referme dans A.
    Holonomie P = produit des transports (repere de depart) ; orientation : sens de rotation de la boucle dans le plan (a, b) du repere de A."""
    N = 1; m = 2; out = {}
    for d in range(3):
        a, b = [c for c in range(3) if c != d]
        for ea, eb in itertools.product((1, -1), repeat=2):
            co = [0, 0, 0]; co[a] = m - 1 if ea > 0 else 0; co[b] = m - 1 if eb > 0 else 0
            co = tuple(co); P = I8.copy(); q = 'A'; c = co
            # etape 1 : traverser la face (a, ea) : transport du spineur de u vers v = W(v<-u) = W(u<-v)^+
            q2, c2, W = neighbour('oreiller', N, q, c, a, ea, signs); P = W.conj().T @ P; q, c = q2, c2
            # dans B : la face b de A est, dans le repere de B, la face (b, -eb) en coordonnees y ; on la traverse dans le sens ou elle mene a A
            # (dans le repere de B le sens de traversee de la face b est -eb : y_b = m-1-x_b)
            q2, c2, W = neighbour('oreiller', N, q, c, b, -eb, signs); P = W.conj().T @ P; q, c = q2, c2
            assert q == 'A'
            # sens de la boucle dans le plan (a, b) de A : sortie par +ea selon a, retour par la face b (venant de l'exterieur, sens -eb)
            # angle de sortie : direction ea e_a ; angle d'entree : direction eb e_b ; la boucle tourne de la direction (ea e_a) vers (eb e_b) par l'exterieur
            # (le quadrant exterieur au cube) ; on l'ecrit comme rotation dans le plan (a, b) : sens + si det[[ea, 0], [0, eb]] < 0 (voir texte)
            orient = ea * eb                     # +1 : boucle dans le sens e_a -> e_b (trigonometrique dans le plan (a, b)) ; -1 : sens inverse
            Splus = S_rot(a, b, orient * np.pi)   # relevement continu de la rotation de +pi dans le sens de la boucle (deficit pi)
            if np.allclose(P, Splus): tau = +1
            elif np.allclose(P, -Splus): tau = -1
            else: tau = 0
            out[(d, ea, eb)] = dict(tau=tau, faces=((a, ea), (b, eb)))
    return out

def signes():
    base = {(d, e): 1 for d in range(3) for e in (1, -1)}
    H0 = edge_holonomies(base)
    print('holonomies d\'arete (signes tous +1) : tau = +1 si P = S(+pi) dans le sens de la boucle, -1 si P = -S(+pi), 0 si ni l\'un ni l\'autre')
    for k, v in H0.items(): print('   arete', k, 'faces', v['faces'], 'tau =', v['tau'])
    faces = sorted(base); res = {+1: [], -1: []}
    for bits in itertools.product((1, -1), repeat=5):            # signe de la face (0, +1) fixe a +1 (jauge : signe global de B)
        sg = dict(zip(faces, (1,) + bits))
        H = edge_holonomies(sg)
        for sense in (+1, -1):
            if all(v['tau'] == sense for v in H.values()): res[sense].append(sg)
    for sense in (+1, -1):
        print(f'classes ou les 12 holonomies valent {"+" if sense > 0 else "-"}S(+pi) (sens de la boucle) : {len(res[sense])}', [ {f'{k[0]}{"+" if k[1] > 0 else "-"}': v for k, v in sg.items()} for sg in res[sense] ])
    return res

# ------------------------------------------------------------------ mesures
def build_bloch(N, k, s=+1, mu=1.0):
    BLOCH_K0[:] = [np.pi * x for x in k]; H, _ = build(N, 'bloch', s=s, mu=mu); BLOCH_K0[:] = [0.0] * 3; return H

def levels(N, mode):
    cubes = ['A', 'B'] if mode == 'oreiller' else ['R']
    cells, _ = cells_of(N, cubes); return cells, np.array([c[0] for c in cells])

def low_modes(H, lev, N, stop, kmin=48):
    """modes de plus petites |E| (inversion spectrale en 0, k doublant) jusqu'a ce que stop(E, wN, w0) soit vrai ; dense si dim <= 3000.
    wN = poids au niveau fin, w0 = poids a la racine (niveau 0)."""
    n = H.shape[0]; fine = np.repeat(lev == N, 8); root = np.repeat(lev == 0, 8)
    if n <= 3000:
        E, U = np.linalg.eigh(H.toarray())
    else:
        k = kmin
        while True:
            k = min(k, n - 2); E, U = sla.eigsh(H.tocsc(), k=k, sigma=0.0, which='LM')
            o = np.argsort(np.abs(E)); E, U = E[o], U[:, o]
            if stop(E, (np.abs(U[fine]) ** 2).sum(0), (np.abs(U[root]) ** 2).sum(0)) or k >= n - 2: break
            k *= 2
    o = np.argsort(np.abs(E)); E, U = E[o], U[:, o]
    return E, U, (np.abs(U[fine]) ** 2).sum(0), (np.abs(U[root]) ** 2).sum(0)

def mesure(N, classe=None):
    """Amendements DECLARES apres le premier calcul (N = 2, 15/11 ~ 9 h 30) :
    (A1) modes de la contrepartie a la racine (poids >= 0,5 au niveau 0) exclus de E_bulk et des comptes (comme A3, Prop. Heat kernel) ;
    (A2) la chiralite par <Gamma4> est abandonnee (nulle par construction : Gamma4 est la matrice du saut en profondeur) ; voir CRITERES §7 ;
    (A3) supplementaire, non decisif : comptage lisse Z(t) = sum exp(-t E^2) sur les modes du niveau fin, t = 1/Lambda_j^2, R_lisse = 4 Z_or / Z_tore."""
    t0 = time.time(); signs = classe or {(d, e): 1 for d in range(3) for e in (1, -1)}
    res = {'N': N, 'signes': {f'{k[0]}{k[1]:+d}': v for k, v in signs.items()}}
    cR, lR = levels(N, 'bloch'); cO, lO = levels(N, 'oreiller')
    m = 2 ** N; nr = lambda x: x in (0, 1, m - 2, m - 1)
    nl = np.array([c[0] == N and any(nr(c[2][i]) and nr(c[2][j]) for i, j in ((0, 1), (0, 2), (1, 2))) for c in cO])
    frac = nl.sum() / (lO == N).sum(); nearO = np.repeat(nl, 8)
    print(f'N = {N} : dimension tore {8 * len(cR)} par secteur de Bloch (x 8) ; oreiller {8 * len(cO)} ; fraction de volume fine pres des lignes {frac:.3f}', flush=True)
    for mu in (1.0, 0.0):
        for s in (+1, -1):
            Et, wt, w0t = [], [], []
            for k in itertools.product((0, 1), repeat=3):
                H = build_bloch(N, k, s=s, mu=mu)
                E, U, wN, w0 = low_modes(H, lR, N, lambda E, wN, w0: np.any((wN < 0.5) & (w0 < 0.5)) or np.abs(E).max() > 1.5)
                Et.append(E); wt.append(wN); w0t.append(w0)
            Et, wt, w0t = np.abs(np.concatenate(Et)), np.concatenate(wt), np.concatenate(w0t)
            keep_t = w0t < 0.5
            bulk = Et[(wt < 0.5) & keep_t]; Eb = float(bulk.min()) if len(bulk) else float('nan')
            Lmax = 0.5 * Eb; L = [j / 4 * Lmax for j in (2, 3, 4)]
            Ho, _ = build(N, 'oreiller', s=s, mu=mu, signs=signs)
            Eo, Uo, wo, w0o = low_modes(Ho, lO, N, lambda E, wN, w0: np.abs(E).max() > Lmax and len(E) >= 12)
            Eo = np.abs(Eo); keep_o = w0o < 0.5
            rho = ((np.abs(Uo[nearO]) ** 2).sum(0) / np.maximum(wo, 1e-12)) / frac
            key = f'mu={mu:g} s={s:+d}'; r = dict(E_bulk=Eb, Lambda=L)
            r['racine_exclus'] = dict(tore=int(np.sum(~keep_t & (Et < Lmax))), oreiller=int(np.sum(~keep_o & (Eo < Lmax))))
            r['N_tore'] = [int(np.sum(keep_t & (Et < x))) for x in L]; r['N_oreiller'] = [int(np.sum(keep_o & (Eo < x))) for x in L]
            r['R'] = [4 * no / nt if nt else None for no, nt in zip(r['N_oreiller'], r['N_tore'])]
            ft = keep_t & (wt >= 0.5); fo = keep_o & (wo >= 0.5)
            r['R_lisse'] = [float(4 * np.sum(np.exp(-Eo[fo] ** 2 / x ** 2)) / np.sum(np.exp(-Et[ft] ** 2 / x ** 2))) for x in L]
            so = keep_o & (Eo < Lmax)
            r['oreiller'] = dict(n=int(so.sum()))
            if so.any():
                r['oreiller'].update(frac_wN_ge_08=float(np.mean(wo[so] >= 0.8)), rho_moy=float(np.mean(rho[so])), rho_max=float(np.max(rho[so])))
            r['oreiller'].update(plus_basses=np.round(Eo[:16], 4).tolist(), wN=np.round(wo[:16], 3).tolist(), rho=np.round(rho[:16], 2).tolist())
            r['tore'] = dict(n=int(np.sum(keep_t & (Et < Lmax))), plus_basses=np.round(np.sort(Et[keep_t])[:16], 4).tolist())
            res[key] = r
            print(f'== N = {N}, {key} : E_bulk(tore) = {Eb:.4f}, Lambda_j = {np.round(L, 4).tolist()} ; racine exclus {r["racine_exclus"]}', flush=True)
            print(f'   N_tore = {r["N_tore"]}, N_oreiller = {r["N_oreiller"]}, R = {[None if x is None else round(x, 3) for x in r["R"]]}, R_lisse = {np.round(r["R_lisse"], 3).tolist()}')
            print('   tore    :', r['tore']); print('   oreiller:', r['oreiller'], f'({time.time() - t0:.0f} s)', flush=True)
    json.dump(res, open(os.path.join(HERE, f'D2_N{N}.json'), 'w'), indent=1, default=float)
    return res

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: print(__doc__)
    elif a[0] == 'controles': controles()
    elif a[0] == 'signes': signes()
    elif a[0] == 'mesure':
        cl = None
        if len(a) > 2: cl = {(int(k[0]), int(k[1:])): v for k, v in json.loads(a[2]).items()}
        mesure(int(a[1]), cl)
