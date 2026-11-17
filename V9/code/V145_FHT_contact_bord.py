"""V145 : chantier 1 (V144 §4) -- mouvements lateraux, arbre augmente et bord quotient.
E1 : substitution sigma_F de K26 (article IV) avec les regles laterales R0a (freres), R0b (+ cousins germains),
     R1 (contact par la realisation, aretes relevees e x V_F) : composantes et diametres des graphes horizontaux.
E2 : substitution par les sommets de contact (8 coins x, faces stella (B) de V87), construite COMBINATOIREMENT
     (union-find sur les coins, sans coordonnees) : (i) identification diagonale, (ii) recollement des copies voisines,
     (iii) heritage des recollements de faces (cousins de tout ordre). Tests M1 (croissance), M2 (liens, homologie
     du quotient : cube et fermeture T3), M3 (dimension spectrale), coherence de (i) avec (ii)-(iii).
E3 : sigma_F complete (26 sommets substitues) avec contacts par les coins seulement : fraction spatiale.
JJ SALONE, FHT 2026"""
import itertools, time, sys
from collections import defaultdict, Counter, deque
import numpy as np
from V86_FHT_k26 import build_K
from V87_FHT_stella import rules_stella, faces_of

t00 = time.time()
K, I, II = build_K()
V = sorted(set().union(*K)); nV = len(V); idx = {v: i for i, v in enumerate(V)}
print(f'K26 (modele en joins) : {nV} sommets, {len(K)} facettes', flush=True)

# ---------- 1-squelette de K26 (hote de la racine de sigma_F : graine plate) ----------
adj = np.zeros((nV, nV), bool)
for f in K:
    l = [idx[v] for v in f]
    for a in l:
        for b in l:
            if a != b: adj[a, b] = True
def bfs_dist(A, s):
    d = [-1] * len(A); d[s] = 0; q = deque([s])
    while q:
        u = q.popleft()
        for w in np.nonzero(A[u])[0]:
            if d[w] < 0: d[w] = d[u] + 1; q.append(w)
    return d
D1 = np.array([bfs_dist(adj, s) for s in range(nV)])
print(f'1-squelette de K26 : {adj.sum()//2} aretes sur {nV*(nV-1)//2}, diametre {D1.max()}, rayon {D1.max(1).min()}', flush=True)
nonadj = [(V[a], V[b]) for a in range(nV) for b in range(a+1, nV) if not adj[a, b]]
print('  paires non adjacentes :', len(nonadj), nonadj[:12], flush=True)

# ---------- E1 : sigma_F, regles R0a, R0b, R1 ----------
print('\n=== E1 : sigma_F de K26 (26 enfants par noeud, hote = 1-squelette de K26) ===', flush=True)
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
def words(n): return list(itertools.product(range(nV), repeat=n))
for n in (1, 2, 3):
    W = words(n); N = len(W); pos = {w: i for i, w in enumerate(W)}
    for regle in ('R0a', 'R0b'):
        rows, cols = [], []
        for w in W:
            i = pos[w]
            # freres : meme parent, adjacents dans l'hote
            for b in range(nV):
                if b != w[-1] and adj[w[-1], b]:
                    rows.append(i); cols.append(pos[w[:-1] + (b,)])
            if regle == 'R0b' and n >= 2:   # cousins germains : meme grand-parent, parents distincts (tous, comme l'article)
                for a in range(nV):
                    if a == w[-2]: continue
                    for b in range(nV):
                        rows.append(i); cols.append(pos[w[:-2] + (a, b)])
        G = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N)).tocsr()
        nc, lab = connected_components(G, directed=False)
        sizes = Counter(Counter(lab).values())
        # diametre de la composante du noeud 0
        comp0 = np.nonzero(lab == lab[0])[0]
        sub = G[comp0][:, comp0]
        dm = shortest_path(sub, unweighted=True, directed=False, indices=list(range(min(30, len(comp0))))).max()
        attendu = {'R0a': 26 ** (n - 1), 'R0b': 26 ** max(n - 2, 0) if n >= 2 else 1}[regle]
        print(f'  n={n} {regle} : {N} noeuds, {nc} composantes (attendu {attendu}), tailles {dict(sizes)}, diametre (composante de la racine, echantillon) {int(dm)}', flush=True)
    # R1 : contact par la realisation : u ~ v ssi les enfants de lca sur les chemins sont adjacents dans l'hote.
    # distance exacte = distance dans l'hote entre ces deux enfants (formule), verifiee par BFS implicite a n = 2.
    if n <= 2:
        rows, cols = [], []
        for w in W:
            for w2 in W:
                if w2 == w: continue
                m = next(j for j in range(n) if w[j] != w2[j])
                if adj[w[m], w2[m]]: rows.append(pos[w]); cols.append(pos[w2])
        G = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N)).tocsr()
        nc, lab = connected_components(G, directed=False)
        dm = shortest_path(G, unweighted=True, directed=False).max()
        print(f'  n={n} R1 : {nc} composante(s), diametre {int(dm)} (borne V83 : diam du 1-squelette = {D1.max()})', flush=True)
    else:
        print(f'  n={n} R1 : diametre = max de la distance dans l\'hote entre enfants de l\'ancetre commun = {D1.max()} (formule, independante de n)', flush=True)

# ---------- E2 : substitution par les coins, construction combinatoire ----------
print('\n=== E2 : substitution par les sommets de contact (coins x, faces stella (B)) ===', flush=True)
RB = rules_stella()                     # RB[k] = [(v dans F_k^+, tau_k(v) dans F_k^-)]
Fam = faces_of(RB)                       # [F_0^+, F_0^-, F_1^+, F_1^-, F_2^+, F_2^-]
C = sorted(set().union(*Fam)); nC = len(C); ci = {c: i for i, c in enumerate(C)}
Fp = [set(Fam[2*k]) for k in range(3)]; Fm = [set(Fam[2*k+1]) for k in range(3)]
tau = [dict(RB[k]) for k in range(3)]
print('  coins :', C)
for k in range(3): print(f'  direction {k} : F+ = {sorted(Fp[k])}, F- = {sorted(Fm[k])}, tau = {tau[k]}')
# signes combinatoires (appartenance aux faces), uniquement pour reperer les coins voisins dans la mere
sgn = {c: tuple(1 if c in Fp[k] else -1 for k in range(3)) for c in C}
assert len(set(sgn.values())) == 8, 'coins non distingues par les faces'
for k in range(3):   # tau_k = translation : envoie F+ sur F- en conservant les autres signes
    for v, w in tau[k].items():
        assert sgn[w] == tuple(-s if j == k else s for j, s in enumerate(sgn[v])), 'tau non compatible'
print('  controle : 8 coins distingues par les faces ; chaque tau_k conserve les deux autres signes (Q3)', flush=True)
voisin = {}   # voisin[(c, k)] = coin c' de la mere situe en direction +k de c (si sgn_k(c) = -1)
for c in C:
    for k in range(3):
        if sgn[c][k] == -1:
            voisin[(c, k)] = next(d for d in C if sgn[d] == tuple(-s if j == k else s for j, s in enumerate(sgn[c])))
class UF:
    def __init__(s): s.p = {}
    def f(s, a):
        p = s.p
        while p.get(a, a) != a:
            p[a] = p.get(p[a], p[a]); a = p[a]
        return a
    def u(s, a, b):
        a, b = s.f(a), s.f(b)
        if a != b: s.p[max(a, b)] = min(a, b)
def construire(nmax, tore=False):
    """retourne, par profondeur n, (cellules, union-find des coins, liste des paires recollees (u, v, k))"""
    res = {}
    glued = [((), (), k) for k in range(3)] if tore else []   # fermeture T3 : la racine recollee a elle-meme
    cells = [()]
    for n in range(1, nmax + 1):
        cells = [w + (c,) for w in cells for c in C]
        new = []
        for w in set(x[:-1] for x in cells):          # (ii) freres voisins
            for (c, k), c2 in voisin.items(): new.append((w + (c,), w + (c2,), k))
        for (u, v, k) in glued:                        # (iii) heritage : enfants sur F_k^+ de u et leurs images
            for c in Fp[k]: new.append((u + (c,), v + (tau[k][c],), k))
        glued = new
        uf = UF()
        for (u, v, k) in glued:
            for x in Fp[k]: uf.u((u, x), (v, tau[k][x]))
        res[n] = (cells, uf, glued)
    return res
for tore in (False, True):
    nom = 'T3 (fermeture periodique)' if tore else 'cube (sans fermeture)'
    print(f'\n  --- {nom} ---', flush=True)
    t0 = time.time(); R = construire(5 if not tore else 4, tore)
    for n, (cells, uf, glued) in R.items():
        cl = {(w, x): uf.f((w, x)) for w in cells for x in C}
        nsom = len(set(cl.values()))
        inj = all(len({cl[(w, x)] for x in C}) == 8 for w in cells)
        # coherence de la diagonale (i) avec la profondeur precedente
        coh = None
        if n >= 2:
            cellsP, ufP, _ = R[n-1]; ok = True
            grp = defaultdict(list)
            for w in cellsP:
                for x in C: grp[ufP.f((w, x))].append((w, x))
            for members in grp.values():
                rep = cl[(members[0][0] + (members[0][1],), members[0][1])]
                for (w, x) in members[1:]:
                    if cl[(w + (x,), x)] != rep: ok = False; break
                if not ok: break
            coh = ok
        att = 8 ** n if tore else (2 ** n + 1) ** 3
        print(f'  n={n} : {len(cells)} cellules, {len(glued)} paires recollees, {nsom} classes de coins (attendu {att}), injectif par cellule {inj}, diagonale coherente {coh}', flush=True)
    # graphes de contact (partage d'au moins un coin) et M1
    for n, (cells, uf, glued) in R.items():
        if n < 2: continue
        pos = {w: i for i, w in enumerate(cells)}
        byc = defaultdict(list)
        for w in cells:
            for x in C: byc[uf.f((w, x))].append(pos[w])
        rows, cols = [], []
        for lst in byc.values():
            for a in lst:
                for b in lst:
                    if a != b: rows.append(a); cols.append(b)
        G = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(cells),) * 2).tocsr(); G.data[:] = 1
        nc, _ = connected_components(G, directed=False)
        deg = Counter(np.diff(G.indptr))
        src = list(range(0, len(cells), max(1, len(cells) // 40)))
        dm = int(shortest_path(G, unweighted=True, directed=False, indices=src).max())
        dlog = np.log(len(cells)) / np.log(dm) if dm > 1 else float('nan')
        print(f'  M1 n={n} : composantes {nc}, N_n = {len(cells)}, diametre {dm} (attendu {2**n - 1 if not tore else 2**(n-1)}), ln N / ln diam = {dlog:.3f}, degres {dict(sorted(deg.items()))}', flush=True)
    # M2 : complexe cubique (n = 2, 3), liens des sommets et homologie Z/2
    for n in [m for m in R if m in (2, 3)]:
        cells, uf, glued = R[n]
        cube = [frozenset(uf.f((w, x)) for x in C) for w in cells]
        aretes_tuile = [(a, b) for a, b in itertools.combinations(C, 2) if sum(sa != sb for sa, sb in zip(sgn[a], sgn[b])) == 1]
        faces_tuile = [tuple(sorted(F_)) for k in range(3) for F_ in (Fp[k], Fm[k])]
        E = set(); Q = set()
        for w in cells:
            for a, b in aretes_tuile: E.add(frozenset((uf.f((w, a)), uf.f((w, b)))))
            for F_ in faces_tuile: Q.add(frozenset(uf.f((w, x)) for x in F_))
        Vs = sorted(set().union(*cube)); vi = {v: i for i, v in enumerate(Vs)}
        El = sorted(E, key=lambda e: tuple(sorted(e))); ei = {e: i for i, e in enumerate(El)}
        Ql = sorted(Q, key=lambda q: tuple(sorted(q))); qi = {q: i for i, q in enumerate(Ql)}
        Cl = cube
        def rank2(rows_list):
            piv = {}; r = 0
            for col in rows_list:
                while col:
                    h = col.bit_length() - 1
                    if h in piv: col ^= piv[h]
                    else: piv[h] = col; r += 1; break
            return r
        Eby = defaultdict(list)
        for e in El:
            for v in e: Eby[v].append(e)
        d1 = [(1 << vi[a]) | (1 << vi[b]) for a, b in (tuple(e) for e in El)]
        d2 = []
        for q in Ql:
            es = {e for v in q for e in Eby[v] if e <= q}
            assert len(es) == 4, 'face carree attendue'
            d2.append(sum(1 << ei[e] for e in es))
        Qby = defaultdict(list)
        for q in Ql:
            for v in q: Qby[v].append(q)
        d3 = []
        for cbe in Cl:
            fs = {q for v in cbe for q in Qby[v] if q <= cbe}
            assert len(fs) == 6, 'cube a 6 faces attendu'
            d3.append(sum(1 << qi[q] for q in fs))
        r1, r2, r3 = rank2(d1), rank2(d2), rank2(d3)
        b0 = len(Vs) - r1; b1 = len(El) - r1 - r2; b2 = len(Ql) - r2 - r3; b3 = len(Cl) - r3
        chi = len(Vs) - len(El) + len(Ql) - len(Cl)
        # liens : pour chaque sommet, un triangle par cube contenant le sommet (les 3 aretes du cube en ce sommet)
        Eat = defaultdict(set)
        for e in El:
            for v in e: Eat[v].add(e)
        Cby = defaultdict(list)
        for cbe in Cl:
            for v in cbe: Cby[v].append(cbe)
        types = Counter(); ok_var = True
        for v in Vs:
            tris = []
            for cbe in Cby[v]:
                tris.append(frozenset(e for e in Eat[v] if e <= cbe))
            assert all(len(t) == 3 for t in tris)
            sv = set().union(*tris); sed = Counter(frozenset(p) for t in tris for p in itertools.combinations(t, 2))
            nv_, ne_, nt_ = len(sv), len(sed), len(tris)
            bord = sum(1 for c in sed.values() if c == 1); interieur = all(c in (1, 2) for c in sed.values())
            chi_l = nv_ - ne_ + nt_
            typ = ('sphere' if bord == 0 and chi_l == 2 and nt_ == 8 else
                   'disque' if bord > 0 and chi_l == 1 and interieur else 'AUTRE')
            if typ == 'AUTRE': ok_var = False
            types[(typ, nt_)] += 1
        print(f'  M2 n={n} : V={len(Vs)} E={len(El)} F={len(Ql)} C={len(Cl)}, chi={chi}, Betti Z/2 = ({b0},{b1},{b2},{b3}) ; liens (type, nb de cubes) : {dict(types)} ; variete combinatoire : {ok_var}', flush=True)
    # M3 : dimension spectrale (tore seulement, n = 4) par loi de comptage du laplacien du graphe de contact
    if tore:
        cells, uf, glued = R[4]
        pos = {w: i for i, w in enumerate(cells)}
        A = np.zeros((len(cells),) * 2)
        byc = defaultdict(list)
        for w in cells:
            for x in C: byc[uf.f((w, x))].append(pos[w])
        for lst in byc.values():
            for a in lst:
                for b in lst:
                    if a != b: A[a, b] = 1
        L = np.diag(A.sum(1)) - A
        ev = np.sort(np.linalg.eigvalsh(L))
        # N(lambda) ~ lambda^{d_s/2} sur les petites valeurs propres (hors zero)
        lam = ev[1:]; Ncount = np.arange(1, len(lam) + 1)
        sel = (lam > 0) & (lam < np.quantile(lam, 0.08))
        pente = np.polyfit(np.log(lam[sel]), np.log(Ncount[sel]), 1)[0]
        print(f'  M3 (T3, n=4, {len(cells)} cellules) : pente de log N(lambda) = {pente:.3f}, d_s = {2*pente:.3f} (attendu 3) ; noyau {int((ev < 1e-9).sum())}', flush=True)

# ---------- E3 : sigma_F complete (26 sommets substitues), contacts par les coins seulement ----------
print('\n=== E3 : sigma_F complete avec contacts par les coins (sans P-contact) ===', flush=True)
Cset = set(C)
for n in (1, 2, 3, 4):
    # un noeud est spatial ssi tous ses ancetres stricts et lui-meme sont des coins ; les autres n'ont aucun coin identifie.
    # avec adjacence fraternelle (hote) en plus, la composante geante = noeuds dont les n-1 premieres lettres sont des coins.
    tot = 26 ** n; spat = 8 ** n; geant = 8 ** (n - 1) * 26
    print(f'  n={n} : noeuds {tot}, spatiaux {spat} (fraction {spat/tot:.4f}), composante geante avec freres {geant} (fraction {geant/tot:.4f} = (8/26)^{n-1})', flush=True)
# verification explicite a n = 2 et 3 : graphe = contacts par coins (E2) + adjacence fraternelle dans l'hote
for n in (2, 3):
    W = list(itertools.product(range(nV), repeat=n)); pos = {w: i for i, w in enumerate(W)}
    Rn = construire(n)[n]; cells, uf, _ = Rn
    rows, cols = [], []
    byc = defaultdict(list)
    for w in cells:
        for x in C: byc[uf.f((w, x))].append(pos[tuple(idx[c] for c in w)])
    for lst in byc.values():
        for a in lst:
            for b in lst:
                if a != b: rows.append(a); cols.append(b)
    for w in W:
        for b in range(nV):
            if b != w[-1] and adj[w[-1], b]: rows.append(pos[w]); cols.append(pos[w[:-1] + (b,)])
    G = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(W),) * 2).tocsr()
    nc, lab = connected_components(G, directed=False)
    big = Counter(lab).most_common(1)[0][1]
    print(f'  explicite n={n} : {nc} composantes, geante {big} sur {len(W)} (fraction {big/len(W):.4f}) ; ilots de taille {Counter(Counter(lab).values()).most_common(3)}', flush=True)
print(f'\nduree totale {time.time()-t00:.0f} s')
