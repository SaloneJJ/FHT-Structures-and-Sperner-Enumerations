#!/usr/bin/env python3
"""AF23 -- Dirac natif de A3 (V352.Hnat, point R-unif, secteur s = +1) sur l'arbre 8-aire referme sur S^3 (bord du 4-cube).
Criteres : 05_audit/CONTINU_S3_14-11.md, section « Dirac natif ... » (14/11/2026, 19 h 55, avant calcul).
Racine : les 8 cubes du bord de [0,1]^4 (facettes (i, s), i = axe normal, s = 0/1) ; niveau n : grille 2^n par facette ; chaque facette a son repere
(ses 3 axes tangents, dans l'ordre croissant) ; un saut qui traverse le bord d'une facette porte la matrice de transport V (automorphisme interieur de Cl(6) :
X_c -> X_c', Y_c -> sigma Y_c', la composante Y changeant de signe avec l'orientation).  Mode « tore » : une facette periodique = V352.Hnat(N, k = 0) (validation).
Usage : python3 AF23_dirac_natif_S3.py valider            (tore, N = 1, 2 : ecart au spectre de V352.Hnat)
        python3 AF23_dirac_natif_S3.py spectre N [k]      (S^3, profondeur N ; k plus petites |E| ; N <= 2 dense, N >= 3 creux)
"""
import sys, os, itertools, time, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, '..', 'V352'), os.path.join(HERE, '..', '..', 'V352')): sys.path.insert(0, p)
import V352_arbre as M
Xg, Yg, G4, G5, FG = M.Xg, M.Yg, M.G4, M.G5, M.FG
S_, MU, R_, TF, TC = +1, 1.0, 1.0, 1.0, 1.0

def h_fwd(d):
    """bloc H[u, v] pour v = voisin avant de u selon la direction locale d (comme Hnat)"""
    return -TC / 2j * Yg[d] - S_ * TC / 2 * Xg[d]

def transport(cA, cB, sig):
    """V unitaire 8x8 avec V X_c V^+ = X_cB[c], V Y_c V^+ = sig[c] Y_cB[c] (c = indices de A)"""
    rows = []
    for c in range(3):
        for P, Q in ((Xg[c], Xg[cB[c]]), (Yg[c], sig[c] * Yg[cB[c]])):
            rows.append(np.kron(np.eye(8), P.T) - np.kron(Q, np.eye(8)))      # vec(V P - Q V) = 0 (vectorisation par lignes)
    A = np.vstack(rows); _, s, Vh = np.linalg.svd(A)
    v = Vh[-1].conj(); V = v.reshape(8, 8)
    V = V / np.sqrt(np.trace(V @ V.conj().T).real / 8)
    assert s[-1] < 1e-8 and np.allclose(V @ V.conj().T, np.eye(8), atol=1e-8), 'transport introuvable'
    return V

def build(N, mode='S3'):
    """renvoie (H creuse ou dense, liste des cellules, profondeur de chaque cellule)"""
    import scipy.sparse as sp
    facets = [(i, s) for i in range(4) for s in (0, 1)] if mode == 'S3' else [(3, 0)]
    tang = {f: [a for a in range(4) if a != f[0]] for f in facets}
    cells = []; index = {}
    for n in range(N + 1):
        m = 2 ** n
        for f in facets:
            for co in itertools.product(range(m), repeat=3):
                index[(n, f, co)] = len(cells); cells.append((n, f, co))
    nc = len(cells); blocks = {}
    def add(u, v, B):
        blocks[(u, v)] = blocks.get((u, v), 0) + B
    Vcache = {}
    links = set()
    for (n, f, co) in cells:
        m = 2 ** n; u = index[(n, f, co)]; tA = tang[f]
        add(u, u, S_ * TF * sum(Xg) + MU * FG)                                   # termes sur site du repere local
        for d in range(3):
            for e in ((+1,) if mode == 'tore' else (+1, -1)):
                c2 = list(co); c2[d] += e
                if 0 <= c2[d] < m or mode == 'tore':
                    c2[d] %= m; v = index[(n, f, tuple(c2))]; W = np.eye(8)
                else:
                    j = tA[d]; g = (j, 1 if e > 0 else 0); tB = tang[g]; i, s = f
                    cB = {}
                    for a, val in zip(tA, co): cB[a] = val
                    del cB[j]; cB[i] = 0 if s == 0 else m - 1
                    v = index[(n, g, tuple(cB[a] for a in tB))]
                    key = (f, g)
                    if key not in Vcache:
                        sigma_j = e * (1 if s == 0 else -1)
                        cmap = [tB.index(i) if a == j else tB.index(a) for a in tA]
                        sig = [sigma_j if a == j else 1 for a in tA]
                        Vcache[key] = transport(None, cmap, sig)
                    W = Vcache[key].conj().T                                       # spineur de v (repere B) -> repere A
                if mode == 'tore':
                    add(u, v, h_fwd(d)); add(v, u, h_fwd(d).conj().T)
                else:
                    lk = (min(u, v), max(u, v), 'lat')
                    if lk in links: continue
                    links.add(lk)
                    B = (h_fwd(d) if e > 0 else h_fwd(d).conj().T) @ W
                    add(u, v, B); add(v, u, B.conj().T)
    # vertical : enfant <- parent (1/sqrt 8), meme facette, meme repere
    for (n, f, co) in cells:
        if n == 0: continue
        ch = index[(n, f, co)]; pa = index[(n - 1, f, tuple(c // 2 for c in co))]
        s8 = 1 / np.sqrt(8)
        # Hnat : kron((S - S^+)/2j, G4) + s (r/2) kron(2 Id - S - S^+, G5), S[ch, pa] = 1/sqrt 8
        add(ch, pa, s8 / 2j * G4 - S_ * R_ / 2 * s8 * G5); add(pa, ch, -s8 / 2j * G4 - S_ * R_ / 2 * s8 * G5)
    for u in range(nc): add(u, u, S_ * R_ / 2 * 2 * G5)
    I, J, Vv = [], [], []
    for (u, v), B in blocks.items():
        r, c = np.nonzero(np.abs(B) > 1e-14)
        I += list(8 * u + r); J += list(8 * v + c); Vv += list(B[r, c])
    H = sp.csr_matrix((Vv, (I, J)), shape=(8 * nc, 8 * nc))
    return H, cells

def valider():
    ok = True
    for N in (1, 2):
        H, _ = build(N, 'tore'); E1 = np.sort(np.linalg.eigvalsh(H.toarray())); E0 = np.sort(np.linalg.eigvalsh(M.Hnat(N)))
        err = np.max(np.abs(E1 - E0)); print(f'tore N = {N} : ecart max au spectre de V352.Hnat = {err:.2e}', flush=True); ok &= err < 1e-8
        print('   hermitien :', np.allclose(H.toarray(), H.toarray().conj().T))
    print('VALIDATION', 'OK' if ok else 'ECHEC'); return ok

def spectre(N, k=60):
    import scipy.sparse.linalg as sla
    t = time.time(); H, cells = build(N, 'S3'); n = H.shape[0]
    print(f'S^3, N = {N} : {len(cells)} cellules, dimension {n}, hermitien {abs(H - H.conj().T).max():.1e} ({time.time() - t:.0f} s)', flush=True)
    if n <= 6000: E = np.linalg.eigvalsh(H.toarray())
    else: E = sla.eigsh(H.tocsc(), k=k, sigma=1e-6, which='LM', return_eigenvectors=False)
    E = np.sort(E[np.argsort(np.abs(E))][:k])
    a = np.sort(np.abs(E)); packs = []; cur = [a[0]]
    for x in a[1:]:
        if x <= 1.2 * cur[-1] and x - cur[0] <= 0.2 * cur[0]: cur.append(x)
        else: packs.append(cur); cur = [x]
    packs.append(cur)
    print('  plus petites |E| :', np.round(a[:24], 5))
    print('  paquets (taille, moyenne, rapport au premier) :', [(len(p), round(float(np.mean(p)), 5), round(float(np.mean(p) / np.mean(packs[0])), 3)) for p in packs[:6]])
    print('  signes du premier paquet :', int(np.sum(E[np.abs(E) <= 1.2 * packs[0][-1]] > 0)), '+ /', int(np.sum(E[np.abs(E) <= 1.2 * packs[0][-1]] < 0)), '-')
    print(f'  ({time.time() - t:.0f} s)', flush=True)
    json.dump({'N': N, 'E': [float(x) for x in E]}, open(os.path.join(HERE, f'AF23_S3_N{N}.json'), 'w'))

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: print(__doc__)
    elif a[0] == 'valider': valider()
    elif a[0] == 'spectre': spectre(int(a[1]), int(a[2]) if len(a) > 2 else 60)
