#!/usr/bin/env python3
"""D2, version 2 (15/11/2026, 14 h 30) -- meme modele et memes mesures que D2_oreiller.py (importe), calcul accelere pour N = 4.
Changements DECLARES (apres le calcul N <= 3, avant N = 4 ; aucun critere de lecture modifie) :
  (a) une seule factorisation LU de (H - sigma) par matrice (sigma = 1e-4), reutilisee quand on demande plus de valeurs propres ;
      k initial dimensionne d'avance (tore : 160 par secteur ; oreiller : 64) ;
  (b) cas complet seulement pour (mu = 1, s = +1) ; (mu = 1, s = -1) : controle que les plus basses |E| de l'oreiller coincident avec s = +1
      (identiques a 1e-10 pres a N = 2 et 3) ; (mu = 0, s = +1) : controle C4 sur l'oreiller seul, avec le Lambda_max du cas (mu = 1, s = +1) ;
  (c) avancement imprime secteur par secteur ;
  (d) dans un groupe de valeurs propres degenerees, rotation qui diagonalise le poids a la racine avant le tri racine / niveau fin (A1 applique au sous-espace).
Usage : python3 D2_oreiller_v2.py N | tee D2_N<N>_v2.log     (sortie : D2_N<N>_v2.json)"""
import sys, os, itertools, json, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import D2_oreiller as D

CLASSE = {(0, -1): 1, (0, 1): -1, (1, -1): -1, (1, 1): 1, (2, -1): 1, (2, 1): -1}     # classe fixee par le principe de signes (D2_signes.log)
SIGMA = 1e-4

def modes(H, lev, N, k0, stop, tag):
    """plus petites |E| par inversion spectrale ; LU calcule une fois ; k double jusqu'a stop(E, wN, w0)"""
    t0 = time.time(); n = H.shape[0]; fine = np.repeat(lev == N, 8); root = np.repeat(lev == 0, 8)
    lu = sla.splu((H - SIGMA * sp.identity(n, dtype=complex, format='csc')).tocsc())
    op = sla.LinearOperator((n, n), matvec=lu.solve, dtype=complex)
    k = k0
    while True:
        k = min(k, n - 2)
        E, U = sla.eigs(H, k=k, sigma=SIGMA, OPinv=op, which='LM')
        E = E.real; o = np.argsort(E); E, U = E[o], U[:, o]
        # groupes degeneres (|E_i - E_j| < 1e-8) : rotation qui diagonalise le poids a la racine (separe la contrepartie de la racine, amendement A1)
        i = 0
        while i < len(E):
            j = i + 1
            while j < len(E) and abs(E[j] - E[i]) < 1e-8: j += 1
            if j - i > 1:
                Uc, _ = np.linalg.qr(U[:, i:j]); Mr = Uc[root].conj().T @ Uc[root]; _, V = np.linalg.eigh(Mr); U[:, i:j] = Uc @ V
            i = j
        o = np.argsort(np.abs(E)); E, U = E[o], U[:, o]
        wN = (np.abs(U[fine]) ** 2).sum(0); w0 = (np.abs(U[root]) ** 2).sum(0)
        if stop(np.abs(E), wN, w0) or k >= n - 2: break
        k *= 2
    print(f'      {tag} : n = {n}, k = {k}, |E| max = {np.abs(E).max():.4f}, {time.time() - t0:.0f} s', flush=True)
    return np.abs(E), U, wN, w0

def run(N):
    t0 = time.time(); res = {'N': N, 'version': 2, 'signes': {f'{d}{e:+d}': v for (d, e), v in CLASSE.items()}}
    cR, lR = D.levels(N, 'bloch'); cO, lO = D.levels(N, 'oreiller')
    m = 2 ** N; nr = lambda x: x in (0, 1, m - 2, m - 1)
    nl = np.array([c[0] == N and any(nr(c[2][i]) and nr(c[2][j]) for i, j in ((0, 1), (0, 2), (1, 2))) for c in cO])
    frac = nl.sum() / (lO == N).sum(); nearO = np.repeat(nl, 8)
    print(f'N = {N} (version 2) : tore {8 * len(cR)} par secteur (x 8) ; oreiller {8 * len(cO)} ; fraction fine pres des lignes {frac:.3f}', flush=True)
    # ---- cas complet : mu = 1, s = +1
    print('== mu = 1, s = +1 : tore (8 secteurs de Bloch)', flush=True)
    Et, wt, w0t = [], [], []
    for k in itertools.product((0, 1), repeat=3):
        H = D.build_bloch(N, k, s=+1, mu=1.0)
        E, U, wN, w0 = modes(H, lR, N, 160, lambda E, wN, w0: np.any((wN < 0.5) & (w0 < 0.5)), f'secteur {k}')
        Et.append(E); wt.append(wN); w0t.append(w0)
    Et, wt, w0t = np.concatenate(Et), np.concatenate(wt), np.concatenate(w0t); keep_t = w0t < 0.5
    Eb = float(Et[(wt < 0.5) & keep_t].min()); Lmax = 0.5 * Eb; L = [j / 4 * Lmax for j in (2, 3, 4)]
    # completude du tore sous Lmax : chaque secteur a ete calcule au-dela d'un mode de bulk (>= Eb > Lmax)
    print(f'   E_bulk(tore) = {Eb:.4f} ; Lambda_j = {np.round(L, 4).tolist()} ({time.time() - t0:.0f} s)', flush=True)
    print('== mu = 1, s = +1 : oreiller', flush=True)
    Ho, _ = D.build(N, 'oreiller', s=+1, mu=1.0, signs=CLASSE)
    Eo, Uo, wo, w0o = modes(Ho, lO, N, 64, lambda E, wN, w0: E.max() > Lmax and len(E) >= 12, 'oreiller')
    keep_o = w0o < 0.5
    rho = ((np.abs(Uo[nearO]) ** 2).sum(0) / np.maximum(wo, 1e-12)) / frac
    r = dict(E_bulk=Eb, Lambda=L)
    r['racine_exclus'] = dict(tore=int(np.sum(~keep_t & (Et < Lmax))), oreiller=int(np.sum(~keep_o & (Eo < Lmax))))
    r['N_tore'] = [int(np.sum(keep_t & (Et < x))) for x in L]; r['N_oreiller'] = [int(np.sum(keep_o & (Eo < x))) for x in L]
    r['R'] = [4 * no / nt if nt else None for no, nt in zip(r['N_oreiller'], r['N_tore'])]
    ft = keep_t & (wt >= 0.5); fo = keep_o & (wo >= 0.5)
    r['R_lisse'] = [float(4 * np.sum(np.exp(-Eo[fo] ** 2 / x ** 2)) / np.sum(np.exp(-Et[ft] ** 2 / x ** 2))) for x in L]
    so = keep_o & (Eo < Lmax); r['oreiller'] = dict(n=int(so.sum()))
    if so.any(): r['oreiller'].update(frac_wN_ge_08=float(np.mean(wo[so] >= 0.8)), rho_moy=float(np.mean(rho[so])), rho_max=float(np.max(rho[so])))
    r['oreiller'].update(plus_basses=np.round(Eo[:16], 4).tolist(), wN=np.round(wo[:16], 3).tolist(), rho=np.round(rho[:16], 2).tolist())
    r['tore'] = dict(n=int(np.sum(keep_t & (Et < Lmax))), plus_basses=np.round(np.sort(Et[keep_t])[:16], 4).tolist())
    res['mu=1 s=+1'] = r
    print(f'   N_tore = {r["N_tore"]}, N_oreiller = {r["N_oreiller"]}, R = {[None if x is None else round(x, 3) for x in r["R"]]}, R_lisse = {np.round(r["R_lisse"], 3).tolist()}')
    print('   tore    :', r['tore']); print('   oreiller:', r['oreiller'], f'({time.time() - t0:.0f} s)', flush=True)
    json.dump(res, open(os.path.join(HERE, f'D2_N{N}_v2.json'), 'w'), indent=1, default=float)
    # ---- controle s = -1 (oreiller)
    print('== controle mu = 1, s = -1 : oreiller', flush=True)
    Hm, _ = D.build(N, 'oreiller', s=-1, mu=1.0, signs=CLASSE)
    Em, _, _, _ = modes(Hm, lO, N, 64, lambda E, wN, w0: len(E) >= 24, 'oreiller s = -1')
    nn = min(24, len(Em), len(Eo)); ecart = float(np.max(np.abs(np.sort(Em)[:nn] - np.sort(Eo)[:nn])))
    res['controle s=-1'] = dict(ecart_24_plus_basses=ecart, plus_basses=np.round(np.sort(Em)[:16], 4).tolist())
    print(f'   ecart des 24 plus basses |E| (s = -1 contre s = +1) : {ecart:.2e} ({time.time() - t0:.0f} s)', flush=True)
    json.dump(res, open(os.path.join(HERE, f'D2_N{N}_v2.json'), 'w'), indent=1, default=float)
    # ---- controle C4 : mu = 0, oreiller, Lambda_max du cas mu = 1
    print('== controle C4 : mu = 0, s = +1, oreiller', flush=True)
    Hz, _ = D.build(N, 'oreiller', s=+1, mu=0.0, signs=CLASSE)
    Ez, Uz, wz, w0z = modes(Hz, lO, N, 64, lambda E, wN, w0: E.max() > Lmax, 'oreiller mu = 0')
    bad = int(np.sum((Ez < Lmax) & (wz >= 0.5) & (w0z < 0.5)))
    res['C4 mu=0'] = dict(Lambda_max=Lmax, modes_fins_sous_Lambda_max=bad, plus_basses=np.round(Ez[:12], 4).tolist(), wN=np.round(wz[:12], 3).tolist())
    print(f'   C4 : modes de poids fin >= 0,5 sous Lambda_max = {bad} ({"passe" if bad == 0 else "ECHEC"}) ; plus basses {np.round(Ez[:8], 4).tolist()} ({time.time() - t0:.0f} s)', flush=True)
    json.dump(res, open(os.path.join(HERE, f'D2_N{N}_v2.json'), 'w'), indent=1, default=float)
    print('FIN', flush=True)

if __name__ == '__main__':
    run(int(sys.argv[1]))
