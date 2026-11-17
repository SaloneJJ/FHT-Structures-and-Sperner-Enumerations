#!/usr/bin/env python3
"""PROPOSITION P2 -- coefficient d'Einstein-Hilbert induit par le Dirac natif de l'arbre, avec controle nul.
Criteres : 05_audit/AUDIT_FABLE_PROPOSITIONS.md, P2 (13/11, 20 h 30) et son AMENDEMENT du 14/11/2026, 7 h 55 (avant tout calcul).

Modele : D natif de la piste B (V352_arbre.Hnat, point R-unif t_f = t_c = r = mu = 1, secteur gamma_F = +1), arbre augmente de profondeur N,
foret d'arbres (Bloch k sur la periode de la racine, P = 2^N cellules fines).  Modulation t(x) = 1 + a cos(2 pi x / L) des termes LATERAUX :
saut i-j x sqrt(t_i t_j), terme sur site t_f X_d x t_i ; iso = les trois directions ; nul = direction x seule ; vertical, r, mu inchanges.
z(t) = Tr exp(-t D^2) par cellule fine ;  c2(t) = [ln z(a) - ln z(0)] / a^2  (limite a -> 0).
Calcul exact du coefficient de a^2 (formule de Daleckii-Krein) : avec H_m la partie modulee,
  H(a) = H0 + (a/2){C, H_m} - (a^2/8)[C,[C,H_m]] + O(a^3),  C = cos(qx), q = 2 pi / L ;
  d^2/da^2 Tr f(H) / 2 = Tr f'(H0) V2_diag + sum_k sum_{m in k, n in k+qP} f'[1](E_m, E_n) |<n| B_k |m>|^2,
  V2_diag(k) = -(1/8) [H_m(k) - (H_m(k+qP) + H_m(k-qP))/2],  B_k = (H_m(k) + H_m(k+qP)) / 4,  f(E) = exp(-t E^2).
Prediction (continu, Weyl a deux composantes, a1 = -R/12) : c2(t) = c_vol - (t/12)(2 pi / L)^2 + ...  (iso) ; pente nulle (controle nul).

Usage :
  python3 P2_EH_induit.py                               SEQUENCE COMPLETE : validation, puis N = 2 (grilles 8 et 12), bilans (~ 30 min a 1 h)
  python3 P2_EH_induit.py tout N3                       idem, plus N = 3 (grille 4 ; ~ 1 a 2 h de plus, ~ 3 Go)
  python3 P2_EH_induit.py bis                           P2-BIS (criteres du 14/11, 13 h 46) : q commensurables, tores 96 et 128, N = 2 (~ 3 h 30)
  python3 P2_EH_induit.py valider                       validation du code (N=2, L=16, a=0.01, super-cellule de 4 arbres) : ecart < 2 %
  python3 P2_EH_induit.py calcul N m [i/n] [sortie.npz] sommes partielles (L = 64 et 128 ; iso et nul) ; i/n = morceau i (0..n-1) de la grille
  python3 P2_EH_induit.py bilan fichiers.npz ...        combine, ajuste c2 = alpha + beta/t + gamma t sur [15, 40], applique les criteres
Exemples (Mac, venv) : depuis 03_code/CalculCodes/AUDIT_FABLE :
  python3 P2_EH_induit.py valider
  python3 P2_EH_induit.py calcul 2 8 0/1 P2_N2.npz && python3 P2_EH_induit.py bilan P2_N2.npz
  python3 P2_EH_induit.py calcul 3 4 0/1 P2_N3.npz && python3 P2_EH_induit.py bilan P2_N3.npz     (N = 3 : ~ 1 a 2 h, ~ 3 Go)
numpy seul.  Ecrit seulement le fichier .npz demande.
"""
import sys, os, itertools, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'V352')); sys.dont_write_bytecode = True
import V352_arbre as M

TS = np.array([1, 2, 3, 5, 7, 10, 12, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40.])
LS = (64, 128)
MODES = {'iso': (0, 1, 2), 'nul': (0,)}
S_, MU, R_, TF, TC = +1, 1.0, 1.0, 1.0, 1.0

def lateral(T, d, Nc, W=None):
    """terme lateral de la direction d (sauts + sur site t_f X_d), comme dans Hnat ; W = diag(sqrt t) ou None"""
    Td = T[d] if W is None else (W[:, None] * T[d] * W[None, :])
    on = np.eye(Nc) if W is None else np.diag(W ** 2)
    H = -TC * np.kron((Td - Td.conj().T) / 2j, M.Yg[d])
    H += S_ * (TF * np.kron(on, M.Xg[d]) - TC * np.kron((Td + Td.conj().T) / 2, M.Xg[d]))
    return H

def rest(S, Nc):
    return np.kron((S - S.conj().T) / 2j, M.G4) + S_ * (R_ / 2) * np.kron(2 * np.eye(Nc) - S - S.conj().T, M.G5) + MU * np.kron(np.eye(Nc), M.FG)

def parts(N, k):
    T, S = M.cell_ops(N, k); Nc = T[0].shape[0]
    lat = [lateral(T, d, Nc) for d in range(3)]
    return lat, rest(S, Nc)

def fprime(E, t): return -2 * t * E * np.exp(-t * E ** 2)
def fsecond(E, t): return (-2 * t + 4 * t * t * E * E) * np.exp(-t * E ** 2)

def dd(Em, En, t):
    """differences divisees f'[1](E_m, E_n), matrice (n, m)"""
    X = En[:, None]; Y = Em[None, :]; D = X - Y; near = np.abs(D) < 1e-9
    out = np.where(near, 0.0, (fprime(X, t) - fprime(Y, t)) / np.where(near, 1.0, D))
    if near.any(): out = np.where(near, fsecond(0.5 * (X + Y), t), out)
    return out

def grid(N, m):
    return [2 * np.pi * np.array(j) / m for j in itertools.product(range(m), repeat=3)]

def calcul(N, m, part='0/1', out=None, drop_zero=True, Ls=LS, ks=None, verbose=True):
    i, n = map(int, part.split('/'))
    ks = grid(N, m) if ks is None else ks
    ks = ks[i::n]; P = 2 ** N
    z0 = np.zeros(len(TS)); acc = {(mo, L): np.zeros(len(TS)) for mo in MODES for L in Ls}; nzero = 0; t0 = time.time()
    for ik, k in enumerate(ks):
        lat, rs = parts(N, k); H0 = sum(lat) + rs
        E, U = np.linalg.eigh(H0)
        Ek, Uk = E, U
        if drop_zero and np.allclose(np.asarray(k) % (2 * np.pi), 0):
            # retire les 2 modes nuls de la racine (comme V352) : dans l'espace nul, les 2 vecteurs de plus grand poids sur la racine
            z = np.abs(E) < 1e-9; Z = U[:, z]
            root = np.zeros(len(E)); root[:8] = 1                      # cellule de profondeur 0 = 8 premieres composantes
            w, R = np.linalg.eigh(Z.conj().T @ (root[:, None] * Z))
            nd = 2 if Z.shape[1] >= 2 else 0
            Zk = Z @ R[:, :Z.shape[1] - nd]                              # on garde les vecteurs nuls de plus faible poids sur la racine
            Ek = np.concatenate([E[~z], np.zeros(Zk.shape[1])]); Uk = np.hstack([U[:, ~z], Zk]); nzero += nd
            if verbose: print(f'  k = 0 : {int(z.sum())} modes nuls, poids racine {np.round(w, 3)}, {nd} retires', flush=True)
        z0 += np.array([np.exp(-t * Ek ** 2).sum() for t in TS])
        for L in Ls:
            qP = 2 * np.pi * P / L
            sh = {}
            for sgn in (+1, -1):
                kk = np.array(k, float); kk[0] += sgn * qP; sh[sgn] = parts(N, kk)[0]
            Ep, Up = np.linalg.eigh(sum(sh[+1]) + parts(N, np.array(k) + np.array([qP, 0, 0]))[1])
            for mo, dirs in MODES.items():
                Hm = sum(lat[d] for d in dirs); Hp = sum(sh[+1][d] for d in dirs); Hmn = sum(sh[-1][d] for d in dirs)
                V2 = -(Hm - 0.5 * (Hp + Hmn)) / 8
                B = (Hm + Hp) / 4
                v2d = np.real(np.einsum('ij,ij->j', Uk.conj(), V2 @ Uk))         # <m|V2|m>
                Bnm = np.abs(Up.conj().T @ B @ Uk) ** 2                            # (n in k+qP, m in k)
                acc[(mo, L)] += np.array([(fprime(Ek, t) * v2d).sum() + (dd(Ek, Ep, t) * Bnm).sum() for t in TS])
        if verbose and (ik % 16 == 0 or ik == len(ks) - 1):
            print(f'  k {ik + 1}/{len(ks)}  ({time.time() - t0:.0f} s)', flush=True)
    res = dict(N=N, m=m, part=part, nk=len(ks), z0=z0, nzero=nzero, TS=TS, Ls=np.array(Ls),
               **{f'd_{mo}_{L}': v for (mo, L), v in acc.items()})
    if out: np.savez(out, **res)
    return res

def bilan(files):
    R = [dict(np.load(f, allow_pickle=True)) for f in files]
    N = int(R[0]['N']); nk = sum(int(r['nk']) for r in R); z0 = sum(r['z0'] for r in R); Ls = [int(x) for x in R[0]['Ls']]
    assert nk == int(R[0]['m']) ** 3, f'grille incomplete : {nk} points sur {int(R[0]["m"]) ** 3}'
    zc = z0 / nk / 8 ** N
    print(f'=== P2 : N = {N}, grille {int(R[0]["m"])}^3, modes nuls de la racine retires : {sum(int(r["nzero"]) for r in R)}')
    print('  t      n_eff = z (4 pi t)^{3/2}   ' + '   '.join(f'c2_{mo}(L={L})' for L in Ls for mo in MODES))
    c2 = {(mo, L): sum(r[f'd_{mo}_{L}'] for r in R) / z0 for mo in MODES for L in Ls}
    for j, t in enumerate(TS):
        print(f'  {t:5.1f}  {zc[j] * (4 * np.pi * t) ** 1.5:8.4f}          ' + '   '.join(f'{c2[(mo, L)][j]:+.6f}' for L in Ls for mo in MODES))
    sel = (TS >= 15) & (TS <= 40); A = np.vstack([np.ones(sel.sum()), 1 / TS[sel], TS[sel]]).T
    ok_all = True; valid = True; null_result = True
    for L in Ls:
        pred = -(1 / 12) * (2 * np.pi / L) ** 2
        g = {mo: np.linalg.lstsq(A, c2[(mo, L)][sel], rcond=None)[0][2] for mo in MODES}
        r_iso = g['iso'] / pred; r_nul = abs(g['nul']) / max(abs(g['iso']), 1e-300)
        v_ok = r_nul < 0.2; p_ok = g['iso'] < 0 and abs(r_iso - 1) <= 0.3; z_ok = abs(g['iso']) < 0.2 * abs(pred)
        valid &= v_ok; ok_all &= p_ok; null_result &= z_ok
        print(f'  L = {L} : pente iso gamma = {g["iso"]:+.3e}, prediction {pred:+.3e} (rapport {r_iso:.3f}) ; pente nul = {g["nul"]:+.3e} '
              f'(|nul|/|iso| = {r_nul:.3f}) -> controle nul {"conforme" if v_ok else "ECHOUE"} ; iso {"conforme" if p_ok else "non conforme"}')
    if not valid: v = 'TEST INVALIDE (le controle nul ne s\'annule pas : la mesure n\'isole pas la courbure)'
    elif ok_all: v = 'le Dirac natif induit le terme d\'Einstein-Hilbert d\'un Weyl de reseau [E] (nature R : Sakharov)'
    elif null_result: v = 'pas de terme induit a cet ordre (resultat negatif propre a l\'arbre)'
    else: v = 'ni conforme ni nul : non concluant, le dire'
    print('VERDICT P2 :', v)

# ---------------------------------------------------------------- validation par diagonalisation directe d'une super-cellule
def cell_ops_super(N, Mx, k):
    """comme M.cell_ops, mais Mx arbres le long de x (periode Mx 2^N en x, 2^N en y, z) ; renvoie T, S, positions x des cellules"""
    L = 2 ** N; per = np.array([Mx * L, L, L], float); k = np.asarray(k, float)
    sides = [2 ** n for n in range(N + 1)]; dims = [(Mx * s, s, s) for s in sides]
    offs = np.concatenate([[0], np.cumsum([a * b * c for a, b, c in dims])]).astype(int); Nc = int(offs[-1])
    T = [np.zeros((Nc, Nc), complex) for _ in range(3)]; S = np.zeros((Nc, Nc), complex); X = np.zeros(Nc)
    idx = lambda n, I, J, K: offs[n] + (I * dims[n][1] + J) * dims[n][2] + K
    for n, s in enumerate(sides):
        d = L // s; nx = dims[n][0]
        for I, J, K in itertools.product(range(nx), range(s), range(s)):
            i = idx(n, I, J, K); X[i] = (I + 0.5) * d
            for jd, sh in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
                nb = idx(n, (I + sh[0]) % nx, (J + sh[1]) % s, (K + sh[2]) % s)
                T[jd][i, nb] += np.exp(-1j * k[jd] * d / per[jd])
    for n in range(1, N + 1):
        s = sides[n]; d = L // s
        for I, J, K in itertools.product(range(Mx * s), range(s), range(s)):
            ch = idx(n, I, J, K); pa = idx(n - 1, I // 2, J // 2, K // 2)
            dx = np.array([(I % 2) - 0.5, (J % 2) - 0.5, (K % 2) - 0.5]) * d
            S[ch, pa] = np.exp(1j * np.sum(k * dx / per)) / np.sqrt(8)
    return T, S, X

def zdirect(N, Mx, Lm, a, mode, kgrid):
    zs = np.zeros(len(TS))
    for k in kgrid:
        T, S, X = cell_ops_super(N, Mx, k); Nc = T[0].shape[0]
        W = np.sqrt(1 + a * np.cos(2 * np.pi * X / Lm))
        H = sum(lateral(T, d, Nc, W if d in MODES[mode] else None) for d in range(3)) + rest(S, Nc)
        E = np.linalg.eigvalsh(H); zs += np.array([np.exp(-t * E ** 2).sum() for t in TS])
    return zs

def valider():
    N, Lm, a, my = 2, 16, 0.01, 2; P = 2 ** N; Mx = Lm // P; mx = 2
    # controle de coherence : super-cellule a 1 arbre = Hnat
    T, S, X = cell_ops_super(N, 1, (0.3, 0.2, 0.1)); Nc = T[0].shape[0]
    assert np.allclose(sum(lateral(T, d, Nc) for d in range(3)) + rest(S, Nc), M.Hnat(N, (0.3, 0.2, 0.1), +1)), 'constructeur != Hnat'
    ksup = [2 * np.pi * np.array([i / mx, j / my, l / my]) for i, j, l in itertools.product(range(mx), range(my), range(my))]
    # grille equivalente pour un arbre : kx = (kx_s + 2 pi j) / Mx
    kone = [np.array([(ks[0] + 2 * np.pi * j) / Mx, ks[1], ks[2]]) for ks in ksup for j in range(Mx)]
    ok = True
    for mode in MODES:
        t0 = time.time()
        z0 = zdirect(N, Mx, Lm, 0.0, mode, ksup); za = zdirect(N, Mx, Lm, a, mode, ksup)
        fd = (np.log(za) - np.log(z0)) / a ** 2
        r = calcul(N, 1, '0/1', None, drop_zero=False, Ls=(Lm,), ks=kone, verbose=False)
        pt = r[f'd_{mode}_{Lm}'] / r['z0']
        rel = np.abs(fd - pt) / np.maximum(np.abs(pt), 1e-12)
        sel = TS <= 40
        print(f'  validation {mode} (N = 2, L = 16, a = 0.01) : ecart relatif max {rel[sel].max():.2e} ; c2(t=10) direct {fd[5]:+.5f}, perturbatif {pt[5]:+.5f} ({time.time() - t0:.0f} s)')
        ok &= rel[sel].max() < 0.02
    print('VALIDATION :', 'conforme (< 2 %)' if ok else 'ECHEC : calcul invalide')
    return ok

# ---------------------------------------------------------------- P2-bis : q commensurables (criteres du 14/11/2026, 13 h 46)
TSB = TS[TS >= 10]

def eig_sector(N, k, drop_zero, verbose=False):
    """valeurs et vecteurs propres de H0(k) ; si k = 0 mod 2 pi (toutes composantes), retire les 2 modes nuls de plus grand poids sur la racine"""
    lat, rs = parts(N, k); E, U = np.linalg.eigh(sum(lat) + rs); nd = 0
    if drop_zero and np.allclose(np.mod(np.asarray(k, float) + 1e-12, 2 * np.pi), 0, atol=1e-9):
        z = np.abs(E) < 1e-9; Z = U[:, z]; root = np.zeros(len(E)); root[:8] = 1
        w, R = np.linalg.eigh(Z.conj().T @ (root[:, None] * Z)); nd = 2 if Z.shape[1] >= 2 else 0
        Zk = Z @ R[:, :Z.shape[1] - nd]
        E = np.concatenate([E[~z], np.zeros(Zk.shape[1])]); U = np.hstack([U[:, ~z], Zk])
    return lat, E, U, nd

def bis_calcul(N, m, part='0/1', out=None, drop_zero=True, Ecut=1.5, verbose=True):
    """P2-bis : tore T = m 2^N, q_j = 2 pi j / T (j = 1, 2, 3) ; par ligne (ky, kz), une seule diagonalisation par point kx (k + q_j sur la ligne) ;
    paires (m, n) restreintes a celles ou |E_m| < Ecut ou |E_n| < Ecut (erreur < exp(-10 Ecut^2) pour t >= 10)."""
    i0, n0 = map(int, part.split('/')); T = m * 2 ** N; J = (1, 2, 3); Ls = tuple(T / j for j in J)
    lines = [(a, b) for a in range(m) for b in range(m)][i0::n0]
    z0 = np.zeros(len(TSB)); acc = {(mo, L): np.zeros(len(TSB)) for mo in MODES for L in Ls}; nzero = 0; t0 = time.time()
    for il, (iy, iz) in enumerate(lines):
        sec = {}; Hm = {}
        for ix in range(-3, m + 3):
            k = 2 * np.pi * np.array([ix, iy, iz]) / m
            if 0 <= ix < m + 3:
                lat, E, U, nd = eig_sector(N, k, drop_zero); sec[ix] = (E, U)
                if 0 <= ix < m: nzero += nd
            else:
                lat = parts(N, k)[0]
            Hm[ix] = {mo: sum(lat[d] for d in dirs) for mo, dirs in MODES.items()}
        for ix in range(m):
            E, U = sec[ix]; z0 += np.array([np.exp(-t * E ** 2).sum() for t in TSB])
            lo = np.abs(E) < Ecut
            for j, L in zip(J, Ls):
                Ep, Up = sec[ix + j]; lop = np.abs(Ep) < Ecut
                for mo in MODES:
                    V2 = -(Hm[ix][mo] - 0.5 * (Hm[ix + j][mo] + Hm[ix - j][mo])) / 8
                    v2d = np.real(np.einsum('ij,ij->j', U[:, lo].conj(), V2 @ U[:, lo]))
                    B = (Hm[ix][mo] + Hm[ix + j][mo]) / 4
                    M1 = np.abs(Up.conj().T @ B @ U[:, lo]) ** 2                   # (n tous, m bas)
                    M2 = np.abs(Up[:, lop].conj().T @ B @ U[:, ~lo]) ** 2          # (n bas, m hauts)
                    acc[(mo, L)] += np.array([(fprime(E[lo], t) * v2d).sum() + (dd(E[lo], Ep, t) * M1).sum()
                                              + (dd(E[~lo], Ep[lop], t) * M2).sum() for t in TSB])
        if verbose and (il % 8 == 0 or il == len(lines) - 1):
            print(f'  ligne {il + 1}/{len(lines)}  ({time.time() - t0:.0f} s)', flush=True)
    res = dict(N=N, m=m, part=part, nk=len(lines) * m, z0=z0, nzero=nzero, TS=TSB, Ls=np.array(Ls), **{f'd_{mo}_{L}': v for (mo, L), v in acc.items()})
    if out: np.savez(out, **res)
    return res

def bis_valider():
    """bis_calcul (lignes, troncature) == calcul (generique) sur un petit tore, a t >= 10"""
    N, m = 2, 6; T = m * 4
    r1 = bis_calcul(N, m, '0/1', None, drop_zero=False, verbose=False)
    ks = grid(N, m); r2 = calcul(N, m, '0/1', None, drop_zero=False, Ls=tuple(T / j for j in (1, 2, 3)), ks=ks, verbose=False)
    sel = TS >= 10; err = 0
    for mo in MODES:
        for L in r1['Ls']:
            a = r1[f'd_{mo}_{L}'] / r1['z0']; b = r2[f'd_{mo}_{L}'][sel] / r2['z0'][sel]
            err = max(err, np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-12)))
    print(f'  validation de P2-bis (lignes + troncature contre calcul generique, tore 24) : ecart relatif max {err:.2e}', '-> conforme' if err < 1e-6 else '-> ECHEC')
    return err < 1e-6

def bis_bilan(files):
    R = [dict(np.load(f, allow_pickle=True)) for f in files]
    N = int(R[0]['N']); m = int(R[0]['m']); Ls = [float(x) for x in R[0]['Ls']]; T = m * 2 ** N
    nk = sum(int(r['nk']) for r in R); assert nk == m ** 3, f'grille incomplete : {nk} / {m ** 3}'
    z0 = sum(r['z0'] for r in R)
    c2 = {(mo, L): sum(r[f'd_{mo}_{L}'] for r in R) / z0 for mo in MODES for L in Ls}
    q2 = np.array([(2 * np.pi / L) ** 2 for L in Ls])
    V = np.vstack([np.ones(3), q2, q2 ** 2]).T
    B = {mo: np.array([np.linalg.solve(V, np.array([c2[(mo, L)][i] for L in Ls]))[1] for i in range(len(R[0]['TS']))]) for mo in MODES}
    print(f'=== P2-bis : N = {N}, tore T = {T} (grille {m}^3), q_j = 2 pi j / T, j = 1, 2, 3')
    TSl = R[0]['TS']
    print('  t      n_eff     B_iso(t)     B_nul(t)')
    for i, t in enumerate(TSl):
        print(f'  {t:5.1f}  {z0[i] / nk / 8 ** N * (4 * np.pi * t) ** 1.5:7.4f}  {B["iso"][i]:+10.4f}  {B["nul"][i]:+10.4f}')
    sel = (TSl >= 15) & (TSl <= 40); A = np.vstack([np.ones(sel.sum()), TSl[sel]]).T
    b1 = {mo: np.linalg.lstsq(A, B[mo][sel], rcond=None)[0][1] for mo in MODES}
    print(f'  pente b1 sur [15, 40] : iso {b1["iso"]:+.4f}, nul {b1["nul"]:+.4f} ; prediction iso {-1 / 12:+.4f}, nul 0')
    return b1

def bis_verdict(b1s):
    """b1s : liste des dictionnaires b1 des deux tores"""
    iso = [b['iso'] for b in b1s]; nul = [b['nul'] for b in b1s]; pred = -1 / 12
    res = all(abs(x - y) <= 0.2 * max(abs(x), abs(y)) for x, y in ((iso[0], iso[1]), (nul[0], nul[1]))) if len(b1s) == 2 else None
    valid = all(abs(n) < 0.2 * abs(i) for i, n in zip(iso, nul))
    ok = all(i < 0 and 0.7 <= i / pred <= 1.3 for i in iso)
    neg = all(abs(i) < 0.2 / 12 for i in iso)
    print(f'  resolution (accord des deux tores a 20 %) : {res} ; validite (|nul| < 0,2 |iso|) : {valid}')
    if res is False: v = 'NON RESOLU (les deux tores ne s\'accordent pas) : le dire'
    elif not valid: v = 'TEST INVALIDE (le controle nul ne s\'annule pas)'
    elif ok: v = 'le Dirac natif induit le terme d\'Einstein-Hilbert d\'un Weyl de reseau [E] (nature R : gravite induite)'
    elif neg: v = 'pas de terme induit a cet ordre (resultat negatif propre a l\'arbre)'
    else: v = 'non concluant (ni la prediction a 30 %, ni zero) : le dire'
    print('VERDICT P2-bis :', v)

def bis_tout():
    print('=== P2-bis, sequence complete', time.strftime('%d/%m/%Y %H:%M'), flush=True)
    if not bis_valider(): print('arret : validation echouee'); return
    b1s = []
    for m in (24, 32):
        out = os.path.join(HERE, f'P2bis_N2_m{m}.npz')
        if not os.path.exists(out):
            print(f'--- calcul N = 2, tore {4 * m} (grille {m}^3) -> {out}', flush=True); bis_calcul(2, m, '0/1', out)
        b1s.append(bis_bilan([out]))
    bis_verdict(b1s)

def tout(N3=False):
    print('=== P2, sequence complete', time.strftime('%d/%m/%Y %H:%M'), flush=True)
    if not valider(): print('arret : validation echouee'); return
    for N, m in ([(2, 8), (2, 12)] + ([(3, 4)] if N3 else [])):
        out = os.path.join(HERE, f'P2_N{N}_m{m}.npz')
        print(f'--- calcul N = {N}, grille {m}^3 -> {out}', flush=True)
        calcul(N, m, '0/1', out); bilan([out])
    print('Lecture : la grille 8 et la grille 12 doivent donner les memes pentes (sinon la grille est trop grossiere : le dire).')

if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'tout'
    if cmd == 'tout': tout(N3=('N3' in sys.argv[2:]))
    elif cmd == 'valider': valider()
    elif cmd == 'calcul':
        N, m = int(sys.argv[2]), int(sys.argv[3]); part = sys.argv[4] if len(sys.argv) > 4 else '0/1'
        out = sys.argv[5] if len(sys.argv) > 5 else f'P2_N{N}_m{m}_{part.replace("/", "de")}.npz'
        calcul(N, m, part, out); print('->', out)
    elif cmd == 'bilan': bilan(sys.argv[2:])
    elif cmd == 'bis': bis_tout()
    elif cmd == 'bis_calcul':
        N, m = int(sys.argv[2]), int(sys.argv[3]); part = sys.argv[4] if len(sys.argv) > 4 else '0/1'
        out = sys.argv[5] if len(sys.argv) > 5 else f'P2bis_N{N}_m{m}_{part.replace("/", "de")}.npz'
        bis_calcul(N, m, part, out); print('->', out)
    elif cmd == 'bis_bilan': bis_verdict([bis_bilan(sys.argv[2:])])
