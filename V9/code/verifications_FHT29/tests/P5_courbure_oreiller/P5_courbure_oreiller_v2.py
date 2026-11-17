#!/usr/bin/env python3
"""
VERSION 2 (16/11/2026, 13 h 50 AST) : seul changement, le solveur de low() passe de eigs (matrices quelconques) a eigsh
(hermitien, meme inversion spectrale), apres un blocage de plus de 5 h a N = 4 ; criteres et mesures inchanges.
P5 -- couplage a la courbure du Dirac natif de A3 sur la racine fermee (oreiller), compare a un Wilson-Dirac ordinaire (4 composantes)
sur la meme geometrie (niveau fin de l'oreiller, meme atlas oriente, meme principe de signes).  Criteres : CRITERES_P5.md (ecrits avant).
Usage :
  python3 P5_courbure_oreiller.py valider            controles V1-V3 (N <= 2, < 1 min)
  python3 P5_courbure_oreiller.py mesure N [K]       N = 3 ou 4 ; K modes natifs (defaut 160), 2K modes Wilson -> P5_N<N>.json
  python3 P5_courbure_oreiller.py wilson N [K]       Wilson seul sur l'oreiller (controle de reseau, N = 5 : L = 32) -> P5_W<N>.json
  Fenetre de comparaison : niveaux natifs sous 0,75 E_bulk (E_bulk lu dans le JSON de D2) ; amendement de CRITERES_P5.md.
Depend de FHT_28/09_racine_S3/D2_oreiller.py (lecture seule).  numpy + scipy."""
import os, sys, json, time, itertools
import numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
D2DIR = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..', 'FHT_28', '09_racine_S3'))
sys.path.insert(0, D2DIR)
import D2_oreiller as D
try:
    import scipy.sparse as sp, scipy.sparse.linalg as sla
    HAVE_SP = True
except ImportError:
    HAVE_SP = False

s0 = np.eye(2); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex)
SIG = [sx, sy, sz]
ALPHA = [np.kron(sx, SIG[d]) for d in range(3)]          # Dirac 4 composantes : alpha_d = sigma_x (x) sigma_d, beta = sigma_z (x) 1
BETA = np.kron(sz, s0)
I4 = np.eye(4)


def V_face_W(d, sign):
    """demi-tour autour de l'axe d sur les spineurs (analogue de D2.V_face = Y_a Y_b sur l'octet)"""
    a, b = [c for c in range(3) if c != d]
    return sign * np.kron(s0, SIG[a] @ SIG[b])


def neighbour_W(mode, N, q, co, d, e, signs):
    q2, co2, _ = D.neighbour(mode, N, q, co, d, e, signs)
    if mode != 'oreiller' or q2 == q:
        return q2, co2, I4
    face = (d, e) if q == 'A' else (d, -e)
    V = V_face_W(d, signs[face])
    return q2, co2, (V.conj().T if q == 'A' else V)


def build_W(N, mode='oreiller', signs=None, m=0.0, r=1.0):
    """Wilson-Dirac (4 comp.) sur le niveau fin : cubes A, B (oreiller) ou 8 cubes periodiques (tore)"""
    if signs is None:
        signs = {(d, e): 1 for d in range(3) for e in (1, -1)}
    cubes = ['A', 'B'] if mode == 'oreiller' else list(itertools.product((0, 1), repeat=3))
    L = 2 ** N
    cells = [(q, co) for q in cubes for co in itertools.product(range(L), repeat=3)]
    idx = {c: i for i, c in enumerate(cells)}
    blocks = {}

    def add(u, v, B):
        blocks[(u, v)] = blocks.get((u, v), 0) + B
    for (q, co) in cells:
        u = idx[(q, co)]
        add(u, u, (m + 3 * r) * BETA)
        for d in range(3):
            q2, co2, W = neighbour_W(mode, N, q, co, d, +1, signs)
            v = idx[(q2, co2)]
            B = (ALPHA[d] / 2j - r / 2 * BETA) @ W
            add(u, v, B); add(v, u, B.conj().T)
    n = 4 * len(cells)
    I, J, Vv = [], [], []
    for (u, v), B in blocks.items():
        rr, cc = np.nonzero(np.abs(B) > 1e-14)
        I += list(4 * u + rr); J += list(4 * v + cc); Vv += list(B[rr, cc])
    if HAVE_SP:
        return sp.csr_matrix((Vv, (I, J)), shape=(n, n))
    H = np.zeros((n, n), complex); np.add.at(H, (np.array(I), np.array(J)), np.array(Vv)); return H


def dense(H):
    return H.toarray() if hasattr(H, 'toarray') else np.asarray(H)


def S_rot_W(a, b, th):
    return np.cos(th / 2) * I4 - np.sin(th / 2) * np.kron(s0, SIG[a] @ SIG[b])


def edge_holonomies_W(signs):
    """meme boucle que D2.edge_holonomies, transports de spineurs"""
    N = 1; m = 2; out = {}
    for d in range(3):
        a, b = [c for c in range(3) if c != d]
        for ea, eb in itertools.product((1, -1), repeat=2):
            co = [0, 0, 0]; co[a] = m - 1 if ea > 0 else 0; co[b] = m - 1 if eb > 0 else 0
            P = I4.copy(); q, c = 'A', tuple(co)
            q, c, W = neighbour_W('oreiller', N, q, c, a, ea, signs); P = W.conj().T @ P
            q, c, W = neighbour_W('oreiller', N, q, c, b, -eb, signs); P = W.conj().T @ P
            Sp = S_rot_W(a, b, ea * eb * np.pi)
            out[(d, ea, eb)] = +1 if np.allclose(P, Sp) else (-1 if np.allclose(P, -Sp) else 0)
    return out


def classes_W():
    faces = sorted({(d, e) for d in range(3) for e in (1, -1)}); ok = []
    for bits in itertools.product((1, -1), repeat=5):
        sg = dict(zip(faces, (1,) + bits))
        if all(v == +1 for v in edge_holonomies_W(sg).values()):
            ok.append(sg)
    return ok


def ebulk_D2(N):
    for f in ('D2_N%d_v2.json' % N, 'D2_N%d.json' % N):
        p = os.path.join(D2DIR, f)
        if os.path.exists(p):
            d = json.load(open(p)).get('mu=1 s=+1', {})
            if 'E_bulk' in d:
                return float(d['E_bulk'])
    return 0.7742                                             # D2_N4_v2.log (tore, mu = 1, s = +1)


def classe_D2(N):
    for f in ('D2_N%d_v2.json' % N, 'D2_N%d.json' % N, 'D2_N3.json'):
        p = os.path.join(D2DIR, f)
        if os.path.exists(p):
            s = json.load(open(p))['signes']
            return {(int(k[0]), int(k[1:])): v for k, v in s.items()}
    raise SystemExit('classe de signes de D2 introuvable')


def niveaux(E, tol=1e-4):
    """regroupe des |E| tries en niveaux (valeur moyenne, multiplicite)"""
    E = np.sort(np.abs(E)); out = []
    for x in E:
        if out and abs(x - out[-1][0]) <= tol * max(x, 1e-12):
            v, k = out[-1]; out[-1] = ((v * k + x) / (k + 1), k + 1)
        else:
            out.append((x, 1))
    return out


def low(H, k, sigma=1e-4):
    n = H.shape[0]
    if not HAVE_SP or n <= 3000:
        E, U = np.linalg.eigh(dense(H)); o = np.argsort(np.abs(E))[:k]; return E[o], U[:, o]
    lu = sla.splu((H - sigma * sp.identity(n, dtype=complex, format='csc')).tocsc())
    op = sla.LinearOperator((n, n), matvec=lu.solve, dtype=complex)
    t1 = time.time(); print('      inversion spectrale (eigsh, hermitien) : n = %d, k = %d ...' % (n, min(k, n - 2)), flush=True)
    E, U = sla.eigsh(H.tocsc(), k=min(k, n - 2), sigma=sigma, OPinv=op, which='LM')
    print('      ... fait en %.0f s' % (time.time() - t1), flush=True)
    E = E.real; o = np.argsort(np.abs(E)); return E[o], U[:, o]


def valider():
    print('V1 : Wilson sur le tore (8 cubes periodiques) contre le spectre analytique', flush=True)
    for N in (1, 2):
        E = np.linalg.eigvalsh(dense(build_W(N, 'tore')))
        L = 2 * 2 ** N; ks = 2 * np.pi * np.arange(L) / L; ref = []
        for kx, ky, kz in itertools.product(ks, repeat=3):
            k = np.array([kx, ky, kz]); e = np.sqrt(np.sum(np.sin(k) ** 2) + np.sum(1 - np.cos(k)) ** 2); ref += [e, e, -e, -e]
        print('   N = %d : ecart max %.1e' % (N, np.max(np.abs(np.sort(E) - np.sort(ref)))))
    print('V3 : classes de signes satisfaisant le principe pour les spineurs (sigma_a sigma_b)', flush=True)
    cw = classes_W(); cd = classe_D2(3)
    print('   %d classe(s) ; classe de D2 parmi elles : %s' % (len(cw), any(c == cd for c in cw)))
    print('V2 : hermiticite, symetrie +-, jauge (signe global de B)', flush=True)
    for N in (1, 2):
        H = dense(build_W(N, 'oreiller', cd)); E = np.linalg.eigvalsh(H)
        cd2 = {k: -v for k, v in cd.items()}; E2 = np.linalg.eigvalsh(dense(build_W(N, 'oreiller', cd2)))
        Hn = dense(D.build(N, 'oreiller', s=+1, mu=1.0, signs=cd)[0])
        print('   N = %d : herm %.1e ; +- %.1e ; jauge %.1e ; natif herm %.1e ; plus petites |E| Wilson %s' % (
            N, np.max(np.abs(H - H.conj().T)), np.max(np.abs(np.sort(E) - np.sort(-E))), np.max(np.abs(np.sort(E) - np.sort(E2))),
            np.max(np.abs(Hn - Hn.conj().T)), np.round(np.sort(np.abs(E))[:8], 4)))


def mesure(N, K=160):
    t0 = time.time(); cd = classe_D2(N); res = {'N': N, 'K': K}
    print('P5 : N = %d (L = %d), classe de signes %s' % (N, 2 ** N, cd), flush=True)
    # reference Wilson
    HW = build_W(N, 'oreiller', cd)
    EW, _ = low(HW, 2 * K); EW = np.abs(EW)
    print('   Wilson : dimension %d, %d modes, |E| max %.4f (%.0f s)' % (HW.shape[0], len(EW), EW.max(), time.time() - t0), flush=True)
    # natif
    Hn, cells = D.build(N, 'oreiller', s=+1, mu=1.0, signs=cd)
    lev = np.array([c[0] for c in cells]); fine = np.repeat(lev == N, 8); root = np.repeat(lev == 0, 8)
    En, Un = low(Hn, K)
    wN = (np.abs(Un[fine]) ** 2).sum(0); w0 = (np.abs(Un[root]) ** 2).sum(0)
    keep = (wN >= 0.5) & (w0 < 0.5); Es = np.abs(En[keep])
    print('   natif : dimension %d, %d modes dont %d survivants du niveau fin, |E| max %.4f (%.0f s)' % (
        Hn.shape[0], len(En), keep.sum(), np.abs(En).max(), time.time() - t0), flush=True)
    # comparaison : seulement sous la plus petite des deux bornes atteintes
    Eb = ebulk_D2(N); cut = min(EW.max(), np.abs(En).max() * 0.98, 0.75 * Eb)
    print('   fenetre de comparaison : |E| < %.4f (0,75 E_bulk = %.4f ; bornes atteintes Wilson %.4f, natif %.4f)' % (cut, 0.75 * Eb, EW.max(), np.abs(En).max()))
    lw = [x for x in niveaux(EW) if x[0] < cut]; ln = [x for x in niveaux(Es) if x[0] < cut]
    print('   niveaux Wilson (|E|, mult.) :', [(round(a, 5), k) for a, k in lw[:10]])
    print('   niveaux natif  (|E|, mult.) :', [(round(a, 5), k) for a, k in ln[:10]])
    nmax = min(6, len(lw), len(ln))
    if nmax >= 2:
        rw = np.array([lw[i][0] / lw[0][0] for i in range(nmax)]); rn = np.array([ln[i][0] / ln[0][0] for i in range(nmax)])
        dev = (rn - rw) / rw
        A = np.vstack([np.array([lw[i][0] ** 2 for i in range(nmax)]), np.ones(nmax)]).T
        (v2, c), *_ = np.linalg.lstsq(A, np.array([ln[i][0] ** 2 for i in range(nmax)]), rcond=None)
        print('   rapports rho_n Wilson :', np.round(rw, 4)); print('   rapports rho_n natif  :', np.round(rn, 4))
        print('   ecarts relatifs       :', np.round(dev, 4), ' ; max |ecart| = %.4f' % np.max(np.abs(dev[1:])))
        print('   ajustement E_nat^2 = v^2 E_W^2 + c : v = %.4f, xi_eff = c / (v^2 E_W0^2) = %.4f' % (np.sqrt(max(v2, 0)), c / (v2 * lw[0][0] ** 2)))
        print('   multiplicites natif / Wilson :', [round(ln[i][1] / lw[i][1], 3) for i in range(nmax)])
        res.update(rho_W=rw.tolist(), rho_nat=rn.tolist(), ecart=dev.tolist(), v=float(np.sqrt(max(v2, 0))), xi_eff=float(c / (v2 * lw[0][0] ** 2)))
    ts = np.linspace(5, 60, 12); zW = [float(np.exp(-t * EW ** 2).sum()) for t in ts]; zn = [float(np.exp(-t * Es ** 2).sum()) for t in ts]
    trunc = [float(np.exp(-t * cut ** 2)) for t in ts]
    print('   trace de la chaleur (t, z_natif, z_Wilson, rapport, borne de troncature relative) :')
    for t, a, b, tr in zip(ts, zn, zW, trunc):
        print('      %5.1f  %10.4f  %10.4f  %7.4f  %.1e' % (t, a, b, a / b if b else float('nan'), tr))
    res.update(niveaux_W=lw, niveaux_nat=ln, E_W=EW.tolist(), E_nat_survivants=Es.tolist(), wN=wN.tolist(), w0=w0.tolist(), E_nat=np.abs(En).tolist(),
               t=ts.tolist(), z_nat=zn, z_W=zW, duree_s=time.time() - t0)
    json.dump(res, open(os.path.join(HERE, 'P5_N%d.json' % N), 'w'))
    print('FIN (%.0f s)' % (time.time() - t0), flush=True)


def wilson_seul(N, K=200):
    t0 = time.time(); cd = classe_D2(4); HW = build_W(N, 'oreiller', cd)
    EW, _ = low(HW, K); EW = np.abs(EW); lw = niveaux(EW)
    lw = [x for x in lw if x[0] < EW.max() * 0.98]
    print('Wilson seul, N = %d (L = %d) : dimension %d, %d modes (%.0f s)' % (N, 2 ** N, HW.shape[0], len(EW), time.time() - t0))
    print('   niveaux (|E|, mult.) :', [(round(a, 5), k) for a, k in lw[:10]])
    rw = [a / lw[0][0] for a, _ in lw[:8]]; print('   rapports rho_n :', np.round(rw, 4))
    json.dump(dict(N=N, niveaux=lw, rho=rw, E=EW.tolist()), open(os.path.join(HERE, 'P5_W%d.json' % N), 'w'))
    print('FIN', flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'valider':
        valider()
    elif sys.argv[1] == 'wilson':
        wilson_seul(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 200)
    else:
        mesure(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 160)
