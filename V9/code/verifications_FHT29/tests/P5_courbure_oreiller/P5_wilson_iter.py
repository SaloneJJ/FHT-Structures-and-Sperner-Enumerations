#!/usr/bin/env python3
"""P5 -- contrôle de la référence Wilson sur l'oreiller sans factorisation LU (amendement 3 de CRITERES_P5.md).
Itération de sous-espace filtrée par Tchebychev sur A = H^2 (H = Wilson-Dirac de la v2, inchangé).
Usage :
  python3 P5_wilson_iter.py essai        N = 4 : validation contre P5_N4.json (bloquante) et estimation de durée pour N = 5
  python3 P5_wilson_iter.py n5           N = 5 (L = 32) : contrôle ; écrit P5_W5_iter.json
JJ SALONE, FHT 2026."""
import os, sys, json, time
import numpy as np
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import P5_courbure_oreiller_v2 as P
CUT = 0.5806
def grappes(E, tol=2e-3, cut=CUT):
    E = np.sort(np.abs(np.asarray(E))); E = E[E < cut]; out = []; cur = [E[0]]
    for x in E[1:]:
        if (x - cur[-1]) / x < tol: cur.append(x)
        else: out.append((float(np.mean(cur)), len(cur))); cur = [x]
    out.append((float(np.mean(cur)), len(cur))); return out
def chfsi(H, nwant=120, b=160, deg=40, tol=1e-7, maxit=60, seed=1):
    n = H.shape[0]; H = H.tocsr()
    A = lambda X: H @ (H @ X)
    gb = n * b * 16 * 5 / 1e9
    print(f'   dimension {n}, bloc {b}, mémoire estimée ≈ {gb:.1f} Go', flush=True)
    if gb > 8: raise SystemExit('mémoire estimée > 8 Go : arrêt')
    rng = np.random.default_rng(seed)
    v = rng.normal(size=n) + 1j * rng.normal(size=n); v /= np.linalg.norm(v)
    for _ in range(60): w = A(v); lam = np.linalg.norm(w); v = w / lam
    lmax = 1.10 * lam
    X = rng.normal(size=(n, b)) + 1j * rng.normal(size=(n, b)); X, _ = np.linalg.qr(X)
    a = lmax / 40; t0 = time.time()
    for it in range(1, maxit + 1):
        e = (lmax - a) / 2; c = (lmax + a) / 2
        Y0 = X; Y1 = (A(X) - c * X) / e
        for _ in range(2, deg + 1):
            Y2 = 2 * (A(Y1) - c * Y1) / e - Y0; Y0, Y1 = Y1, Y2
            s = np.linalg.norm(Y1, axis=0).max(); Y0 = Y0 / s; Y1 = Y1 / s
        Q, _ = np.linalg.qr(Y1); AQ = A(Q); G = Q.conj().T @ AQ; G = (G + G.conj().T) / 2
        th, W = np.linalg.eigh(G); X = Q @ W; AX = AQ @ W
        R = np.linalg.norm(AX[:, :nwant] - X[:, :nwant] * th[:nwant], axis=0)
        nconv = int(np.sum(R < tol)); a = max(th[-1], th[nwant] * 1.0001)
        print(f'   itération {it:2d} : convergées {nconv}/{nwant}, résidu max {R.max():.1e}, |E| voulue max {np.sqrt(max(th[nwant-1],0)):.4f}, '
              f'borne a -> |E| {np.sqrt(a):.4f}, λmax {lmax:.1f}, {time.time() - t0:.0f} s', flush=True)
        if nconv >= nwant: break
    return np.sqrt(np.clip(th[:nwant], 0, None)), time.time() - t0, it
def essai():
    cd = P.classe_D2(4); H = P.build_W(4, 'oreiller', cd)
    E, dt, it = chfsi(H)
    ref = np.sort(np.abs(json.load(open(os.path.join(HERE, 'P5_N4.json')))['E_W']))
    ref = ref[ref < CUT]; m = min(len(ref), len(E)); E = np.sort(E)
    err = np.max(np.abs(E[:m] - ref[:m]))
    print(f'ESSAI N = 4 : {m} valeurs comparées, écart max {err:.1e} -> {"VALIDÉ" if err < 1e-6 and m >= 112 else "ÉCHEC (ne pas lancer n5)"}')
    print(f'   durée {dt:.0f} s en {it} itérations ; estimation grossière pour N = 5 : {dt * 8 * 2.5 / 60:.0f} min (×8 en taille, ×2,5 en itérations), à prendre avec prudence')
def n5():
    cd = P.classe_D2(4); H = P.build_W(5, 'oreiller', cd)
    E, dt, it = chfsi(H)
    g32 = grappes(E, cut=E.max() * 0.999); d16 = json.load(open(os.path.join(HERE, 'P5_N4.json'))); g16 = grappes(d16['E_W'])
    print('   L = 32 :', [(round(a, 5), k) for a, k in g32[:12]]); print('   L = 16 :', [(round(a, 5), k) for a, k in g16[:12]])
    P_ = []
    for (a16, k16), (a32, k32) in zip(g16, g32):
        if k16 != k32: print(f'   multiplicités différentes ({k16} contre {k32}) : arrêt'); break
        P_.append((a16, a32))
    res = dict(E=E.tolist(), niveaux_32=g32, niveaux_16=g16, duree_s=dt, iterations=it)
    if len(P_) >= 3:
        P_ = np.array(P_[:10]); r16 = P_[:, 0] / P_[0, 0]; r32 = P_[:, 1] / P_[0, 1]; dv = (r32 - r16) / r16
        print(f'   ρ L=16 {np.round(r16, 4)}\n   ρ L=32 {np.round(r32, 4)}\n   variation {np.round(dv, 4)} ; max {np.max(np.abs(dv)):.4f} -> '
              f'{"CONTRÔLE PASSÉ" if np.max(np.abs(dv)) < 0.03 else "CONTRÔLE ÉCHOUÉ"}')
        res.update(rho16=r16.tolist(), rho32=r32.tolist(), variation=dv.tolist())
    else: print('   moins de 3 niveaux appariés : contrôle non concluant')
    json.dump(res, open(os.path.join(HERE, 'P5_W5_iter.json'), 'w')); print(f'FIN ({dt:.0f} s)', flush=True)
if __name__ == '__main__':
    {'essai': essai, 'n5': n5}[sys.argv[1]]()
