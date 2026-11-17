# Amendement 2 de CRITERES_P5.md : appariement par multiplicité
import json, sys, numpy as np
CUT = 0.5806
def grappes(E, tol):
    E = np.sort(np.abs(np.asarray(E))); E = E[E < CUT]; out = []; cur = [E[0]]
    for x in E[1:]:
        if (x - cur[-1]) / x < tol: cur.append(x)
        else: out.append((float(np.mean(cur)), len(cur))); cur = [x]
    out.append((float(np.mean(cur)), len(cur))); return out
def lire(fn, tol):
    d = json.load(open(fn)); W = grappes(d['E_W'], tol); Nn = grappes(d['E_nat_survivants'], tol)
    print(f'--- {fn}, seuil de fusion {tol:g}')
    print('  Wilson :', [(round(a, 5), k) for a, k in W]); print('  natif  :', [(round(a, 5), k) for a, k in Nn])
    P = []
    for (ew, kw), (en, kn) in zip(W, Nn):
        if 2 * kn != kw: print(f'  paire ({ew:.5f}, ×{kw}) / ({en:.5f}, ×{kn}) : rapport {kn/kw:.3f} ≠ 1/2, arrêt'); break
        P.append((ew, en))
    if len(P) < 3: print('  moins de 3 paires valides : NON RÉSOLU'); return
    P = np.array(P); rw = P[:, 0] / P[0, 0]; rn = P[:, 1] / P[0, 1]; dev = (rn - rw) / rw
    A = np.column_stack([P[:, 0] ** 2, np.ones(len(P))]); (v2, c), *_ = np.linalg.lstsq(A, P[:, 1] ** 2, rcond=None)
    m = np.max(np.abs(dev[1:]))
    verdict = 'STANDARD' if m <= 0.05 else ('ANOMAL' if m >= 0.10 else 'NON RÉSOLU')
    print(f'  {len(P)} paires valides ; ρ Wilson {np.round(rw, 4)} ; ρ natif {np.round(rn, 4)}')
    print(f'  écarts {np.round(dev, 4)} ; max |écart| = {m:.4f} -> {verdict} (sous réserve du contrôle Wilson L = 32)')
    print(f'  v = {np.sqrt(v2):.4f} ; ξ_eff = {c / (v2 * P[0, 0] ** 2):+.4f}')
for tol in (2e-3, 1e-3, 5e-3): lire('P5_N4.json', tol)
