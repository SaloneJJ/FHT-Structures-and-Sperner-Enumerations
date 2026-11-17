"""Analyse des balayages N = 3 : convergence en pas de theta, demi-grilles, exclusion de theta = 0 mod 2pi/q, fenetre eps/2,
structure pres de theta = 0. Usage : python3 rv3_an_N3.py tag   (balayage sur [0, pi/q], symetries f(-th) = f(th), periode 2pi/q)"""
import sys, numpy as np
sys.path.insert(0, '/home/claude/fht/relecture3'); from rv3_an import integrand
z = np.load(f'/home/claude/fht/relecture3/scan_{sys.argv[1]}.npz'); ths = z['ths']; q = float(z['q']); N = int(z['N']); eps = float(z['eps'])
P = 2 * np.pi / q; nper = int(round(q))             # periode ; le balayage couvre [0, P/2]
h = ths[1] - ths[0]; M = len(ths) - 1                # M intervalles sur [0, P/2]
fmt = lambda F: '[' + ' '.join(f'{x:+.3f}' for x in F) + f'] total {F.sum():+.3f}'
def full_grid_sum(f, NTH, phase=0):
    """somme de Riemann de V166 sur la grille 2 pi j / NTH (j = 0..NTH-1), par symetrie et periodicite, a partir des points du balayage"""
    tot = np.zeros(f.shape[1]); dth = 2 * np.pi / NTH
    for j in range(NTH):
        t = (2 * np.pi * j / NTH) % P; t = min(t, P - t)        # reduit dans [0, P/2]
        k = t / h
        if abs(k - round(k)) > 1e-6: return None
        tot += f[int(round(k))] * dth
    return tot
for e_ in (eps, eps / 2):
    f = integrand(z, e_)
    w = np.full(M + 1, h); w[0] = w[-1] = h / 2
    exact = 2 * nper * (w[:, None] * f).sum(0)       # trapeze sur [0, P/2], x 2 (parite) x q (periodes)
    print(f'=== {sys.argv[1]} : N = {N}, q = {q}, fenetre {e_:.4f}, balayage {M + 1} points sur [0, {ths[-1]:.4f}] (pas {h:.5f}) ===')
    print(f'  integrale (trapeze, equivaut a {int(round(2 * np.pi / h))} pas) : {fmt(exact)}')
    for NTH in (32, 48, 64, 128, 256, 512, 1024):
        s = full_grid_sum(f, NTH)
        if s is not None: print(f'  grille V166 a {NTH:4d} pas : {fmt(s)}')
    if e_ == eps:
        for NTH in (64, 128):
            d = 2 * np.pi / NTH; js = np.arange(NTH); th = 2 * np.pi * js / NTH
            red = lambda t: np.minimum(t % P, P - t % P)
            idx = np.round(red(th) / h).astype(int)
            if np.any(np.abs(red(th) / h - idx) > 1e-6): continue
            ev, od = js % 2 == 0, js % 2 == 1; zero = np.isclose(red(th), 0)
            print(f'  {NTH} pas : pairs {fmt(f[idx[ev]].sum(0) * d)} ; impairs {fmt(f[idx[od]].sum(0) * d)}')
            print(f'  {NTH} pas : points theta = 0 mod 2pi/q seuls {fmt(f[idx[zero]].sum(0) * d)} ; complete hors ces points {fmt(f[idx[~zero]].sum(0) * d)}')
        # structure pres de 0
        E, S, W = z['E'], z['S'], z['W']
        print('  pres de theta = 0 : les 4 niveaux de plus petit |E| (E, pente, poids couche fine)')
        for i in list(range(0, 9)) + [12, 16, 24, 32]:
            if i >= len(ths): break
            o = np.argsort(np.abs(E[i]))[:4]
            print(f'   th = {ths[i]:.4f} : ' + ' ; '.join(f'{E[i, k]:+.4f}/{S[i, k]:+.3f}/{W[i, N, k]:.2f}' for k in o) + f' ; f_fine = {f[i, N]:+.3f}, f_tot = {f[i].sum():+.3f}')
        mins = np.nanmin(np.abs(E), axis=1); i0 = np.argmin(mins)
        print(f'  min |E| sur le balayage : {mins[i0]:.2e} a theta = {ths[i0]:.4f} ; |pente| max dans la fenetre : {np.nanmax(np.where(np.abs(E) < eps, np.abs(S), 0)):.3f}')
