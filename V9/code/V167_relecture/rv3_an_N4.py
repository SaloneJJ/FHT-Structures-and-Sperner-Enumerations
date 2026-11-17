"""Analyse N = 4 : (1) reproduction des sommes de V166 (grille 64 pas, pairs/impairs, theta = 0 mod pi) a partir des balayages ;
(2) integrale continue par couche (trapeze non uniforme sur la reunion des balayages), fenetres eps, eps/2 et fenetre lisse cos^2 ;
(3) structure pres de theta = 0. Symetries utilisees (verifiees) : f(-th) = f(th), periode 2pi/q.
Usage : python3 rv3_an_N4.py q tags_uniformes(sep ,) [tags_zoom(sep ,)]"""
import sys, numpy as np
q = float(sys.argv[1]); P = 2 * np.pi / q; N = 4
def load(tag):
    z = np.load(f'/home/claude/fht/relecture3/scan_{tag}.npz'); ok = ~np.isnan(z['E'][:, 0])
    return z['ths'][ok], z['E'][ok], z['S'][ok], z['W'][ok], float(z['eps'])
def integ(E, S, W, eps, kind='box', absval=False):
    m = np.abs(E) < eps
    if absval: S = np.abs(S)
    if kind == 'box': g = np.where(m, 1 / (2 * eps), 0.0)
    else: g = np.where(m, np.cos(np.pi * E / (2 * eps)) ** 2 / eps, 0.0)
    return np.einsum('tnk,tk->tn', W, S * g)
uni = [load(t) for t in sys.argv[2].split(',')]
zoom = [load(t) for t in sys.argv[3].split(',')] if len(sys.argv) > 3 else []
eps = uni[0][4]
fmt = lambda F: '[' + ' '.join(f'{x:+.3f}' for x in F) + f'] total {F.sum():+.3f}'
# (1) grille V166
th_u = np.concatenate([u[0] for u in uni]); E_u = np.concatenate([u[1] for u in uni]); S_u = np.concatenate([u[2] for u in uni]); W_u = np.concatenate([u[3] for u in uni])
o = np.argsort(th_u); th_u, E_u, S_u, W_u = th_u[o], E_u[o], S_u[o], W_u[o]
_, iu = np.unique(np.round(th_u, 10), return_index=True); th_u, E_u, S_u, W_u = th_u[iu], E_u[iu], S_u[iu], W_u[iu]
print(f'=== N = 4, q = {q} : {len(th_u)} points uniformes sur [{th_u[0]:.4f}, {th_u[-1]:.4f}] ; zooms : {[ (round(z[0][0],4), round(z[0][-1],4), len(z[0])) for z in zoom]} ; eps = {eps:.4f} ===')
f_u = integ(E_u, S_u, W_u, eps)
def at(t):
    t = t % P; t = min(t, P - t); i = np.argmin(np.abs(th_u - t)); return i if abs(th_u[i] - t) < 1e-9 else None
for NTH in (64, 128):
    js = np.arange(NTH); idx = [at(2 * np.pi * j / NTH) for j in js]
    if any(i is None for i in idx): print(f'  grille {NTH} pas : non couverte'); continue
    idx = np.array(idx); d = 2 * np.pi / NTH; red = np.array([min((2*np.pi*j/NTH) % P, P - (2*np.pi*j/NTH) % P) for j in js]); zero = np.isclose(red, 0)
    print(f'  grille {NTH} pas (Riemann, comme V166) : complete {fmt(f_u[idx].sum(0) * d)}')
    print(f'      pairs {fmt(f_u[idx[js % 2 == 0]].sum(0) * d)} ; impairs {fmt(f_u[idx[js % 2 == 1]].sum(0) * d)}')
    print(f'      theta = 0 mod 2pi/q seuls {fmt(f_u[idx[zero]].sum(0) * d)} ; hors ces points {fmt(f_u[idx[~zero]].sum(0) * d)}')
    if NTH == 64:
        print(f'      sous-grilles 1/4 et 3/4 : {fmt(f_u[idx[js % 4 == 1]].sum(0) * d)} / {fmt(f_u[idx[js % 4 == 3]].sum(0) * d)}')
# (2) integrale continue : reunion des points
TH = [th_u]; EE = [E_u]; SS = [S_u]; WW = [W_u]
for z in zoom: TH.append(z[0]); EE.append(z[1]); SS.append(z[2]); WW.append(z[3])
TH = np.concatenate(TH); EE = np.concatenate(EE); SS = np.concatenate(SS); WW = np.concatenate(WW)
o = np.argsort(TH, kind='stable'); TH, EE, SS, WW = TH[o], EE[o], SS[o], WW[o]
_, iu = np.unique(np.round(TH, 10), return_index=True); TH, EE, SS, WW = TH[iu], EE[iu], SS[iu], WW[iu]
fac = 2 * round(P / (2 * np.pi) ** 0 * 2 * np.pi / P) if False else 2 * int(round(2 * np.pi / P))   # parite x nombre de periodes
cover = TH[-1] - TH[0]
for e_, kind, lab in ((eps, 'box', 'fenetre eps'), (eps / 2, 'box', 'fenetre eps/2'), (eps, 'cos2', 'fenetre lisse cos^2, demi-largeur eps')):
    f = integ(EE, SS, WW, e_, kind); I = fac * np.trapezoid(f, TH, axis=0)
    fu = integ(E_u, S_u, W_u, e_, kind); Iu = fac * np.trapezoid(fu, th_u, axis=0)
    print(f'  integrale continue, {lab} : {fmt(I)}   (uniforme seule, pas {th_u[1]-th_u[0]:.4f} : {fmt(Iu)})')
fb = integ(EE, SS, WW, eps, absval=True); print(f'  flot brut continu (fenetre eps) : [' + ' '.join(f'{x:.3f}' for x in fac * np.trapezoid(fb, TH, axis=0)) + ']')
fbu = integ(E_u, S_u, W_u, eps, absval=True)
for NTH in (64,):
    js = np.arange(NTH); idx = [at(2 * np.pi * j / NTH) for j in js]
    if all(i is not None for i in idx): print(f'  flot brut, grille 64 pas : [' + ' '.join(f'{x:.3f}' for x in fbu[np.array(idx)].sum(0) * 2 * np.pi / NTH) + ']')
# contribution de la zone d hybridation : |theta| < 0.01, comparee a l integrande regulier prolonge
f = integ(EE, SS, WW, eps); m = TH <= 0.01
if m.sum() > 3:
    reg = f[(TH > 0.01) & (TH < 0.03)].mean(0)
    print(f'  zone |theta| < 0,01 : integrale {fmt(fac * np.trapezoid(f[m], TH[m], axis=0))} ; avec l integrande regulier (moyenne sur 0,01-0,03) : {fmt(fac * 0.01 * reg)}')
# (3) structure
print('  theta : niveaux de plus petit |E| (E/pente/poids fin) | f_fin f_total')
for i in range(len(TH)):
    if TH[i] in (0.0,) or np.isclose(TH[i], [0.0005, 0.001, 0.002, 0.003, 0.005, 0.01, 0.05, 0.1, 0.15, 0.19, 0.2]).any() or (i % 4 == 0 and TH[i] > 0.2):
        k = np.argsort(np.abs(EE[i]))[:4]
        print(f'   {TH[i]:.4f} : ' + ' ; '.join(f'{EE[i, j]:+.4f}/{SS[i, j]:+.3f}/{WW[i, N, j]:.2f}' for j in k) + f' | {f[i, N]:+.3f} {f[i].sum():+.3f} (n fenetre {np.sum(np.abs(EE[i]) < eps)})')
print(f'  min |E| : {np.nanmin(np.abs(EE)):.2e} ; max |pente| dans la fenetre : {np.max(np.where(np.abs(EE) < eps, np.abs(SS), 0)):.3f} ; fenetres saturees (>= 15 etats) : {(np.sum(np.abs(EE) < eps, axis=1) >= 15).sum()}')
