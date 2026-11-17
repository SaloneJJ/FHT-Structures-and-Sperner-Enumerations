# Lecture exacte des tables NuFIT 6.1 (critères : CRITERES_VERIF_EDSU.md + addendum)
import sys, os, math
import numpy as np
D = sys.argv[1]  # dossier des tables décompressées
WANT = ['T23/DCP', 'T12/DCP', 'T12', 'T23', 'DCP', 'T13']
def read(path):
    secs, cur, order = {}, None, None
    with open(path) as fh:
        for line in fh:
            if line.startswith('###'): order = line.strip('# \n'); continue
            if line.startswith('#'):
                tag = line[1:].split('projection')[0].strip()
                cur = tag if tag in WANT else None
                if cur: secs[cur] = []
                continue
            if cur and line.strip(): secs[cur].append([float(x) for x in line.split()])
    return order, {k: np.array(v) for k, v in secs.items()}
F = {}
for f in sorted(os.listdir(D)):
    o, s = read(os.path.join(D, f)); F[f] = (o, s)
    t23 = s['T23']; i = t23[:, 1].argmin()
    print(f'{f}: {o} ; min Δχ²(T23) = {t23[i,1]:.3f} en sin²θ23 = {t23[i,0]:.3f}')
def interp1(tab, x):
    return float(np.interp(x, tab[:, 0], tab[:, 1]))
def grid2(tab):
    xs = np.unique(tab[:, 0]); ys = np.unique(tab[:, 1])
    Z = np.full((len(xs), len(ys)), np.nan)
    ix = {v: k for k, v in enumerate(xs)}; iy = {v: k for k, v in enumerate(ys)}
    for x, y, z in tab[:, :3]: Z[ix[x], iy[y]] = z
    return xs, ys, Z
def interp2(g, x, y):
    xs, ys, Z = g
    i = np.searchsorted(xs, x) - 1; j = np.searchsorted(ys, y) - 1
    i = min(max(i, 0), len(xs)-2); j = min(max(j, 0), len(ys)-2)
    u = (x-xs[i])/(xs[i+1]-xs[i]); v = (y-ys[j])/(ys[j+1]-ys[j])
    return float((1-u)*(1-v)*Z[i,j] + u*(1-v)*Z[i+1,j] + (1-u)*v*Z[i,j+1] + u*v*Z[i+1,j+1])
v1 = lambda d: 'compatible' if d <= 4 else ('en tension' if d <= 9 else 'défavorisé')
v2 = lambda d: 'compatible (≤2σ)' if d <= 6.18 else ('en tension (≤3σ)' if d <= 11.83 else 'défavorisé (>3σ)')
GJ = [(0.348, 77, 80), (0.321, 87, 91), (0.309, 92, 95), (0.295, 97, 101), (0.276, 107, 109)]
def gj_curve():
    t = sorted(GJ); pts = []
    for (a, l1, h1), (b, l2, h2) in zip(t, t[1:]):
        for k in range(51):
            u = k/50; s = a + u*(b-a)
            for w in range(6):
                d = (l1 + u*(l2-l1)) + w/5*((h1 + u*(h2-h1)) - (l1 + u*(l2-l1)))
                pts.append((s, d))
    return pts
out = {}
for f, (o, s) in F.items():
    out[f] = o
for f, (o, s) in F.items():
    # minimum global = min sur les deux fichiers de la même variante : identifié plus bas ; ici valeurs brutes
    pass
# Regroupement en variantes : paire (NO, IO) ; identification par la signature (min IO, position du min θ23 en NO)
NO = [f for f in F if F[f][0].startswith('Normal')]; IO = [f for f in F if F[f][0].startswith('Inverted')]
print('\nNO :', NO, ' IO :', IO)
for f in IO: print(f'  {f}: min Δχ² IO (T23 1D) = {F[f][1]["T23"][:,1].min():.3f}')
for f in NO: print(f'  {f}: min Δχ² NO (T23 1D) = {F[f][1]["T23"][:,1].min():.3f}, argmin sin²θ23 = {F[f][1]["T23"][F[f][1]["T23"][:,1].argmin(),0]:.3f}')
def report(label, f):
    o, s = F[f]
    m0 = s['T23'][:, 1].min()
    print(f'\n### {label} — {o} ({f}) ; offset de l\'ordre (min Δχ² dans ce fichier) = {m0:.3f}')
    d05 = interp1(s['T23'], 0.5); dgj = min(interp1(s['T23'], x) for x in np.linspace(0.493, 0.510, 35))
    print(f'  1D T23 : Δχ²(0.5) = {d05:.2f} [rel. ordre {d05-m0:.2f}] {v1(d05)} ; min GJ [0.493,0.510] = {dgj:.2f} [{dgj-m0:.2f}]')
    dp = interp1(s['DCP'], 90); dm = interp1(s['DCP'], -90)
    print(f'  1D DCP : Δχ²(+90°) = {dp:.2f} [{dp-m0:.2f}] {v1(dp)} ; Δχ²(−90°) = {dm:.2f} [{dm-m0:.2f}] {v1(dm)}')
    i = s['DCP'][:,1].argmin(); print(f'           min DCP en δ = {s["DCP"][i,0]:.0f}°')
    d12 = interp1(s['T12'], 0.318); print(f'  1D T12 : Δχ²(0.318) = {d12:.2f} [{d12-m0:.2f}] {v1(d12-m0)} (rel. ordre)')
    g = grid2(s['T23/DCP'])
    a = interp2(g, 0.5, 90); b = interp2(g, 0.5, -90)
    gjr = min(interp2(g, x, sg*d) for x in np.linspace(0.493, 0.510, 18) for d in np.linspace(89, 98, 10) for sg in (1, -1))
    print(f'  2D T23/DCP : (0.5, +90°) = {a:.2f} [{a-m0:.2f}] {v2(a)} ; (0.5, −90°) = {b:.2f} [{b-m0:.2f}] {v2(b)} ; GJ rectangle min = {gjr:.2f} [{gjr-m0:.2f}] {v2(gjr)}')
    g12 = grid2(s['T12/DCP'])
    best = min((interp2(g12, x, sg*d), x, sg*d) for x, d in gj_curve() for sg in (1, -1))
    tm1 = min(interp2(g12, 0.318, sg*90) for sg in (1, -1))
    print(f'  2D T12/DCP : courbe GJ min = {best[0]:.2f} [{best[0]-m0:.2f}] {v2(best[0])} en (sin²θ12, δ) = ({best[1]:.3f}, {best[2]:.1f}°) ; TM1 (0.318, ±90°) min = {tm1:.2f} [{tm1-m0:.2f}] {v2(tm1)}')
    return m0
import itertools
glob = {}
for no, io in itertools.product(NO, IO):
    pass
# identification par Δm²3ℓ (addendum 2 ; id_dma.py) : TBoff = (144289, 196693), TByes = (506694, 460581)
no_yes, io_yes, no_off, io_off = '506694.txt', '460581.txt', '144289.txt', '196693.txt'
print(f'\nIdentification (Δm²3ℓ) : TByes = ({no_yes}, {io_yes}) ; TBoff = ({no_off}, {io_off})')
for lab, (n, i) in [('TByes (IC24 + SK-atm)', (no_yes, io_yes)), ('TBoff (IC23, sans SK-atm)', (no_off, io_off))]:
    report(lab, n); report(lab, i)
