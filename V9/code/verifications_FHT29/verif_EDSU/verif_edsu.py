# Vérifications EDSU (critères : CRITERES_VERIF_EDSU.md, écrits avant ce script)
import math, itertools
def dchi2(x, bf, lo1, hi1, lo3, hi3, periodic=False):
    if periodic:  # angle en degrés : on ramène x dans [lo3, lo3+360)
        while x < lo3: x += 360
        while x >= lo3 + 360: x -= 360
    if x == bf: return 0.0
    if x > bf:
        if x > hi3: return float('inf')
        s = (x-bf)/(hi1-bf) if x <= hi1 else 1 + 2*(x-hi1)/(hi3-hi1)
    else:
        if x < lo3: return float('inf')
        s = (bf-x)/(bf-lo1) if x >= lo1 else 1 + 2*(lo1-x)/(lo1-lo3)
    return s*s
def fmt(v): return '>9 (hors 3σ)' if v == float('inf') else f'{v:.2f}'
def verdict(v): return 'compatible' if v <= 4 else ('en tension' if v <= 9 else 'défavorisé')
D = {
 'NO': dict(t23=(0.470, 0.457, 0.487, 0.435, 0.585), d=(212, 171, 238, 124, 364)),
 'IO': dict(t23=(0.550, 0.535, 0.562, 0.440, 0.584), d=(274, 249, 296, 201, 335)),
}
print('=== V1 : θ23 et δ (Δχ² 1D estimés par interpolation en √Δχ²)')
for o, p in D.items():
    t = p['t23']; d = p['d']
    v05 = dchi2(0.5, *t)
    vGJ = min(dchi2(0.493 + k*0.017/100, *t) for k in range(101))
    v90 = dchi2(90, *d, periodic=True); v270 = dchi2(270, *d, periodic=True)
    # GJ : |δ| dans [89, 99] (sin²θ12 à ±1σ NuFIT 6.1), signe libre -> δ = 360-|δ| aussi
    gj = [dchi2(a, *d, periodic=True) for a in range(89, 100)] + [dchi2(360-a, *d, periodic=True) for a in range(89, 100)]
    print(f'{o}: sin²θ23=0.5 : Δχ²≈{fmt(v05)} ({verdict(v05)}) ; GJ [0.493,0.510] : min {fmt(vGJ)} ({verdict(vGJ)})')
    print(f'    δ=90° : {fmt(v90)} ({verdict(v90)}) ; δ=270° : {fmt(v270)} ({verdict(v270)}) ; GJ |δ|∈[89°,99°], signe libre : min {fmt(min(gj))} ({verdict(min(gj))}), à δ = {[*range(89,100), *[360-a for a in range(89,100)]][gj.index(min(gj))]}°')
print('    ordre : IO défavorisé, Δχ² = 4.6 (sans SK) / 9.4 (avec SK-atm, IC24) [JUNO + NuFIT 6.1]')

print('\n=== V3 : sin²θ12')
s13, ds13 = 0.02248, 0.00057
tm1 = 1 - 2/(3*(1-s13)); dtm1 = 2/3*ds13/(1-s13)**2
print(f'TM1 : sin²θ12 = {tm1:.4f} ± {dtm1:.4f} (s13² = {s13} ± {ds13})')
for nom, c, e in [('JUNO 59 j', 0.3092, 0.0087), ('NuFIT 6.1', 0.3088, 0.0067), ('Capozzi et al.', 0.3085, 0.0073), ('JUNO+NuFIT6.1 (cnf2)', 0.3096, 0.0057)]:
    n = (tm1 - c)/math.hypot(e, dtm1)
    print(f'  contre {nom:22s} {c} ± {e} : écart {n:+.2f} σ')
tab = [(0.348, 77, 80), (0.321, 87, 91), (0.309, 92, 95), (0.295, 97, 101), (0.276, 107, 109)]
def gj_delta(s):
    t = sorted(tab)
    for (a, l1, h1), (b, l2, h2) in zip(t, t[1:]):
        if a <= s <= b:
            u = (s-a)/(b-a); return l1 + u*(l2-l1), h1 + u*(h2-h1)
for s in (0.3088-0.0066, 0.3088, 0.3088+0.0067):
    l, h = gj_delta(s); print(f'  GJ : sin²θ12 = {s:.4f} -> |δ| ≈ {l:.1f}°–{h:.1f}°')

print('\n=== V2 : m_ββ en NO (m1 -> 0) et Σmν')
def box(c, lo, hi): return (lo, hi)
dm21 = (7.530e-5, 7.530e-5 - 3*0.097e-5, 7.530e-5 + 3*0.096e-5)
dm31 = (2.511e-3, 2.511e-3 - 3*0.020e-3, 2.511e-3 + 3*0.021e-3)
s13r = (s13, s13 - 3*ds13, s13 + 3*ds13)
def mbb_ii(m1, a, b, s, e2, e3):
    m2 = math.sqrt(m1**2 + a); m3 = math.sqrt(m1**2 + b)
    return abs(2/3*m1 + (1/3 - s)*m2*e2 + s*m3*e3)
for e2, e3 in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
    vals = [mbb_ii(0, a, b, s, e2, e3) for a in dm21[1:] for b in dm31[1:] for s in s13r[1:]]
    vals.append(mbb_ii(0, dm21[0], dm31[0], s13, e2, e3))
    print(f'  (ii) TM1, η2={e2:+d}, η3={e3:+d} : {1e3*min(vals):.2f}–{1e3*max(vals):.2f} meV (central {1e3*vals[-1]:.2f})')
s12r = (0.3088 - 3*0.0066, 0.3088 + 3*0.0067)
vals = []
for a in dm21[1:]:
    for b in dm31[1:]:
        for s in s13r[1:]:
            for t in s12r:
                ue2 = t*(1-s); m2 = math.sqrt(a); m3 = math.sqrt(b)
                for e2, e3 in itertools.product((1, -1), repeat=2):
                    vals.append(abs(ue2*m2*e2 + s*m3*e3))
print(f'  (i) réalité seule, m1 -> 0 : {1e3*min(vals):.2f}–{1e3*max(vals):.2f} meV')
sNO = math.sqrt(dm21[0]) + math.sqrt(dm31[0]); m2io = math.sqrt(2.484e-3); m1io = math.sqrt(2.484e-3 - dm21[0])
print(f'  Σmν min : NO {1e3*sNO:.1f} meV ; IO {1e3*(m1io+m2io):.1f} meV ; DESI DR2+CMB (ΛCDM) < 64.2 meV')

print('\n=== V2 bis (ajout déclaré après le premier calcul) : NO, m1 borné par DESI (ΛCDM)')
def sigma(m1): return m1 + math.sqrt(m1**2+dm21[0]) + math.sqrt(m1**2+dm31[0])
lo, hi = 0.0, 0.02
for _ in range(60):
    mid = (lo+hi)/2
    lo, hi = (mid, hi) if sigma(mid) < 0.0642 else (lo, mid)
m1max = lo
print(f'  Σmν < 64.2 meV  <=>  m1 < {1e3*m1max:.2f} meV (valeurs centrales)')
for e2, e3 in [(1, 1), (1, -1)]:
    for e1 in (1, -1):
        v = [abs(2/3*m1*e1 + (1/3-s13)*math.sqrt(m1**2+dm21[0])*e2 + s13*math.sqrt(m1**2+dm31[0])*e3) for m1 in [m1max*k/400 for k in range(401)]]
        print(f'  TM1, η=({e1:+d},{e2:+d},{e3:+d}) : m_ββ ∈ [{1e3*min(v):.2f}, {1e3*max(v):.2f}] meV pour m1 ∈ [0, {1e3*m1max:.2f}] meV')
