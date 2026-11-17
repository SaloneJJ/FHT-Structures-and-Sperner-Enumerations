# (c2) m_bb at the 8 corners of the 3-sigma box, NuFIT 6.0 "IC24 with SK-atm" (arXiv:2410.05380v2, Table 1)
from itertools import product
from math import sqrt
IC24 = {'NO': dict(d21=(6.92e-5, 8.05e-5), d3l=(2.451e-3, 2.578e-3), s13=(0.02030, 0.02388)),
        'IO': dict(d21=(6.92e-5, 8.05e-5), d3l=(2.421e-3, 2.547e-3), s13=(0.02060, 0.02409))}
IC19 = {'NO': dict(d21=(6.92e-5, 8.05e-5), d3l=(2.463e-3, 2.606e-3), s13=(0.02023, 0.02376)),
        'IO': dict(d21=(6.92e-5, 8.05e-5), d3l=(2.438e-3, 2.584e-3), s13=(0.02053, 0.02397))}
BF24 = {'NO': (7.49e-5, 2.513e-3, 0.02215), 'IO': (7.49e-5, 2.484e-3, 0.02231)}

def masses(order, m0, d21, d3l):
    if order == 'NO':
        m1 = m0; m2 = sqrt(m0*m0 + d21); m3 = sqrt(m0*m0 + d3l)
    else:  # d3l = |Dm2_32|
        m3 = m0; m2 = sqrt(m0*m0 + d3l); m1 = sqrt(m2*m2 - d21)
    return m1, m2, m3

def mbb(order, m0, d21, d3l, s13, e2, e3, ue1=2/3):
    m1, m2, m3 = masses(order, m0, d21, d3l)
    return abs(ue1*m1 + e2*(1-ue1-s13)*m2 + e3*s13*m3)*1e3  # meV

def box(order, m0, data, ue1=2/3):
    out = {}
    d = data[order]
    for e2, e3 in product((1, -1), repeat=2):
        vals = [(mbb(order, m0, a, b, c, e2, e3, ue1), (a, b, c)) for a, b, c in product(d['d21'], d['d3l'], d['s13'])]
        out[(e2, e3)] = (min(vals)[0], max(vals)[0])
    return out

for name, data in (('IC24 with SK-atm', IC24), ('IC19 w/o SK-atm', IC19)):
    print('=====', name)
    for order, m0s in (('IO', (0.0, 0.010)), ('NO', (0.0,))):
        for m0 in m0s:
            r = box(order, m0, data)
            print(f'{order}, m_lightest = {m0*1e3:.0f} meV')
            for k, (lo, hi) in r.items():
                print(f'   eta2={k[0]:+d} eta3={k[1]:+d}: {lo:7.3f} .. {hi:7.3f} meV')
            for e2 in (1, -1):
                lo = min(r[(e2, 1)][0], r[(e2, -1)][0]); hi = max(r[(e2, 1)][1], r[(e2, -1)][1])
                print(f'   band eta2={e2:+d} (both eta3): {lo:6.2f} .. {hi:6.2f} meV')
print('===== 8 corners, IC24, IO m3=0 and NO m1=0 (d21, |d3l|, s13 -> mbb for eta2=+,-)')
for order in ('IO', 'NO'):
    d = IC24[order]
    for a, b, c in product(d['d21'], d['d3l'], d['s13']):
        print(order, f'{a:.2e} {b:.3e} {c:.5f}', f'{mbb(order,0,a,b,c,1,1):7.3f} {mbb(order,0,a,b,c,-1,1):7.3f}')
print('===== central values IC24')
for order in ('IO', 'NO'):
    a, b, c = BF24[order]
    print(order, [round(mbb(order, 0, a, b, c, e2, 1), 3) for e2 in (1, -1)])
print('===== monotonicity check on a 21^3 grid (extremes at corners?)')
import itertools
for order, m0 in (('IO', 0.0), ('IO', 0.010), ('NO', 0.0)):
    d = IC24[order]
    grid = [[lo + (hi-lo)*i/20 for i in range(21)] for lo, hi in (d['d21'], d['d3l'], d['s13'])]
    for e2, e3 in product((1, -1), repeat=2):
        vals = [mbb(order, m0, a, b, c, e2, e3) for a, b, c in itertools.product(*grid)]
        cor = [mbb(order, m0, a, b, c, e2, e3) for a, b, c in product(d['d21'], d['d3l'], d['s13'])]
        print(order, m0, e2, e3, 'grid', round(min(vals), 4), round(max(vals), 4), 'corners', round(min(cor), 4), round(max(cor), 4))
print('===== sensitivity to the exactness of TM1: |U_e1|^2 = 2/3 + delta (IO, m3=0, central)')
a, b, c = BF24['IO']
for dl in (-0.02, -0.01, -0.005, -0.002, 0, 0.002, 0.005, 0.01, 0.02):
    print(f'  delta={dl:+.3f}: upper {mbb("IO",0,a,b,c,1,1,2/3+dl):6.2f}  lower {mbb("IO",0,a,b,c,-1,1,2/3+dl):6.2f}')
print('===== NO: zeros of the curves (central), m1 in meV')
a, b, c = BF24['NO']
for e2, e3 in product((1, -1), repeat=2):
    prev = None; zs = []
    for i in range(1, 40001):
        m0 = i*1e-6
        m1, m2, m3 = masses('NO', m0, a, b)
        v = 2/3*m1 + e2*(1/3-c)*m2 + e3*c*m3
        if prev is not None and v*prev < 0: zs.append(round(m0*1e3, 2))
        prev = v
    print(e2, e3, zs)
