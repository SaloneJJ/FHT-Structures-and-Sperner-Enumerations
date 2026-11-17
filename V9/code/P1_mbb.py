#!/usr/bin/env python3
"""P1 -- m_bb sous TM1 exact + mu-tau-reflexion (phases de Majorana 0 ou pi). Criteres : CRITERES_P1_mbb.md (ecrits avant calcul)."""
import itertools, json, numpy as np
D = {'NO': dict(dm21=7.49e-5, dm21s=0.19e-5, dm3=2.534e-3, dm3s=0.024e-3, s13=0.02215, s13s=0.00057),
     'IO': dict(dm21=7.49e-5, dm21s=0.19e-5, dm3=-2.510e-3, dm3s=0.025e-3, s13=0.02231, s13s=0.00056)}
def masses(o, ml, dm21, dm3):
    if o == 'NO': return ml, np.sqrt(ml**2 + dm21), np.sqrt(ml**2 + dm3)
    m2 = np.sqrt(ml**2 - dm3); return np.sqrt(m2**2 - dm21), m2, ml          # IO : m3 = ml, dm3 = dm32 < 0
def mbb(o, ml, eta, dm21, dm3, s13, ue1=None, ue2=None):
    m1, m2, m3 = masses(o, ml, dm21, dm3)
    if ue1 is None: ue1, ue2 = 2/3, 1/3 - s13
    return abs(ue1 * m1 * eta[0] + ue2 * m2 * eta[1] - s13 * m3 * eta[2])
etas = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)]
out = {}
rng = np.random.default_rng(1)
for o, p in D.items():
    out[o] = {}
    for eta in etas:
        key = ''.join('+' if e > 0 else '-' for e in eta)
        row = {}
        for ml in (0.0, 0.010, 0.030):
            c = mbb(o, ml, eta, p['dm21'], p['dm3'], p['s13'])
            sam = [mbb(o, ml, eta, p['dm21'] + 3 * p['dm21s'] * u[0], p['dm3'] + 3 * p['dm3s'] * u[1] * np.sign(p['dm3']), p['s13'] + 3 * p['s13s'] * u[2])
                   for u in rng.uniform(-1, 1, (4000, 3))]
            row[f'{ml*1e3:.0f} meV'] = (round(c * 1e3, 2), round(min(sam) * 1e3, 2), round(max(sam) * 1e3, 2))
        out[o][key] = row
        print(o, key, row)
# GJ indicatif : s12^2 in [0.30, 0.32], IO, ml = 0
for s12 in (0.30, 0.3092, 0.32):
    p = D['IO']; c13 = 1 - p['s13']
    for eta in [(1, 1, 1), (1, -1, 1)]:
        v = mbb('IO', 0.0, eta, p['dm21'], p['dm3'], p['s13'], ue1=(1 - s12) * c13, ue2=s12 * c13)
        print(f'GJ IO ml=0 s12^2={s12} eta={eta}: {v*1e3:.2f} meV')
json.dump(out, open('P1_mbb.json', 'w'), indent=1)
