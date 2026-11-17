#!/usr/bin/env python3
"""D14 (appariement neutre par subdivision) et D15 (compression : motif A0 sur la piece sans -I, sous S3).
Criteres : 08_descente_Xgrand/CRITERES_D14_D16.md (ecrits avant calcul).  Usage : python3 D14_D15.py | tee D14_D15.log"""
import sys, os, json, itertools, time
import numpy as np
H_ = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(H_, 'D11_D12.py')).read().split("t0 = time.time()")[0]; exec(src)
t0 = time.time()
# ---------------------------------------------------------------- D14
def tdir(t):
    a, b = vecs[t[0]], vecs[t[1]]; d = ((a[0] - b[0]) % 3, (a[1] - b[1]) % 3); return line(vecs.index(d))
partner = {0: 2, 2: 0, 3: 4, 4: 3}
Fsub = []
mids = {}
for f in S0:
    ts = [x for x in f if x.startswith('t')]
    if ts and len(f) == 3:
        t = ts[0]; n = int(t[1:]); v = [x for x in f if x.startswith('v')][0]; i = int(v[1:])
        if line(i) == partner[tdir(T[n])]:
            m = f'm{n}'; mids[t] = m; Fsub += [[t, m, 'c'], [m, v, 'c']]; continue
    Fsub.append(f)
Vs = sorted({x for f in Fsub for x in f})
print(f'D14 : {len(mids)} aretes subdivisees ; Betti {betti_z2(Fsub)} ; |Aut| = {aut_order(Vs, Fsub)}', flush=True)
def limg(p): return {line(i): line(p[i]) for i in range(8)}
pairing = [{0, 2}, {3, 4}]
SD = [p for p in GL if sorted(sorted({limg(p)[x] for x in b}) for b in pairing) == sorted(sorted(b) for b in pairing)]
_vp = vperm
def vperm(p):
    d = _vp(p); d.update({f'm{n}': f'm{T.index(tuple(sorted(p[i] for i in t)))}' for n, t in enumerate(T)}); return d
irrS, ehS, _, _ = tools(SD)
print('    homologie isotypique sous SD16 (degre, irreductible, mult) :', sorted(ehS(Fsub).items()), '; irr', irrS, flush=True)
B0m = [['z', 'x0', 'x1', 'x2', 'x3'], ['w', 'z'], ['w', 'x0'], ['w', 'x1'], ['w', 'x2'], ['w', 'x3']]
ren = lambda f, side: [side + ('w' if x == 'c' else x) for x in f]
G0 = ['g0', 'g1', 'g2', 'g3']; K4 = [list(e) for e in itertools.combinations(G0, 2)]
def assemble(S, Bc):
    X = [ren(a, 'L_') + ren(b, 'R_') + [g] for a in S for b in Bc for g in G0] + [ren(b, 'L_') + ren(a, 'R_') + e for b in Bc for a in S for e in K4]
    X = [list(f) for f in maximal([frozenset(f) for f in X])]; return X, sorted({x for f in X for x in f})
for name, Bc in (("X_PS'", B0m), ("X_LR'", B0m + [['z', 'n']])):
    X, V = assemble(Fsub, Bc); t = time.time()
    print(f'    {name} : sommets {len(V)}, Betti {betti_z2(X)}, |Aut| = {aut_order(V, X)} ({time.time() - t:.0f} s)', flush=True)
json.dump(Fsub, open(os.path.join(H_, 'S_PS_neutre.json'), 'w'))
# ---------------------------------------------------------------- D15
vperm = _vp
B = [p for p in GL if line(p[0]) == 0]
AGL = [p for p in B if vecs[p[vecs.index((1, 0))]][0] == 1]
ta = 't%d' % T.index(tuple(sorted(vecs.index((1, j)) for j in range(3))))
base = S0 + [[ta, 'o']]
irr3, eh3, orbits3, faces3 = tools(AGL)
C_, cl, idx, reps, Xc = char_table(AGL); e0 = [l for l in range(len(reps)) if reps[l] == tuple(range(8))][0]
nm = {}
for a in range(len(reps)):
    d = int(round(Xc[a][e0].real)); nm[a] = 'E' if d == 2 else ('1' if all(np.isclose(Xc[a][l], 1) for l in range(len(reps))) else 'sgn')
print('D15 base (S + pendant sur t_a) :', [(d, nm[a], int(m)) for (d, a), m in sorted(eh3(base).items())], '; |Aut| =', aut_order(sorted({x for f in base for x in f}), base), flush=True)
targets = {'type 0 (A0)': {(0, '1'), (1, 'E'), (2, 'sgn')}, 'type 1': {(1, '1'), (2, 'E'), (3, 'sgn')}}
orbs = orbits3(base); hits = []; tested = 0
for r in range(0, 3):
    for I in itertools.combinations(orbs, r):
        Iset = set().union(*map(set, I)) if I else set()
        Fi = [f for f in base if not (set(f) & Iset)] + [[v] for v in sorted(Iset)]
        Vi = {x for f in Fi for x in f}
        cands = [('rien', [], [])]
        for rw in range(1, len(orbs) + 1):
            for W in itertools.combinations(orbs, rw):
                Wset = set().union(*map(set, W))
                if not Wset <= Vi: continue
                Y = induced(Fi, Wset, faces3)
                if not Y: continue
                cands.append(('cone', [o[0] for o in W], [['N'] + y for y in Y]))
                cands.append(('susp', [o[0] for o in W], [['N'] + y for y in Y] + [['S'] + y for y in Y] + [['N', 'nN']]))
        for op, Wn, add in cands:
            G = Fi + add; b = betti_z2(G)
            if sum(b.values()) != 4: continue
            tested += 1; res = eh3(G)
            pat = {(d, nm[a]) for (d, a), m in res.items() for _ in range(int(m))}
            if sum(int(m) for m in res.values()) != 3: continue
            for tn, tg in targets.items():
                if pat == tg: hits.append((len({x for g in G for x in g}), tn, [o[0] for o in I], op, Wn, G))
print(f'D15 : {tested} testes ; succes : {len(hits)} ({time.time() - t0:.0f} s)', flush=True)
hits.sort(key=lambda h: h[0])
for nv, tn, I, op, Wn, G in hits[:8]:
    print(f'   {tn} : {nv} sommets ; isoles {I} ; {op} sur {Wn} ; |Aut| = {aut_order(sorted({x for g in G for x in g}), G)}', flush=True)
if hits: json.dump(hits[0][5], open(os.path.join(H_, 'A_compresse.json'), 'w'))
