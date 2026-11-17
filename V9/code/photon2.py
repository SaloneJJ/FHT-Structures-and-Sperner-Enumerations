import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('photon.py').read().split("SR = inv_under(thetaR)")[0])
import numpy as np
for k in range(3):
    D4 = stab_pairing(k)
    SL = inv_under([left(gel(g)) for g in D4])
    tau = [g for g in D4 if sgn(g) == -1 and sum(g[i] != i for i in range(4)) == 2][0]
    res = []
    for eps in (+1, -1):
        X = left(1j * (e1 - es) + eps * gprod(e2, u(tau)))
        SCF = real_dim(B1, extra_real=[lambda B: B @ X - X @ B])
        res.append(same_space(SL, SCF))
    # does SL also commute with D4 on the right ?
    errR = max(np.abs(D @ Jop(left(gel(g))) - Jop(left(gel(g))) @ D).max() for D in SL for g in D4)
    print(f"p{k}: D4 gauche -> dim {len(SL)} ; meme espace que C_F (eps=+1, -1) : {res} ; commute aussi avec D4 droit : {errR:.1e}")
import itertools
def closure(gs):
    H = {G[0]}; fr = [G[0]]
    while fr:
        h = fr.pop()
        for g in gs:
            x = mul(h, g)
            if x not in H: H.add(x); fr.append(x)
    return frozenset(H)
subs = sorted({closure([a, b]) for a in G for b in G}, key=len)
def name(H):
    o = len(H)
    cyc = {g: sorted(set(len(c) for c in [[i] for i in range(4)])) for g in H}
    orders = sorted(next(n for n in range(1, 25) if all(p == q for p, q in zip(*[(lambda x: x)(tuple(range(4))), pw(g, n)]))) for g in H)
    return o, tuple(orders)
def pw(g, n):
    x = tuple(range(4))
    for _ in range(n): x = mul(g, x)
    return x
table = {}
for H in subs:
    S = inv_under([left(gel(g)) for g in H])
    o = len(H)
    ords = tuple(sorted(next(n for n in range(1, 13) if pw(g, n) == tuple(range(4))) for g in H))
    fix_pairings = sum(1 for kk in range(3) if all(s3image(g)[kk] == kk for g in H))
    key = (o, ords, fix_pairings)
    table.setdefault(key, [0, len(S)]); table[key][0] += 1
    assert table[key][1] == len(S)
print("H (sous-groupe gauche) : ordre, ordres des elements, appariements fixes -> dim_R des D admissibles invariants")
for (o, ords, fp), (n, dmn) in sorted(table.items()):
    print(f"   |H|={o:2d} ordres={ords} app. fixes={fp} : dim {dmn:2d}  ({n} sous-groupes)")
