"""Independent recomputation (A1, tab:valpha / prop:census6): for the representative group of each
of the 54 classes (v,alpha), v <= 6, read from the table of the article (generators; file valpha54_table_rows.tex, or the .tex of Article I), compute
(1) the real Wedderburn type (r,c,h) of R[G], (2) for every sign character eta of G the decomposition
of R^eta[G], and list the cases with a sector A_F (C, H, M3(C) among the simple factors).
Numerical Wedderburn decomposition (twisted_decomp.decompose). S6 (order 720) is treated by its degrees."""
import re, sys
from twisted_decomp import decompose, group_closure, compose, sign
tex = open(sys.argv[1] if len(sys.argv) > 1 else 'valpha54_table_rows.tex', encoding='utf-8').read()
i = tex.index(r'\label{tab:valpha}'); body = tex[i:i+40000]
rows = [l for l in body.split('\n') if re.match(r'^\d & ', l)][:54]
assert len(rows) == 54, len(rows)
def perm_from_cycles(txt, v):
    p = list(range(v))
    if txt.strip() in ('$e$', 'e'): return tuple(p)
    for cyc in re.findall(r'\(([\d\\, ]+)\)', txt):
        pts = [int(x) - 1 for x in re.findall(r'\d', cyc)]
        for a, b in zip(pts, pts[1:] + pts[:1]): p[a] = b
    return tuple(p)
def characters(G, gens):
    """all homomorphisms G -> {+-1}, by trying sign assignments on the generators"""
    import itertools
    idx = {g: k for k, g in enumerate(G)}; out = []
    e = tuple(range(len(G[0])))
    for signs in itertools.product([1, -1], repeat=len(gens)):
        val = {e: 1}; frontier = [e]; ok = True
        while frontier and ok:
            new = []
            for g in frontier:
                for s, eps in zip(gens, signs):
                    h = compose(s, g); w = val[g] * eps
                    if h in val:
                        if val[h] != w: ok = False; break
                    else: val[h] = w; new.append(h)
                if not ok: break
            frontier = new
        if ok and len(val) == len(G):
            key = tuple(val[g] for g in G)
            if key not in [tuple(c[g] for g in G) for c in out]: out.append(val)
    return out
def rch(parts):
    r = sum(1 for p in parts if p.endswith('(R)') or p == 'R'); c = sum(1 for p in parts if 'C' in p)
    h = sum(1 for p in parts if 'H' in p); return (r, c, h)
nAF = 0; bad = 0
for row in rows:
    cols = [c.strip() for c in row.split('&')]
    v = int(cols[0]); name = cols[1]; order = int(cols[2]); gtxt = cols[3]; rch_tab = cols[8].strip(' \\').strip('$')
    gens = [perm_from_cycles(g, v) for g in gtxt.split(r'\allowbreak')] if gtxt != '$e$' else [tuple(range(v))]
    G = group_closure(gens, v)
    assert len(G) == order, (row, len(G))
    if order == 720:
        print(f"v={v} {name:14s} |G|=720 : degrees 16,10,10,9,9,5,5,5,5,1,1 -> no degree 3, no sector A_F; R[G] (11,0,0) (S_n)")
        continue
    _, parts1 = decompose(G, lambda g: 1)
    t = rch(parts1); okt = (str(t).replace(' ', '') == rch_tab.replace(' ', ''))
    if not okt: bad += 1
    chars = characters(G, gens); found = []
    for eta in chars:
        if all(eta[g] == 1 for g in G): continue
        _, parts = decompose(G, lambda g: eta[g])
        if all(any(p == q for p in parts) for q in ('C', 'H', 'M3(C)')):
            sv = all(eta[g] == sign(g) for g in G)
            found.append(('sgn_V' if sv else 'other') + ': ' + ' + '.join(parts))
    if found: nAF += 1
    print(f"v={v} {name:14s} |G|={order:4d} (r,c,h) {t} table {rch_tab} {'OK' if okt else 'MISMATCH'}; "
          f"{len(chars)} characters; sector A_F: {found if found else 'none'}")
print(f"\n(r,c,h) mismatches: {bad}; classes with a sector A_F for some eta: {nAF}")
