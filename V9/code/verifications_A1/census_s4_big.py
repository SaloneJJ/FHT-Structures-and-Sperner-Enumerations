#!/usr/bin/env python3
"""census_s4_big.py -- census of census_s4.py for v >= 8, action by action, with checkpoints.
|Aut| is computed by the exact sequence of Theorem autrec (aut_order_tree), checked against the brute-force
count on all classes for v = 6, 7.  Each action writes census_v{v}_a{i}.json (canonical forms and |Aut|);
an interrupted run resumes from these files.
Usage: python3 census_s4_big.py 8 [actions...]   (actions: indices in the printed list; default all)"""
import sys, time, json, os
from collections import Counter
from census_s4 import faithful_actions, action_perms, inv_trees, MEMO, SUB, aut_order_tree
from fht_moves import canon
v = int(sys.argv[1]); sel = [int(x) for x in sys.argv[2:]]
V = frozenset(range(1, v + 1)); T0 = time.time()
acts = faithful_actions(v)
labs = []
for ia, chosen in enumerate(acts):
    lab = "+".join(str(24 // len(SUB[i])) for i in chosen) + " (stabilizers of orders " + ",".join(str(len(SUB[i])) for i in chosen) + ")"
    labs.append(lab)
    if sel and ia not in sel: continue
    fn = f"census_v{v}_a{ia}.json"
    if os.path.exists(fn):
        print(f"action {ia}: {lab}: done earlier ({fn})", flush=True); continue
    print(f"action {ia}: {lab}", flush=True)
    G = action_perms(chosen); MEMO.clear()
    trees = inv_trees(V, G); cl = {}
    for t in trees:
        c = canon(t)
        if c not in cl: cl[c] = aut_order_tree(t)
    json.dump({"action": lab, "labelled": len(trees), "classes": cl}, open(fn + ".tmp", "w"))
    os.replace(fn + ".tmp", fn)
    print(f"  {len(trees)} invariant labelled FHTs, {len(cl)} classes ({time.time()-T0:.0f}s)", flush=True)
found = {}
for ia in range(len(acts)):
    fn = f"census_v{v}_a{ia}.json"
    if not os.path.exists(fn): continue
    d = json.load(open(fn))
    for c, n in d["classes"].items():
        if c not in found: found[c] = (n, ia)
done = [ia for ia in range(len(acts)) if os.path.exists(f"census_v{v}_a{ia}.json")]
orders = Counter(n for n, ia in found.values()); s4 = Counter(ia for n, ia in found.values() if n == 24)
print(f"v={v}: actions processed {done} of {len(acts)}; {len(found)} classes with Aut containing S4")
print(f"  |Aut| distribution: {dict(sorted(orders.items()))}")
print(f"  Aut = S4: {orders[24]}, by action: " + "; ".join(f"{labs[i]}: {k}" for i, k in sorted(s4.items())) + f"; |Aut| = 48: {orders[48]}; |Aut| = 144: {orders[144]} ({time.time()-T0:.0f}s)", flush=True)
