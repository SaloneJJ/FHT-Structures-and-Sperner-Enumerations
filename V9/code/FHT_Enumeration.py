#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FHT_Enumeration.py  --  Fractal Hyper-Trees (FHT), CODE 1 (version 5 ; ex V0_FHT_Enum.py + V0_fht_upset.py + V0_FHT_fixcache.json): definition, enumeration, catalogue.
SALONE Jean-Jacques 2026

WHAT THIS PROGRAM DOES
  1. BRUTE FORCE : constructs every FHT up to isomorphism for 1 <= v <= 6 vertices (default), straight from
     the recursive decomposition, and checks the axioms (A1)-(A4), (A5 component-wise) on each realization.
  2. RECURRENCE  : counts the FHT by the recursive algorithm (Burnside on hosts, flat part from A003182),
     for 1 <= v <= 9; a species-equation computation gives an independent cross-check of the totals.
     At every v <= 6 the recurrence is compared with the brute force.
  3. CATALOGUE   : writes the table of Article 1 (V0_FHT_DefinitionAndEnumeration): one line per FHT,
     for 1 <= v <= 5:  F_k, edges, |V|, N (depth), |E|.   Files written in --outdir:
        FHT_Table.csv          (read by CODE 2, FHT_Invariants.py)
        FHT_Table_rows.tex     (body of the longtable of Article 1)
        FHT_Table.tex / .pdf   (standalone table; pdf only if pdflatex is available)
        FHT_Enum_results.json  (counts S, I, f by depth for v <= 9; Fix tables; growth data)

DEFINITION IMPLEMENTED (axiom A5, component-wise)
    (A5)  every component of level n >= 1 contains a component of level n-1.
    Levels are recomputed from the leaves: level(H) = 1 + max level of the children of H.
    Consequence (Proposition "recursive decomposition"): an FHT on V is the same thing
    as   - a partition of V into m >= 2 blocks, each carrying an FHT (trivial for a singleton),
         - ANY connected clutter ("host") with vertex set = the set of blocks.
Consequence for the Burnside count: the number of hosts fixed by a permutation g of the
blocks depends only on the cycle type of g on the m blocks (not on the depth of the blocks).

USAGE
    python3 FHT_Enumeration.py                      # brute force v<=6 (long: 10-12 minutes in pure Python),
                                                #   recurrence v<=9, table v<=5
    python3 FHT_Enumeration.py --brute-max 5        # quicker (about 1 s)
    python3 FHT_Enumeration.py --selfcheck [--deep] # internal consistency tests
    options: --vmax (recurrence, <=9)  --brute-max (<=6)  --table-max (<=brute-max, default 5)  --outdir
Limits (guard-rails): brute force v <= 6; recurrence v <= 9; below v = 6 nothing needs a cache; for
m = 7, 8 host points the Fix tables use V0_fht_fix.c / V0_fht_upset.py (see fix_generic) and the cache
V0_FHT_fixcache.json.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
from math import comb, factorial, prod

# ======================================================================
# Guard-rails (hard limits, also enforced inside the functions)
# ======================================================================
LIMIT_BRUTE = 6        # brute force: v <= 6 (v = 6: long run with pynauty)
LIMIT_RECURRENCE = 9   # recurrence: v <= 9 (hosts on at most 8 points; see below)
GENERIC_MAX_M = 8      # Burnside host counting: at most 8 host points (see fix_generic)
FLAT_INTERNAL_MAX = 6  # S(v) computed internally up to 6, then derived from the exact A003182
# Host sizes.  A nested FHT on v vertices has m <= v-1 blocks, hence a host on m <= v-1 points:
#   v <= 7 : m <= 6, direct enumeration in Python (a few seconds);
#   v = 8  : m = 7, needs the fixed-host numbers for the 14 non-trivial cycle types of S_7
#            (C helper V0_fht_fix.c, ~1 min) and the identity (labelled count, from Dedekind D(7));
#   v = 9  : m = 8.  The fixed-host numbers of the involutive types (2^j,1^r) -- (2,1^6), (2,2,1^4),
#            (2,2,2,1,1), (2,2,2,2) -- are NOT enumerated (they exceed 10^12 hosts): they come from the
#            Dedekind-type recurrence of V0_fht_upset.py (invariant antichains by splitting off a 2-cycle, then
#            plethystic logarithm), as is (3,1^5); the identity comes from the Dedekind number D(8); the remaining types
#            have few orbits and are counted by the C helper.  Burnside's lemma applied to the flat case,
#            sum_g Fix(g) = 8! * S(8), with S(8) derived from the exact A003182(8), is then an independent check of the whole
#            m = 8 table (verified: equality).
# The former shortcuts (Lemmas singletop / twotop) are NOT valid any more: they assumed that
# every host edge meets the set T of deepest blocks.
# Precomputed values of the m = 7, 8 tables are stored in V0_FHT_fixcache.json (recomputable).

# Flat FHTs = connected clutters (unlabeled).  OEIS A003182 (Dedekind numbers up to equivalence) is EXACT
# up to n = 9 (n = 8, 9: B. Pawelski); the unlabeled antichains are the Euler transform of the connected ones
# (two kinds of one-point components: an uncovered point and a point carrying the edge {x}), so
#     sum_{n>=0} (A003182(n) - 1) x^n = prod_{k>=1} (1 - x^k)^(-c_k),   c_1 = 2,  c_k = S(k) for k >= 2,
# and S(k) is obtained by the inverse Euler transform.  NOTE: the terms n = 8, 9 printed in OEIS A261006
# (1392195547909966848 and 789204635842035012967539870068113408) are double-precision roundings of the true values;
# the exact values are the ones derived here (S(8) is also recomputed independently by Burnside, see selfcheck).
A003182 = {0: 2, 1: 3, 2: 5, 3: 10, 4: 30, 5: 210, 6: 16353, 7: 490013148,
           8: 1392195548889993358, 9: 789204635842035040527740846300252680}


def flat_from_A003182(nmax=9):
    b = {0: 1}
    b.update({n: A003182[n] - 1 for n in range(1, nmax + 1)})
    d, c = {}, {}
    for n in range(1, nmax + 1):
        d[n] = n * b[n] - sum(d[k] * b[n - k] for k in range(1, n))
        c[n] = (d[n] - sum(j * c[j] for j in range(1, n) if n % j == 0)) // n
    return c


FLAT_S = flat_from_A003182(9)
FLAT_S[1] = 1      # S(1): the one-vertex FHT (c_1 = 2 counts an uncovered point and a point with its edge)
assert [FLAT_S[v] for v in range(2, 8)] == [1, 3, 14, 157, 15942, 489980450]
# Dedekind numbers A000372 (antichains of subsets of an n-set, n = 0..6), used in self-checks
DEDEKIND = {0: 2, 1: 3, 2: 6, 3: 20, 4: 168, 5: 7581, 6: 7828354,
            7: 2414682040998, 8: 56130437228687557907788}


class GuardError(Exception):
    """Raised when a computation would exceed a guard-rail."""


def log(*args):
    print(*args, flush=True)


# ======================================================================
# Combinatorial & Algebraic helpers
# ======================================================================
def bits(mask):
    i = 0
    while mask:
        if mask & 1:
            yield i
        mask >>= 1
        i += 1


def popcount(x):
    return bin(x).count("1")


def partitions(n, maxpart=None):
    """Integer partitions of n, parts in non-increasing order."""
    maxpart = n if maxpart is None else maxpart
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def part_counts(lam):
    r = {}
    for c in lam:
        r[c] = r.get(c, 0) + 1
    return r


def z(lam):
    """Order of the centralizer of a permutation of cycle type lam."""
    return prod(c ** k * factorial(k) for c, k in part_counts(lam).items())


def mult_partitions(F, k):
    """(lam, w): ways to choose k blocks among F classes with multiplicities lam."""
    for lam in partitions(k):
        l = len(lam)
        if l > F:
            continue
        w = prod(F - i for i in range(l))
        w //= prod(factorial(x) for x in part_counts(lam).values())
        yield lam, w


def perm_from_cycles(cycles, start, perm):
    p = start
    for c in cycles:
        for i in range(c):
            perm[p + i] = p + (i + 1) % c
        p += c
    return p


def image(mask, perm):
    out, i = 0, 0
    while mask:
        if mask & 1:
            out |= 1 << perm[i]
        mask >>= 1
        i += 1
    return out


# ======================================================================
# Homology & Betti numbers (N+1 Betti numbers for depth N)
# ======================================================================


# ======================================================================
# Host counting: generic Burnside
# A host on m points is ANY connected clutter whose vertex set is the whole set of m
# points (each edge then has >= 2 points automatically).
# ======================================================================
def connected(edges, full):
    reach, rest = edges[0], list(edges[1:])
    changed = True
    while changed and rest:
        changed, keep = False, []
        for e in rest:
            if e & reach:
                reach |= e
                changed = True
            else:
                keep.append(e)
        rest = keep
    return reach == full


_fixg_cache = {}


def _fix_python_dfs(m, cycles):
    """Number of connected clutters on m labelled points, covering them, invariant under a
    permutation of cycle type `cycles` (a partition of m)."""
    perm = {}
    perm_from_cycles(cycles, 0, perm)
    full = (1 << m) - 1
    seen, orbits = set(), []
    for s in range(1, full + 1):
        if s in seen or popcount(s) < 2:
            continue
        orb, x = [], s
        while x not in seen:
            seen.add(x)
            orb.append(x)
            x = image(x, perm)
        orbits.append(orb)
    n = len(orbits)
    union = [0] * n
    for i, o in enumerate(orbits):
        for a in o:
            union[i] |= a
    # two orbit-unions are compatible iff no member of one is comparable with a member of the other
    # (and, inside an orbit, members are pairwise incomparable)
    ok = [all((a & b) not in (a, b) for a in o for b in o if a != b) for o in orbits]
    compat = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if not any((a & b) in (a, b) for a in orbits[i] for b in orbits[j]):
                compat[i] |= 1 << j
    cand0 = sum(1 << i for i in range(n) if ok[i])

    def dfs(cand, cover, edges):
        total = 1 if (cover == full and connected(edges, full)) else 0
        while cand:
            lb = cand & -cand
            i = lb.bit_length() - 1
            cand ^= lb
            total += dfs(cand & compat[i], cover | union[i], edges + orbits[i])
        return total

    return dfs(cand0, 0, []) if n else 0


# ---- tables for m = 7, 8 : precomputed cache + C helper + OEIS deduction -------------------
import json
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_FILE = os.path.join(HERE, "catalogues", "FHT_fixcache.json")
C_SOURCE = os.path.join(HERE, "V0_fht_fix.c")
C_BINARY = os.path.join(HERE, "fht_fix")
_disk_cache = None


# Fix tables for m = 7, 8 (precomputed; formerly V0_FHT_FIXCACHE.json). The C helpers V0_fht_fix.c and
# V0_fht_upset.c are no longer available: these values are therefore embedded here.
EMBEDDED_FIX_CACHE = {"7:2,1,1,1,1,1": 2206396178, "7:2,2,1,1,1": 67774938, "7:2,2,2,1": 11996266, "7:3,1,1,1,1": 2041595, "7:3,2,1,1": 57383, "7:3,2,2": 9951, "7:3,3,1": 68093, "7:4,1,1,1": 57610, "7:4,2,1": 26150, "7:4,3": 265, "7:5,1,1": 1325, "7:5,2": 233, "7:6,1": 679, "7:7": 98, "8:2,1,1,1,1,1,1": 101627841297079440, "8:2,2,1,1,1,1": 182754897612838, "8:2,2,2,1,1": 7377622645180, "8:2,2,2,2": 2038188149802, "8:3,1,1,1,1,1": 262788280879, "8:3,2,1,1,1": 401256501, "8:3,2,2,1": 16356331, "8:3,3,1,1": 535150858, "8:3,3,2": 3605736, "8:4,1,1,1,1": 423760086, "8:4,2,1,1": 93886248, "8:4,2,2": 37832074, "8:4,3,1": 24459, "8:4,4": 3211190, "8:5,1,1,1": 522922, "8:5,2,1": 20530, "8:5,3": 829, "8:6,1,1": 141388, "8:6,2": 69954, "8:7,1": 3657, "8:8": 2356}

def _load_disk_cache():
    global _disk_cache
    if _disk_cache is None:
        _disk_cache = dict(EMBEDDED_FIX_CACHE)
        try:
            with open(CACHE_FILE, encoding="utf-8") as f:
                _disk_cache.update(json.load(f))
        except (OSError, ValueError):
            pass
    return _disk_cache


def _save_disk_cache():
    if _disk_cache is not None:
        os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
        with open(CACHE_FILE, "w", encoding="utf-8") as f:
            json.dump(_disk_cache, f, indent=0, sort_keys=True)


def labelled_connected_clutters(n):
    """Number of labelled connected covering clutters on n points (n <= 8), from Dedekind numbers.
    a(k): antichains with all edges of size >= 2 (empty one included); b(k): covering ones;
    exponential formula b = exp(c) with c = connected ones."""
    a = {k: sum((-1) ** s * comb(k, s) * (DEDEKIND[k - s] - 1) for s in range(k + 1))
         for k in range(0, n + 1)}
    b = {k: sum((-1) ** (k - j) * comb(k, j) * a[j] for j in range(k + 1)) for k in range(0, n + 1)}
    c = {}
    for k in range(1, n + 1):
        c[k] = b[k] - sum(comb(k - 1, j - 1) * c[j] * b[k - j] for j in range(1, k))
    return c[n]


def _c_helper():
    """Path of the compiled C helper (compiled on demand), or None."""
    if os.path.exists(C_BINARY):
        return C_BINARY
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if cc and os.path.exists(C_SOURCE):
        r = subprocess.run([cc, "-O2", "-o", C_BINARY, C_SOURCE])
        if r.returncode == 0:
            return C_BINARY
    return None


def _fix_c(m, cycles):
    exe = _c_helper()
    if exe is None:
        return None
    out = subprocess.run([exe, str(m)] + [str(c) for c in cycles],
                         capture_output=True, text=True)
    return int(out.stdout.strip())


_fixg_cache = {}


def fix_generic(m, cycles):
    """Number of connected clutters on m labelled points, covering them, invariant under a
    permutation of cycle type `cycles` (a partition of m).  Layers: memory cache, disk cache,
    identity by the Dedekind formula (m <= 8), direct enumeration (Python for m <= 6, C helper
    for m = 7 and for the few-orbit types at m = 8), and for the involutive types (2^j,1^r) and (3,1^5) at m = 8 the
    Dedekind-type recurrence of V0_fht_upset.py."""
    if m > GENERIC_MAX_M:
        raise GuardError(f"host counting is limited to {GENERIC_MAX_M} points (requested: {m})")
    cycles = tuple(sorted(cycles, reverse=True))
    key = (m, cycles)
    if key in _fixg_cache:
        return _fixg_cache[key]
    dk = f"{m}:{','.join(map(str, cycles))}"
    disk = _load_disk_cache()
    if dk in disk:
        res = disk[dk]
    elif m <= 6:
        res = _fix_python_dfs(m, cycles)
    elif cycles == (1,) * m:
        res = labelled_connected_clutters(m)
    elif m == 8 and ((2 in cycles and all(c <= 2 for c in cycles)) or (cycles[0] == 3 and all(c == 1 for c in cycles[1:]))):
        res = fix_by_recurrence(m, cycles, fix_generic, z, partitions, DEDEKIND)
        disk[dk] = res
        _save_disk_cache()
    else:
        res = _fix_c(m, cycles)
        if res is None:
            if m == 7:
                log(f"  (no C compiler: Python enumeration for m=7, type {cycles}, may take long)")
                res = _fix_python_dfs(m, cycles)
            else:
                raise GuardError("m = 8 needs the C helper (V0_fht_fix.c) or the precomputed cache")
        disk[dk] = res
        _save_disk_cache()
    _fixg_cache[key] = res
    return res


_H_cache = {}


def host_orbits(groups):
    """
    Number of orbits of hosts on m = sum(groups) points under the Young subgroup
    prod_i S_{groups[i]}: the points are split into groups of identical blocks.
    By Burnside, H = (1/prod k!) * sum_g Fix(g), and Fix(g) depends only on the cycle
    type of g on the m points.
    """
    groups = tuple(sorted(groups, reverse=True))
    if groups in _H_cache:
        return _H_cache[groups]
    m = sum(groups)
    order = prod(factorial(k) for k in groups)
    per_group = [[(lam, factorial(k) // z(lam)) for lam in partitions(k)] for k in groups]
    total = 0
    for combo in product(*per_group):
        cnt = prod(c for _, c in combo)
        cyc = tuple(sorted((c for lam, _ in combo for c in lam), reverse=True))
        total += cnt * fix_generic(m, cyc)
    assert total % order == 0, "Burnside: non-integer orbit count"
    res = total // order
    _H_cache[groups] = res
    return res


# ======================================================================
# Method 1: recurrence (Corollary cor:explicit)
# ======================================================================
def run_recurrence(vmax):
    if vmax > LIMIT_RECURRENCE:
        raise GuardError(f"the recurrence is limited to v <= {LIMIT_RECURRENCE}")
    F = {1: {0: 1}}
    res = {1: {"S": 1, "I": 0, "f": 1, "by_depth": {0: 1}}}
    log("[recurrence] v=1: S=1  I=0  f=1")
    for v in range(2, vmax + 1):
        t0 = time.time()
        if v <= FLAT_INTERNAL_MAX:
            S = host_orbits((v,))
            if S != FLAT_S[v]:
                log(f"  WARNING: computed S({v})={S} differs from the A003182 value ({FLAT_S[v]})")
            src = "Burnside"
        else:
            S, src = FLAT_S[v], "A003182, exact"
        Fv = {1: S}
        types = [(s, d) for s in range(1, v) for d in sorted(F[s])]

        def process(chosen):
            if sum(k for _, _, k in chosen) < 2:
                return
            L = max(d for _, d, _ in chosen)
            if L == 0:
                return          # flat case, handled above
            lists = [list(mult_partitions(F[s][d], k)) for s, d, k in chosen]
            for combo in product(*lists):
                w = prod(x[1] for x in combo)
                groups = tuple(p for lam, _ in combo for p in lam)
                Fv[L + 1] = Fv.get(L + 1, 0) + w * host_orbits(groups)

        def rec(i, remaining, chosen):
            if remaining == 0:
                process(chosen)
                return
            if i == len(types):
                return
            s, d = types[i]
            for k in range(remaining // s, -1, -1):
                rec(i + 1, remaining - k * s, chosen + [(s, d, k)] if k else chosen)

        rec(0, v, [])
        F[v] = Fv
        f = sum(Fv.values())
        res[v] = {"S": S, "I": f - S, "f": f, "by_depth": dict(sorted(Fv.items()))}
        log(f"[recurrence] v={v}: S={S} ({src})  I={f - S}  f={f}  "
            f"depth={dict(sorted(Fv.items()))}  ({time.time() - t0:.1f}s)")
    return res


# ======================================================================
# Isomorphism of realizations: pynauty if present, otherwise a pure-Python fallback
# (colour refinement + exhaustive search inside colour classes; fine for v <= 6)
# ======================================================================
try:
    from pynauty import Graph as _NautyGraph, certificate as _nauty_certificate
    HAVE_NAUTY = True
except ImportError:
    HAVE_NAUTY = False


def nauty_graph(v, edges):
    n = v + len(edges)
    adj = {i: [] for i in range(n)}
    for j, e in enumerate(edges):
        for x in bits(e):
            adj[x].append(v + j)
            adj[v + j].append(x)
    return _NautyGraph(n, directed=False, adjacency_dict=adj,
                       vertex_coloring=[set(range(v)), set(range(v, n))])


def _refine(v, E):
    """Label-invariant vertex colours (ranks) by colour refinement on the incidence structure."""
    E = list(E)
    vcol = [tuple(sorted(popcount(e) for e in E if (e >> x) & 1)) for x in range(v)]
    for _ in range(v + 1):
        ecol = [tuple(sorted(vcol[x] for x in bits(e))) for e in E]
        new = [(vcol[x], tuple(sorted(ecol[j] for j, e in enumerate(E) if (e >> x) & 1)))
               for x in range(v)]
        ranks = {c: i for i, c in enumerate(sorted(set(new)))}
        nxt = [ranks[c] for c in new]
        old_ranks = {c: i for i, c in enumerate(sorted(set(vcol)))}
        if len(set(nxt)) == len(set(vcol)):
            vcol = nxt
            break
        vcol = nxt
    ranks = {c: i for i, c in enumerate(sorted(set(vcol)))}
    return [ranks[c] for c in vcol]


def _classes(col):
    cl = {}
    for x, c in enumerate(col):
        cl.setdefault(c, []).append(x)
    return [cl[c] for c in sorted(cl)]


def canon_py(v, E):
    """Canonical form of the hypergraph (v, E): minimum over relabelings compatible with the refinement."""
    col = _refine(v, E)
    classes = _classes(col)
    best = None
    for choice in product(*[list(permutations(c)) for c in classes]):
        order = [x for blk in choice for x in blk]          # new label i <- old vertex order[i]
        pos = {old: new for new, old in enumerate(order)}
        cand = tuple(sorted(sum(1 << pos[x] for x in bits(e)) for e in E))
        if best is None or cand < best:
            best = cand
    return (tuple(sorted(col)), best)


def certificate(v, E):
    if HAVE_NAUTY:
        return _nauty_certificate(nauty_graph(v, E))
    return canon_py(v, E)




# ======================================================================
# Method 3 (cross-check): species equation   F(x) = x + Z_cc(p_k <- F(x^k))
# Summing the compact recurrence F_{L+1} = Z_cc[F_{<=L}] - Z_cc[F_{<L}] over L telescopes.
# Z_cc(p_1, p_2, ...) = sum_m sum_{lam |- m} Fix(lam)/z_lam p_lam  is the cycle index of
# connected clutters.  This organizes the count differently from run_recurrence (no profile
# enumeration) and only gives the totals f(v).
# ======================================================================
def run_species(vmax, mcap=None):
    from fractions import Fraction
    mcap = min(vmax, GENERIC_MAX_M) if mcap is None else mcap
    if vmax > LIMIT_RECURRENCE:
        raise GuardError(f"the species equation is limited to v <= {LIMIT_RECURRENCE}")

    def mul(a, b):
        out = [Fraction(0)] * (vmax + 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y and i + j <= vmax:
                        out[i + j] += x * y
        return out

    def subst(a, k):                  # a(x^k)
        out = [Fraction(0)] * (vmax + 1)
        for i, x in enumerate(a):
            if x and i * k <= vmax:
                out[i * k] += x
        return out

    f = [Fraction(0)] * (vmax + 1)
    if vmax >= 1:
        f[1] = Fraction(1)
    for _ in range(vmax):             # fixed-point iteration, one more correct coefficient per pass
        pk = {k: subst(f, k) for k in range(1, vmax + 1)}
        new = [Fraction(0)] * (vmax + 1)
        new[1] = Fraction(1)
        for m in range(2, mcap + 1):
            for lam in partitions(m):
                term = [Fraction(0)] * (vmax + 1)
                term[0] = Fraction(fix_generic(m, lam), z(lam))
                for c in lam:
                    term = mul(term, pk[c])
                new = [a + b for a, b in zip(new, term)]
        f = new
    out = {}
    for v in range(1, vmax + 1):
        assert f[v].denominator == 1, "species equation: non-integer coefficient"
        out[v] = int(f[v])
    for v in range(mcap + 1, vmax + 1):      # flat term (m = v hosts) read from FLAT_S
        out[v] += FLAT_S[v]
    return out


# ======================================================================
# Method 2: brute force (Proposition prop:recursive)
# ======================================================================
def host_clutters(m):
    """Generator of all labelled hosts on m points: connected covering clutters."""
    if m > LIMIT_BRUTE:
        raise GuardError(f"labelled clutters are enumerated on at most {LIMIT_BRUTE} points")
    full = (1 << m) - 1
    cands = [s for s in range(1, full + 1) if popcount(s) >= 2]

    def bt(i, chosen, cover):
        if chosen and cover == full and connected(chosen, full):
            yield tuple(chosen)
        for j in range(i, len(cands)):
            s = cands[j]
            if all((s & e) != s and (s & e) != e for e in chosen):
                chosen.append(s)
                yield from bt(j + 1, chosen, cover | s)
                chosen.pop()

    return bt(0, [], 0)


def edge_index_map(v):
    order = [sum(1 << x for x in cmb)
             for r in range(1, v + 1) for cmb in combinations(range(v), r)]
    return {mask: i for i, mask in enumerate(order)}


def edges_str(v, E):
    return " ".join("".join(str(x + 1) for x in bits(e)) for e in sorted(E))


def run_brute(vmax, flat_max=None):
    """
    Brute-force construction of all FHT classes.
    flat_max: largest v for which the flat partition (all blocks singletons) is enumerated.
              Default: vmax.  (Restricting it is only meant for cross-checks of the nested part.)
    """
    if vmax > LIMIT_BRUTE:
        raise GuardError(f"the brute force is limited to v <= {LIMIT_BRUTE}")
    flat_max = vmax if flat_max is None else flat_max
    fht = {1: [{"edges": (1,), "depth": 0, "key": (0,)}]}
    res = {1: {"S": 1, "I": 0, "f": 1, "by_depth": {0: 1}}}
    log("[brute] v=1: S=1  I=0  f=1")
    hosts_cache = {}
    nb = "pynauty" if HAVE_NAUTY else "pure-Python canonical form"
    log(f"[brute] canonical certificates: {nb}")
    for v in range(2, vmax + 1):
        t0 = time.time()
        idx = edge_index_map(v)
        seen, classes = {}, []
        plist = [p for p in partitions(v) if len(p) >= 2]
        for pi, p in enumerate(plist):
            m = len(p)
            flat = (m == v)
            if flat and v > flat_max:
                continue
            if m <= 5:
                if m not in hosts_cache:
                    hosts_cache[m] = list(host_clutters(m))
                host_iter = hosts_cache[m]
            else:
                host_iter = host_clutters(m)
            sizes = sorted(set(p), reverse=True)
            choice_lists = [list(combinations_with_replacement(
                range(len(fht[s])), p.count(s))) for s in sizes]
            for choice in product(*choice_lists):
                blocks, offset = [], 0
                for s, idxs in zip(sizes, choice):
                    for i in idxs:
                        c = fht[s][i]
                        mask = ((1 << s) - 1) << offset
                        edges = [] if s == 1 else [e << offset for e in c["edges"]]
                        blocks.append((mask, c["depth"], edges))
                        offset += s
                L = max(b[1] for b in blocks)
                base = [e for b in blocks for e in b[2]]
                for host in host_iter:
                    # no further condition on the host: any connected clutter is admissible
                    new = [sum(blocks[j][0] for j in bits(h)) for h in host]
                    E = tuple(sorted(base + new))
                    key = (tuple(sorted(idx[e] for e in E)) if flat
                           else (pi, choice, len(host), host))
                    cert = certificate(v, E)
                    if cert in seen:
                        cl = classes[seen[cert]]
                        assert cl["decomp"] == (p, choice), \
                            f"non-unique decomposition: {edges_str(v, E)}"
                        if key < cl["key"]:
                            cl["key"], cl["edges"] = key, E
                        continue
                    seen[cert] = len(classes)
                    classes.append({"edges": E, "depth": L + 1, "key": key,
                                    "decomp": (p, choice)})
        classes.sort(key=lambda c: (c["depth"], c["key"]))
        fht[v] = classes
        by_depth = dict(sorted(Counter(c["depth"] for c in classes).items()))
        S = by_depth.get(1, 0)
        res[v] = {"S": S, "I": len(classes) - S, "f": len(classes), "by_depth": by_depth}
        note = "" if v <= flat_max else "  [flat part skipped]"
        log(f"[brute] v={v}: S={S}  I={len(classes) - S}  f={len(classes)}  "
            f"depth={by_depth}  ({time.time() - t0:.1f}s){note}")
    return fht, res


# ======================================================================
# Check of the axioms on a realization (independent of the construction)
# ======================================================================
def check_fht_axioms(v, E):
    """
    Recompute the components of the realization (top edges = inclusion-maximal edges,
    blocks = connected components of the remaining edges, recursively) and check
    (A1)-(A4) and the component-wise (A5).  Returns the depth, or raises AssertionError.
    """
    full = (1 << v) - 1

    def comps(edges, universe):
        """connected components (vertex sets) of the hypergraph edges on `universe`"""
        left = list(edges)
        out = []
        covered = 0
        while left:
            reach, grp, changed = left[0], [left[0]], True
            rest = left[1:]
            while changed:
                changed, keep = False, []
                for e in rest:
                    if e & reach:
                        reach |= e
                        grp.append(e)
                        changed = True
                    else:
                        keep.append(e)
                rest = keep
            out.append((reach, grp))
            covered |= reach
            left = rest
        for x in bits(universe & ~covered):
            out.append((1 << x, []))
        return out

    def depth(U, edges):
        """depth of the FHT carried by vertex set U with realization `edges` (edges inside U)"""
        if popcount(U) == 1:
            assert not edges
            return 0
        maxe = [e for e in edges if not any(e != f and (e & f) == e for f in edges)]
        assert maxe, "no top edge"
        assert comps(maxe, U)[0][0] == U and len(comps(maxe, U)) == 1, "top level not connected/covering"
        rest = [e for e in edges if e not in maxe]
        blocks = comps(rest, U)
        assert len(blocks) >= 2, "single block"
        ds = []
        for B, grp in blocks:
            ds.append(depth(B, grp))
        for e in maxe:
            assert popcount(e) >= 2
            for B, _ in blocks:           # (A4): an edge meeting a block contains it strictly
                if e & B:
                    assert (e & B) == B and e != B, "A4 violated"
        return 1 + max(ds)

    assert comps(E, full)[0][0] == full and len(comps(E, full)) == 1
    return depth(full, list(E))


# ======================================================================
# Catalogue of the FHT (table of Article 1) and data of the recurrence
# ======================================================================
def set_str(mask):
    return "".join(str(x + 1) for x in bits(mask))


def edges_space_str(E):
    es = sorted(E, key=lambda e: (popcount(e), tuple(bits(e))))
    return " ".join(set_str(e) for e in es)


def build_rows(fht, vtab):
    """One row per FHT, identifiers F_1, F_2, ... in the order (v, depth, canonical key)."""
    rows, gid = [], 0
    for v in range(1, vtab + 1):
        for c in fht[v]:
            gid += 1
            E = c["edges"]
            if v > 1:           # F_1 (one vertex) is the axiom at level 0
                assert check_fht_axioms(v, E) == c["depth"], f"axioms fail for F_{gid}"
            rows.append({"id": gid, "v": v, "N": c["depth"], "num_edges": len(E),
                         "edges_str": edges_space_str(E), "edges": E})
    return rows


def latex_rows(rows):
    """Body of the longtable of Article 1: \\midrule between the groups (v, N)."""
    out, prev = [], None
    for r in rows:
        cur = (r["v"], r["N"])
        if prev is not None and cur != prev:
            out.append("\\midrule")
        prev = cur
        out.append(f"$\\F_{{{r['id']}}}$ & {r['edges_str']} & {r['v']} & {r['N']} & {r['num_edges']} \\\\")
    return "\n".join(out)


def export_table(rows, outdir):
    os.makedirs(outdir, exist_ok=True)
    p = lambda name: os.path.join(outdir, name)
    with open(p("FHT_Table.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["FHT", "Edges", "V", "N", "E"])
        for r in rows:
            w.writerow([f"F_{r['id']}", r["edges_str"], r["v"], r["N"], r["num_edges"]])
    body = latex_rows(rows)
    with open(p("FHT_Table_rows.tex"), "w", encoding="utf-8") as f:
        f.write(body + "\n")
    with open(p("FHT_Table.tex"), "w", encoding="utf-8") as f:
        f.write("\\documentclass{article}\n\\usepackage{amsmath,amssymb,booktabs,longtable}\n"
                "\\newcommand{\\F}{\\mathcal{F}}\n\\begin{document}\n"
                "\\begin{longtable}{c p{5cm} c c c}\n\\toprule\n"
                "FHT & Edges & $\\lvert V\\rvert$ & $N$ & $\\lvert E\\rvert$ \\\\\n\\midrule\n\\endhead\n"
                + body + "\n\\bottomrule\n\\end{longtable}\n\\end{document}\n")
    log(f"table written: {p('FHT_Table.csv')} ({len(rows)} FHT), FHT_Table_rows.tex, FHT_Table.tex")
    try:
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-output-directory", outdir,
                            p("FHT_Table.tex")], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for ext in (".aux", ".log"):
            if os.path.exists(p("FHT_Table" + ext)):
                os.remove(p("FHT_Table" + ext))
        if r.returncode == 0:
            log("   FHT_Table.pdf compiled")
    except FileNotFoundError:
        log("   (pdflatex not found: no pdf)")


def fix_tables(vmax):
    """Fix(lam) for every cycle type of S_m, 2 <= m <= 7 (vmax <= 8) or 8; flag = how the value was obtained."""
    mmax = min(vmax - 1, GENERIC_MAX_M)
    out = {}
    for m in range(2, mmax + 1):
        L = []
        for lam in partitions(m):
            val = fix_generic(m, lam)
            flag = ""
            if lam == (1,) * m and m >= 7:
                flag = "dedekind"
            if m == 8 and ((2 in lam and all(c <= 2 for c in lam)) or (lam[0] == 3 and all(c == 1 for c in lam[1:]))):
                flag = "rec"
            L.append([list(lam), val, flag])
        out[str(m)] = L
    return out


def growth_table(res, vmax):
    g = {}
    for v in range(4, vmax + 1):
        I, f = res[v]["I"], res[v]["f"]
        r = I / f
        if r >= 0.01:
            txt = sci = f"{r:.2f}"
        else:
            e = int(f"{r:e}".split("e")[1])
            sci = txt = f"{r / 10 ** e:.1f}\\cdot10^{{{e}}}"
        H1 = host_orbits((1, v - 2))        # one block of 2 vertices + (v-2) singletons
        g[str(v)] = {"I_over_f": txt, "I_over_f_sci": sci, "share_pct": f"{100 * H1 / I:.1f}", "H1": H1}
    return g


def export_results(res, rows, vmax, outdir):
    counts = {str(v): {"S": res[v]["S"], "I": res[v]["I"], "f": res[v]["f"],
                       "by_depth": {str(d): x for d, x in res[v]["by_depth"].items()}} for v in res}
    D = {"counts": counts, "fix": fix_tables(vmax), "growth": growth_table(res, vmax),
         "rows_latex1": latex_rows(rows)}
    path = os.path.join(outdir, "FHT_Enum_results.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(D, f, indent=1)
    log(f"results written: {path}")


# ======================================================================
# Self-checks (no pynauty needed)
# ======================================================================
def selfcheck(deep=False):
    log("[selfcheck] identity term of the DFS versus Dedekind formula (m = 2..6)")
    ok = True
    for m in range(2, 7):
        t0 = time.time()
        d = _fix_python_dfs(m, (1,) * m)
        ref = labelled_connected_clutters(m)
        good = (d == ref)
        ok &= good
        log(f"   m={m}: DFS={d}  Dedekind={ref}  {'OK' if good else 'MISMATCH'}  ({time.time() - t0:.1f}s)")
    log("[selfcheck] flat counts S(v) = H((v)) versus S(v) derived from A003182 (v = 2..7)")
    for v in range(2, 8):
        s = host_orbits((v,))
        good = s == FLAT_S[v]
        ok &= good
        log(f"   v={v}: {s}  A003182-derived={FLAT_S[v]}  {'OK' if good else 'MISMATCH'}")
    log("   (v = 7 checks the Dedekind number D(7), the C helper and all 15 cycle types of S_7)")
    log("[selfcheck] Dedekind-type recurrence (formerly V0_fht_upset.py) versus direct enumeration, m = 3..7")
    try:
        for m in range(3, 8):
            for lam in partitions(m):
                lam = tuple(lam)
                if lam == (1,) * m or 2 not in lam or any(c > 2 for c in lam):
                    continue
                got = fix_by_recurrence(m, lam, fix_generic, z, partitions, DEDEKIND)
                ref = fix_generic(m, lam)
                good = got == ref
                ok &= good
                log(f"   m={m} {lam}: recurrence={got}  enumeration={ref}  {'OK' if good else 'MISMATCH'}")
    except RuntimeError as e:
        log(f"   skipped ({e})")
    if deep:
        log("[selfcheck] m = 8 : Burnside sum over S_8 of ALL cycle types (the involutive ones by recurrence, "
            "the others by the C helper) versus 8! * S(8), S(8) derived from A003182(8)")
        s = host_orbits((8,))
        good = s == FLAT_S[8]
        ok &= good
        log(f"   v=8: {s}  A003182-derived={FLAT_S[8]}  {'OK' if good else 'MISMATCH'}"
            "  (independent check: no term is deduced from S(8))")
    log("[selfcheck] " + ("all good" if ok else "FAILURE"))
    return ok


# ======================================================================
# Command line and main
# ======================================================================
def parse_args():
    ap = argparse.ArgumentParser(description="FHT: brute-force catalogue (v<=6), recurrence (v<=9), table of Article 1")
    ap.add_argument("--vmax", type=int, default=LIMIT_RECURRENCE, help="recurrence up to v (<= 9, default 9)")
    ap.add_argument("--brute-max", type=int, default=LIMIT_BRUTE,
                    help="brute force up to v (<= 6, default 6; 0 = skip)")
    ap.add_argument("--table-max", type=int, default=5, help="table of Article 1 up to v (default 5)")
    ap.add_argument("--outdir", default=os.path.join(HERE, "catalogues"), help="directory of the output files")
    ap.add_argument("--selfcheck", action="store_true", help="run the self-checks and exit")
    ap.add_argument("--deep", action="store_true",
                    help="with --selfcheck: also check m = 8 (needs the cache or hours of C computation)")
    return ap.parse_args()


def main():
    a = parse_args()
    log("=" * 70 + "\n V0_FHT_Enum -- FRACTAL HYPER-TREES (axiom A5 component-wise)\n" + "=" * 70)
    if a.selfcheck:
        sys.exit(0 if selfcheck(a.deep) else 1)
    if not (1 <= a.vmax <= LIMIT_RECURRENCE):
        sys.exit(f"Refused: --vmax must lie in 1..{LIMIT_RECURRENCE}")
    if not (0 <= a.brute_max <= LIMIT_BRUTE):
        sys.exit(f"Refused: --brute-max must lie in 0..{LIMIT_BRUTE}")
    try:
        rec = run_recurrence(a.vmax)
        fht, brute = (run_brute(a.brute_max) if a.brute_max else (None, {}))
    except GuardError as err:
        sys.exit(f"Stopped by a guard-rail: {err}")

    ok = True
    if brute:
        log("\n[comparison] brute force versus recurrence")
        for v in brute:
            if v in rec:
                same = brute[v]["f"] == rec[v]["f"] and brute[v]["by_depth"] == rec[v]["by_depth"]
                ok &= same
                log(f"   v={v}: f={brute[v]['f']}  {'OK' if same else 'DIFFERENT'}")
    log("\n[species equation] independent cross-check of the totals f(v)")
    sp = run_species(a.vmax, mcap=min(a.vmax - 1, GENERIC_MAX_M) if a.vmax > 1 else 1)
    for v in range(1, a.vmax + 1):
        same = sp[v] == rec[v]["f"]
        ok &= same
        log(f"   v={v}: f={sp[v]}  {'OK' if same else 'DIFFERENT'}")
    log("\n[recurrence] f(v) = " + ", ".join(str(rec[v]["f"]) for v in range(1, a.vmax + 1)))

    os.makedirs(a.outdir, exist_ok=True)
    vtab = min(a.table_max, a.brute_max)
    rows = None
    if vtab >= 1:
        rows = build_rows(fht, vtab)
        log(f"\n[table] {len(rows)} FHT for v <= {vtab}; axioms (A1)-(A5) verified on every realization")
        export_table(rows, a.outdir)
    else:
        log("\n[table] not produced (needs the brute force, --brute-max >= 1)")
    if rows is not None:
        export_results(rec, rows, a.vmax, a.outdir)
    if not ok:
        sys.exit("FAILURE: inconsistency detected")


# ======================================================================
# Dedekind-type recurrence for Fix(lambda), involutive types (formerly V0_fht_upset.py)
# ======================================================================

_UP_DMAX = 8


def _upcompile(src, exe):
    path_exe, path_src = os.path.join(HERE, exe), os.path.join(HERE, src)
    if os.path.exists(path_exe):
        return path_exe
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if cc is None or not os.path.exists(path_src):
        raise RuntimeError(f"cannot build {exe}: no C compiler or missing {src}")
    r = subprocess.run([cc, "-O2", "-o", path_exe, path_src])
    if r.returncode:
        raise RuntimeError(f"compilation of {src} failed")
    return path_exe


def _upupsets(r):
    """all upsets of the Boolean lattice 2^r, as integer masks over the 2^r subsets"""
    U = [0, 1]
    for k in range(r):
        w = 1 << k
        U = [a | (b << w) for a in U for b in U if a & ~b == 0]
    return U


def _upcount_A_3cycle(r):
    """A((3,1^r)): the rotation (a b c) identifies the subsets of {a,b,c} of equal size, so an invariant upset is a
    chain e <= s <= d <= t of upsets of 2^r (values at the sizes 0, 1, 2, 3):  A = sum_{s<=d} Nsub(s) Nsup(d)."""
    import numpy as np
    U = _upupsets(r)
    arr = np.array(U, dtype=object) if r > 5 else np.array(U, dtype=np.uint64)
    n = len(U)
    nsub = np.zeros(n, dtype=np.int64); nsup = np.zeros(n, dtype=np.int64)
    for i in range(n):
        a = arr[i]
        nsub[i] = int(np.count_nonzero((arr & ~a) == 0))
        nsup[i] = int(np.count_nonzero((a & ~arr) == 0))
    total = 0
    for i in range(n):                       # s = U[i]
        a = arr[i]
        mask = (a & ~arr) == 0               # d >= s
        total += int(nsub[i]) * int(nsup[mask].sum())
    return total


def _upcount_A(lam, dedekind):
    """Number of antichains (families of subsets, empty family and {emptyset} included) invariant under a
    permutation of cycle type lam in {1,2}^*; lam must be (1^m) or contain a 2-cycle."""
    lam = tuple(sorted(lam, reverse=True))
    m = sum(lam)
    if lam == (1,) * m:
        return dedekind[m]
    if lam[0] == 3 and all(c == 1 for c in lam[1:]):            # (3,1^r)
        return _upcount_A_3cycle(m - 3)
    if any(c > 2 for c in lam) or 2 not in lam:
        raise ValueError("_upcount_A: needs an involution type or (3,1^r)")
    rest = list(lam)
    rest.remove(2)
    if sum(rest) > 6:
        raise ValueError("_upcount_A: remainder too large")
    out = subprocess.run([_upcompile("V0_fht_upset.c", "V0_fht_upset_bin")] + [str(c) for c in rest],
                         capture_output=True, text=True).stdout
    return int(out.split("count=")[1].split()[0])


# ---------------- truncated polynomials in p_1..p_8 (dict: partition tuple -> Fraction)
def _uppadd(f, g, c=1):
    h = dict(f)
    for k, v in g.items():
        h[k] = h.get(k, 0) + c * v
    return {k: v for k, v in h.items() if v != 0}


def _uppmul(f, g):
    h = {}
    for a, x in f.items():
        sa = sum(a)
        for b, y in g.items():
            if sa + sum(b) > _UP_DMAX:
                continue
            k = tuple(sorted(a + b, reverse=True))
            h[k] = h.get(k, 0) + x * y
    return h


def _uppscale(f, k):
    return {tuple(sorted((i * k for i in a), reverse=True)): v for a, v in f.items() if sum(a) * k <= _UP_DMAX}


def _uppexp(f):
    g = {}
    for k in range(1, _UP_DMAX + 1):
        g = _uppadd(g, _uppscale(f, k), Fraction(1, k))
    res, term = {(): Fraction(1)}, {(): Fraction(1)}
    for n in range(1, _UP_DMAX + 1):
        term = _uppmul(term, g)
        term = {k: v / n for k, v in term.items()}
        res = _uppadd(res, term)
    return res


def fix_by_recurrence(m, lam, fixfun, zfun, partitions, dedekind):
    """Fix(lam) on m points by the recurrence above.  fixfun(mm, mu) gives the (already known) Fix numbers
    for mm < m; zfun(mu) = z_mu; partitions(mm) lists the partitions of mm."""
    lam = tuple(sorted(lam, reverse=True))
    base = {(1,): Fraction(2)}
    for mm in range(2, m):
        for mu in partitions(mm):
            base[tuple(mu)] = Fraction(fixfun(mm, tuple(mu)), zfun(mu))
    lower = _uppexp(base).get(lam, Fraction(0))
    val = (_upcount_A(lam, dedekind) - 1) - zfun(lam) * lower
    assert val.denominator == 1, "recurrence: non-integer result"
    return int(val)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\nInterrupted.")
