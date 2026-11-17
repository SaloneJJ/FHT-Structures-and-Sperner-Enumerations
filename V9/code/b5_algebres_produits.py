"""b5_algebres_produits.py -- V58 etape 1 (JJ SALONE, FHT 2026)
Decomposition reelle des algebres tordues R^eta[G], G = produit de groupes symetriques,
pour toutes les classes eta de H^2(G; Z/2) : classes de facteur (0, a = w1^2, e = Schur, a+e)
et classes mixtes w1(i) w1(j). Methode : algebre des classes tordues, caracteres projectifs,
indicateur de Frobenius-Schur tordu. Restriction aux sous-groupes S4 (facteur, diagonal, graphe S4->S3).
Usage : python3 -I b5_algebres_produits.py S5xS3 | S4xS4xS3 | S4xS4xS5
"""
import sys, itertools, numpy as np
from collections import deque

def sym(n):
    E = list(itertools.permutations(range(n)))
    ix = {g: i for i, g in enumerate(E)}
    N = len(E)
    M = np.array([[ix[tuple(g[h[k]] for k in range(n))] for h in E] for g in E])
    inv = np.array([ix[tuple(np.argsort(g))] for g in E])
    par = np.array([sum(1 for i in range(n) for j in range(i+1, n) if g[i] > g[j]) % 2 for g in E])
    gens = []
    for i in range(n-1):
        t = list(range(n)); t[i], t[i+1] = t[i+1], t[i]; gens.append(ix[tuple(t)])
    # cocycle de Schur par Clifford : t_i -> (e_i - e_{i+1})/sqrt2
    X = np.array([[0,1],[1,0]],complex); Y = np.array([[0,-1j],[1j,0]]); Z = np.diag([1.,-1]).astype(complex); I = np.eye(2)
    m = max(1, n//2)
    def kron(L):
        out = np.eye(1)
        for A in L: out = np.kron(out, A)
        return out
    gam = []
    for j in range(m):
        gam.append(kron([Z]*j+[X]+[I]*(m-j-1))); gam.append(kron([Z]*j+[Y]+[I]*(m-j-1)))
    gam.append(kron([Z]*m))
    T = [(gam[i]-gam[i+1])/np.sqrt(2) for i in range(n-1)]
    Mat = [None]*N; Mat[ix[tuple(range(n))]] = np.eye(2**m, dtype=complex)
    dq = deque([ix[tuple(range(n))]])
    while dq:
        g = dq.popleft()
        for i, s in enumerate(gens):
            h = M[g, s]
            if Mat[h] is None: Mat[h] = Mat[g] @ T[i]; dq.append(h)
    A = np.array(Mat); d = 2**m
    P = np.einsum('gab,hbc->ghac', A, A)
    R = A[M]  # R[g,h] = Mat(gh)
    c = np.real(np.einsum('ghba,ghbc->gh', R.conj(), P)) / d
    assert np.allclose(np.abs(c), 1), 'Clifford non fidele'
    ce = np.rint(c).astype(int)
    ca = np.where(np.outer(par, par) == 1, -1, 1)
    coc = {'0': np.ones((N, N), int), 'a': ca}
    if n >= 4: coc['e'] = ce; coc['a+e'] = ca*ce
    return dict(n=n, E=E, ix=ix, N=N, M=M, inv=inv, par=par, gens=gens, coc=coc)

class Prod:
    def __init__(self, ns):
        self.F = [sym(n) for n in ns]; self.ns = ns
        self.sz = [f['N'] for f in self.F]
        self.st = [int(np.prod(self.sz[i+1:])) for i in range(len(ns))]
        self.N = int(np.prod(self.sz))
    def dec(self, x): return [(x // s) % n for s, n in zip(self.st, self.sz)]
    def enc(self, ds): return sum(d*s for d, s in zip(ds, self.st))
    def mul(self, x, y):
        a, b = self.dec(x), self.dec(y)
        return self.enc([f['M'][u, v] for f, u, v in zip(self.F, a, b)])
    def inv(self, x): return self.enc([f['inv'][u] for f, u in zip(self.F, self.dec(x))])
    def cocfun(self, cls, mixed):
        def c(x, y):
            a, b = self.dec(x), self.dec(y)
            out = np.ones(np.shape(x), int)
            for f, k, u, v in zip(self.F, cls, a, b): out = out * f['coc'][k][u, v]
            for (i, j) in mixed:
                out = out * np.where(self.F[i]['par'][a[i]] * self.F[j]['par'][b[j]] == 1, -1, 1)
            return out
        return c

def algebra(P, H, gens, coc):
    """H : tableau d'indices de G (sous-groupe), gens : generateurs (indices G). Renvoie liste (type, n)."""
    H = np.asarray(H); nH = len(H)
    pos = -np.ones(P.N, int); pos[H] = np.arange(nH)
    e = P.enc([f['ix'][tuple(range(f['n']))] for f in P.F])
    # conjugaisons par generateurs
    conj, lam = [], []
    for h in gens:
        hh = np.full(nH, h); hi = np.full(nH, P.inv(h))
        hx = P.mul(hh, H); y = P.mul(hx, hi)
        conj.append(pos[y])
        lam.append(coc(hh, H) * coc(hx, hi) * coc(hh, hi))
    cid = -np.ones(nH, int); s = np.zeros(nH, int); reps = []; bad = []
    for x0 in range(nH):
        if cid[x0] >= 0: continue
        k = len(reps); reps.append(x0); cid[x0] = k; s[x0] = 1; ok = True
        dq = deque([x0])
        while dq:
            x = dq.popleft()
            for cm, lm in zip(conj, lam):
                y = cm[x]; sy = s[x] * lm[x]
                if cid[y] < 0: cid[y] = k; s[y] = sy; dq.append(y)
                elif s[y] != sy: ok = False
        bad.append(not ok)
    bad = np.array(bad)
    regs = [k for k in range(len(reps)) if not bad[k]]
    rmap = -np.ones(len(reps), int); rmap[regs] = np.arange(len(regs))
    rc = rmap[cid]  # classe reguliere ou -1
    nR = len(regs); size = np.bincount(rc[rc >= 0], minlength=nR)
    a = np.zeros((nR, nR, nR))
    Hinv = P.inv(H)
    for kk, k in enumerate(regs):
        r = H[reps[k]]
        Y = P.mul(Hinv, np.full(nH, r)); py = pos[Y]
        m = (rc >= 0) & (rc[py] >= 0)
        w = s[m] * s[py[m]] * coc(H[m], Y[m])
        np.add.at(a[:, :, kk], (rc[m], rc[py[m]]), w)
    rng = np.random.default_rng(1)
    B = np.tensordot(rng.normal(size=nR), a, 1)
    ev, V = np.linalg.eig(B)
    i0 = rc[pos[e]]
    out = []; tot = 0; irr = []
    X2 = P.mul(H, H); p2 = pos[X2]; c2 = coc(H, H)
    for j in range(nR):
        w = V[:, j] / V[i0, j]
        d = np.sqrt(nH / np.sum(np.abs(w)**2 / size)); d = int(round(d.real))
        chi = w * d / size
        m = rc[p2] >= 0
        nu = np.sum(c2[m] * s[p2[m]] * chi[rc[p2[m]]]) / nH
        nu = int(round(nu.real)); out.append((nu, d))
        chiv = np.where(rc >= 0, s * chi[np.maximum(rc, 0)], 0)
        irr.append((nu, d, chiv))
    res = []
    for nu, d in out:
        if nu == 1: res.append(('R', d))
        elif nu == -1: res.append(('H', d//2))
        else: res.append(('C', d))
    # paires complexes comptees deux fois
    from collections import Counter
    C = Counter(res); fin = []
    for (t, d), k in sorted(C.items(), key=lambda z: (z[0][1], z[0][0])):
        kk = k//2 if t == 'C' else k
        fin += [(t, d)]*kk
    dim = sum({'R': 1, 'C': 2, 'H': 4}[t]*({'H': 2*d}.get(t, d) if False else d)**2 for t, d in fin)
    dim = sum(d*d*{'R': 1, 'C': 2, 'H': 4}[t] for t, d in fin)
    assert dim == nH, (dim, nH, fin)
    algebra.irr = irr
    return fin

def name(fin):
    from collections import Counter
    C = Counter(fin); parts = []
    for (t, d), k in sorted(C.items(), key=lambda z: (z[0][1], 'RCH'.index(z[0][0]))):
        f = {'R': 'R', 'C': 'C', 'H': 'H'}[t]
        s = f if d == 1 else 'M%d(%s)' % (d, f)
        parts.append(s if k == 1 else '%s^%d' % (s, k))
    return ' + '.join(parts)

def flags(fin):
    has = lambda t, d: (t, d) in fin
    AF = has('C', 1) and has('H', 1) and has('C', 3)
    PS = fin.count(('H', 1)) >= 2 and has('C', 4)
    return AF, PS

def pattern(irrG, Hs, irrH):
    nH = len(Hs)
    def mult(cg, ch): return np.real(np.sum(cg[Hs] * np.conj(ch)) / nH)
    one = [i for i, (nu, d, c) in enumerate(irrH) if d == 1]
    two = [i for i, (nu, d, c) in enumerate(irrH) if d == 2 and nu == -1]
    thr = [i for i, (nu, d, c) in enumerate(irrH) if d == 3]
    HL = HR = M4 = 0
    for nu, d, cg in irrG:
        m = [round(mult(cg, c[2]), 6) for c in irrH]
        if sum(x*irrH[i][1] for i, x in enumerate(m)) != d: continue  # restriction non portee par A_F
        if nu == -1 and d == 2:
            if sum(m[i] for i in two) == 1: HL += 1
            elif sum(m[i] for i in one) == 2: HR += 1
        if nu == 0 and d == 4 and sum(m[i] for i in one) == 1 and sum(m[i] for i in thr) == 1: M4 += 1
    t = []
    if HL and HR and M4: t.append('CCvS')
    if HL >= 2 and M4: t.append('LR')
    return '+'.join(t) or 'af'

def phi_S4_S3(F4, F3):
    pairs = [((0,1),(2,3)), ((0,2),(1,3)), ((0,3),(1,2))]
    key = lambda p: frozenset(frozenset(q) for q in p)
    K = [key(p) for p in pairs]
    out = []
    for g in F4['E']:
        img = [K.index(frozenset(frozenset(g[i] for i in q) for q in pk)) for pk in K]
        out.append(F3['ix'][tuple(img)])
    return out

def subgroups(P):
    F4 = sym(4); res = []
    modes = []
    for f in P.F:
        modes.append(['1', 'id', 'sgn', 'phi'] if f['n'] >= 4 else ['1', 'phi', 'sgn'])
    for combo in itertools.product(*modes):
        if 'id' not in combo: continue
        imgs = []
        for f, md in zip(P.F, combo):
            n = f['n']
            if md == '1': imgs.append([f['ix'][tuple(range(n))]]*24)
            elif md == 'id': imgs.append([f['ix'][tuple(list(g)+list(range(4, n)))] for g in F4['E']])
            elif md == 'sgn':
                t = list(range(n)); t[0], t[1] = 1, 0
                imgs.append([f['ix'][tuple(t)] if F4['par'][k] else f['ix'][tuple(range(n))] for k in range(24)])
            else:
                F3 = sym(3); ph = phi_S4_S3(F4, F3)
                imgs.append([f['ix'][tuple(list(F3['E'][ph[k]]) + list(range(3, n)))] for k in range(24)])
        H = np.array([P.enc([im[k] for im in imgs]) for k in range(24)])
        gens = [H[k] for k in F4['gens']]
        res.append(('/'.join(combo), H, gens))
    return res

def run(ns, lo=0, hi=10**9):
    P = Prod(ns); L = len(ns)
    pairs = list(itertools.combinations(range(L), 2))
    opts = [list(f['coc'].keys()) for f in P.F]
    H = np.arange(P.N)
    G = [P.enc([f['gens'][i] if j == k else f['ix'][tuple(range(f['n']))] for k, f in enumerate(P.F)])
         for j, f in enumerate(P.F) for i in range(len(f['gens']))]
    subs = subgroups(P)
    if lo == 0: print('G =', ' x '.join('S%d' % n for n in ns), '|G| =', P.N, ' classes :', int(np.prod([len(o) for o in opts]))*2**len(pairs)); sys.stdout.flush()
    idx = -1
    for cls in itertools.product(*opts):
        for zs in itertools.product([0, 1], repeat=len(pairs)):
            idx += 1
            if idx < lo or idx >= hi: continue
            mixed = [p for p, z in zip(pairs, zs) if z]
            coc = P.cocfun(cls, mixed)
            fin = algebra(P, H, G, coc); AF, PS = flags(fin); irrG = algebra.irr
            sub = []
            for lab, Hs, gs in subs:
                fs = algebra(P, Hs, gs, coc)
                if flags(fs)[0]: sub.append(lab + ':' + pattern(irrG, Hs, algebra.irr))
            lab = ','.join(cls) + ' | mixte ' + (','.join('%d%d' % p for p in mixed) or '-')
            print('%-30s PS=%d AF=%d  S4->A_F: %-40s %s' % (lab, PS, AF, ' '.join(sub) or '-', name(fin))); sys.stdout.flush()

if __name__ == '__main__':
    a = sys.argv
    run([int(c) for c in a[1].replace('S', '').split('x')], *(int(x) for x in a[2:4]))
