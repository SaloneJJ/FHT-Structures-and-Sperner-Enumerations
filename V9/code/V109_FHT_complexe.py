"""V109 - Bibliotheque : complexe hierarchique (stella, blocs 2^3, tore L0^3) en entier.
Liens (h horizontaux, v verticaux), plaquettes (B magnetiques, T triangles, S carres electriques),
cellules de dimension 3 (cubes horizontaux, prismes verticaux), cellules de dimension 4 (hyperprismes),
et trois definitions du chiffre de temps q = F u F sur l'hypercube (l, x) :
  'v103'  : cup cubique utilise de V101 a V106 (sans test de degenerescence) ;
  'serre' : cup cubique de Serre tire en arriere sur l'hypercube (faces degenerees = 0) ;
  'aw'    : Alexander-Whitney sur la triangulation de Freudenthal de l'hypercube (diagonales par chemin canonique).
JJ SALONE, FHT 2026"""
import itertools, math
import numpy as np, scipy.sparse as sp

PAIRS = [(0, 1), (0, 2), (1, 2)]


class Complexe:
    def __init__(self, L0):
        self.L0 = L0
        self.NLEV = int(math.log2(L0)) - 1
        self.Lv = [L0 >> l for l in range(self.NLEV + 1)]
        off = {}; n = 0
        for l in range(self.NLEV + 1):
            off[('h', l)] = n; n += 3 * self.Lv[l] ** 3
            if l < self.NLEV: off[('v', l)] = n; n += self.Lv[l] ** 3
        self.off = off; self.nl = n
        self._plaquettes()
        self._cellules3()
        self._cellules4()

    # --- indices ---
    def sidx(self, l, x):
        L = self.Lv[l]; return ((x[0] % L) * L + (x[1] % L)) * L + (x[2] % L)
    def h(self, l, x, k): return self.off[('h', l)] + 3 * self.sidx(l, x) + k
    def v(self, l, x): return self.off[('v', l)] + self.sidx(l, x)
    @staticmethod
    def sh(x, k, s=1): y = list(x); y[k] += s; return tuple(y)
    def par(self, l, x): return tuple((c % self.Lv[l]) // 2 for c in x)
    def sites(self, l): return itertools.product(range(self.Lv[l]), repeat=3)
    def link_type(self, a):
        """('h', l) ou ('v', l) du lien a"""
        best = None
        for key, o in self.off.items():
            if o <= a and (best is None or o > best[1]): best = (key, o)
        return best[0]

    def _plaquettes(self):
        rows, cols, typ, lev = [], [], [], []
        self.pB, self.pE = {}, {}
        for l in range(self.NLEV + 1):
            for x in self.sites(l):
                for i, j in PAIRS:
                    p = len(typ); self.pB[(l, self.sidx(l, x), i, j)] = p; typ.append('B'); lev.append(l)
                    for a in (self.h(l, x, i), self.h(l, self.sh(x, i), j), self.h(l, self.sh(x, j), i), self.h(l, x, j)):
                        rows.append(p); cols.append(a)
                if l < self.NLEV:
                    for k in range(3):
                        y = self.sh(x, k); P, Pp = self.par(l, x), self.par(l, y)
                        ls = [self.h(l, x, k), self.v(l, y), self.v(l, x)]
                        t = 'T' if P == Pp else 'S'
                        if t == 'S': ls.append(self.h(l + 1, P, k))
                        p = len(typ); self.pE[(l, self.sidx(l, x), k)] = p; typ.append(t); lev.append(l)
                        for a in ls: rows.append(p); cols.append(a)
        self.np_ = len(typ)
        A = sp.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=(self.np_, self.nl))
        A.data %= 2; A.eliminate_zeros()
        self.A = A; self.typ = np.array(typ); self.lev = np.array(lev)

    def B(self, l, x, i, j): return self.pB[(l, self.sidx(l, x), i, j)]
    def E(self, l, x, k): return self.pE[(l, self.sidx(l, x), k)]

    def _cellules3(self):
        """cubes horizontaux (type 'C') et prismes verticaux (type 'P', 'P6' si face haute non degeneree)"""
        rows, cols, typ, lev = [], [], [], []
        self.c3 = {}
        for l in range(self.NLEV + 1):
            for x in self.sites(l):
                c = len(typ); self.c3[('C', l, self.sidx(l, x))] = c; typ.append('C'); lev.append(l)
                for (i, j) in PAIRS:
                    k = 3 - i - j
                    for p in (self.B(l, x, i, j), self.B(l, self.sh(x, k), i, j)): rows.append(c); cols.append(p)
                if l < self.NLEV:
                    for (i, j) in PAIRS:
                        faces = [self.B(l, x, i, j), self.E(l, x, i), self.E(l, self.sh(x, j), i),
                                 self.E(l, x, j), self.E(l, self.sh(x, i), j)]
                        top = (x[i] % 2 == 1) and (x[j] % 2 == 1)
                        if top: faces.append(self.B(l + 1, self.par(l, x), i, j))
                        c = len(typ); self.c3[('P', l, self.sidx(l, x), i, j)] = c
                        typ.append('P6' if top else 'P5'); lev.append(l)
                        for p in faces: rows.append(c); cols.append(p)
        self.n3 = len(typ)
        G = sp.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=(self.n3, self.np_))
        G.data %= 2; G.eliminate_zeros()
        self.G = G; self.typ3 = np.array(typ); self.lev3 = np.array(lev)

    def _cellules4(self):
        rows, cols = [], []; n = 0
        for l in range(self.NLEV):
            for x in self.sites(l):
                cells = [self.c3[('C', l, self.sidx(l, x))]]
                for (i, j) in PAIRS:
                    k = 3 - i - j
                    cells.append(self.c3[('P', l, self.sidx(l, x), i, j)])
                    cells.append(self.c3[('P', l, self.sidx(l, self.sh(x, k)), i, j)])
                if all(c % 2 == 1 for c in x): cells.append(self.c3[('C', l + 1, self.sidx(l + 1, self.par(l, x)))])
                for c in cells: rows.append(n); cols.append(c)
                n += 1
        self.n4 = n
        H4 = sp.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=(n, self.n3))
        H4.data %= 2; H4.eliminate_zeros(); self.H4 = H4

    # --- hypercubes : liens tires en arriere ---
    def hypercube_edges(self, l, x):
        """dict (u, w) -> indice de lien ou -1 (degenere), pour les 32 aretes de I^4 ;
        sommets u = (e0, e1, e2, t) ; direction 3 = verticale."""
        E = {}
        for u in itertools.product((0, 1), repeat=4):
            for d in range(4):
                if u[d] == 1: continue
                w = list(u); w[d] = 1; w = tuple(w)
                xs = tuple(x[m] + u[m] for m in range(3))
                if d == 3:
                    E[(u, w)] = self.v(l, xs)
                elif u[3] == 0:
                    E[(u, w)] = self.h(l, xs, d)
                else:
                    P = self.par(l, xs); Q = self.par(l, self.sh(xs, d))
                    E[(u, w)] = -1 if P == Q else self.h(l + 1, P, d)
        return E


def _vert_add(u, d):
    w = list(u); w[d] = 1; return tuple(w)


def faces2_I4():
    """faces de dimension 2 de I^4 : (directions (a,b), point de base u)"""
    out = []
    for a, b in itertools.combinations(range(4), 2):
        others = [m for m in range(4) if m not in (a, b)]
        for vals in itertools.product((0, 1), repeat=2):
            u = [0] * 4
            for m, vv in zip(others, vals): u[m] = vv
            out.append(((a, b), tuple(u)))
    return out


def face_edges(dirs, u):
    a, b = dirs
    ua = _vert_add(u, a); ub = _vert_add(u, b)
    return [(u, ua), (ua, _vert_add(ua, b)), (u, ub), (ub, _vert_add(ub, a))]


def serre_terms():
    """cup de Serre de deux 2-cochaines sur I^4 : liste de paires (face avant A, face arriere B)"""
    terms = []
    for A_ in itertools.combinations(range(4), 2):
        Bd = tuple(m for m in range(4) if m not in A_)
        front = (A_, (0, 0, 0, 0))
        ub = [0] * 4
        for m in A_: ub[m] = 1
        back = (Bd, tuple(ub))
        terms.append((front, back))
    return terms


def freudenthal():
    """24 simplexes de I^4 : suites de sommets v0..v4 (chemins monotones)"""
    out = []
    for perm in itertools.permutations(range(4)):
        u = (0, 0, 0, 0); seq = [u]
        for d in perm: u = _vert_add(u, d); seq.append(u)
        out.append(seq)
    return out


def canonical_path(u, w):
    """aretes du chemin canonique de u a w (u <= w), coordonnees croissantes dans l'ordre 0,1,2,3"""
    path = []; cur = u
    for d in range(4):
        if w[d] == 1 and cur[d] == 0:
            nxt = _vert_add(cur, d); path.append((cur, nxt)); cur = nxt
    return path


class Cups:
    """precalcule, pour chaque hypercube, les listes d'indices de liens permettant d'evaluer
    les trois versions de q a partir du vecteur b (bits de liens) ; et la version V103 a partir de F."""
    def __init__(self, K):
        self.K = K
        f2 = faces2_I4(); fidx = {f: n for n, f in enumerate(f2)}
        sterms = [(fidx[(fr[0], fr[1])], fidx[(bk[0], bk[1])]) for fr, bk in serre_terms()]
        simp = freudenthal()
        self.per_level = []
        for l in range(K.NLEV):
            faceL = []   # pour chaque hypercube : 24 faces x 4 liens (indice ou -1)
            awL = []     # pour chaque hypercube : 24 simplexes x (tri1 liens, tri2 liens) en chemins
            for x in K.sites(l):
                E = K.hypercube_edges(l, x)
                faceL.append([[E[e] for e in face_edges(dirs, u)] for dirs, u in f2])
                sims = []
                for s in simp:
                    def tri(a, b, c):
                        ed = canonical_path(a, b) + canonical_path(b, c) + canonical_path(a, c)
                        return [E[e] for e in ed]
                    t1 = tri(s[0], s[1], s[2]); t2 = tri(s[2], s[3], s[4])
                    sims.append((t1, t2))
                awL.append(sims)
            self.per_level.append(dict(face=np.array(faceL), aw=awL))
        self.sterms = np.array(sterms)
        # AW : padding a longueur fixe (chemins de longueur <= 2+2+... : tri <= 2+2+4 = 8 liens)
        for D in self.per_level:
            n = len(D['aw']); M1 = -np.ones((n, 24, 8), dtype=np.int64); M2 = -np.ones((n, 24, 8), dtype=np.int64)
            for h_, sims in enumerate(D['aw']):
                for s, (t1, t2) in enumerate(sims):
                    M1[h_, s, :len(t1)] = t1; M2[h_, s, :len(t2)] = t2
            D['aw1'] = M1; D['aw2'] = M2; del D['aw']

    @staticmethod
    def _gather(b, idx):
        bb = np.concatenate([b, [0]]).astype(np.int64)
        return bb[np.where(idx < 0, len(b), idx)]

    def q_serre(self, b, l):
        Fh = self._gather(b, self.per_level[l]['face']).sum(-1) % 2   # (n, 24)
        return (Fh[:, self.sterms[:, 0]] * Fh[:, self.sterms[:, 1]]).sum(-1) % 2

    def q_aw(self, b, l):
        D = self.per_level[l]
        t1 = self._gather(b, D['aw1']).sum(-1) % 2; t2 = self._gather(b, D['aw2']).sum(-1) % 2
        return (t1 * t2).sum(-1) % 2

    def q_serre_z4(self, b, l):
        """relevement entier naif (somme des produits dans Z, reduite mod 4) ; sert au chapitre theta"""
        Fh = self._gather(b, self.per_level[l]['face']).sum(-1) % 2
        return (Fh[:, self.sterms[:, 0]] * Fh[:, self.sterms[:, 1]]).sum(-1) % 4


def qpairs_v103(K):
    """paires de plaquettes du cup V103 (pour comparaison)"""
    out = []
    for l in range(K.NLEV):
        P1, P2 = [], []
        for x in K.sites(l):
            a1, a2 = [], []
            for k in range(3):
                i, j = [t for t in range(3) if t != k]
                Pp = K.par(l, K.sh(x, k))
                a1.append(K.E(l, x, k)); a2.append(K.B(l + 1, Pp, i, j))
                a1.append(K.B(l, x, i, j)); a2.append(K.E(l, K.sh(K.sh(x, i), j), k))
            P1.append(a1); P2.append(a2)
        out.append((np.array(P1), np.array(P2)))
    return out
