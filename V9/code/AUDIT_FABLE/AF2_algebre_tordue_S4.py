#!/usr/bin/env python3
"""AUDIT FABLE -- A1, registre F6 (Theoreme « A_F as a twisted group algebra ») : verification independante, par les constantes de structure.
On construit l'algebre reelle R^eta[G] de la Definition « Twisted group algebra » (u_g u_h = (-1)^{e(g)e(h)} u_{gh}) dans sa representation
reguliere (matrices reelles |G| x |G|), on calcule son centre, ses idempotents centraux primitifs et, pour chaque composante simple,
(dimension, dimension du centre, signature de la forme trace) -- ce qui la determine :
    centre de dim 2                       -> M_n(C),  dim = 2 n^2
    centre de dim 1, dim = 4 m^2, forme trace de signature (m(2m-1)... ) : on distingue M_n(R) de M_m(H) par la signature de
    la forme (x, y) -> tr(L_x L_y) : pour M_n(R), n(n+1)/2 carres positifs ; pour M_m(H) (dim 4m^2), m(2m-1)... on utilise simplement
    le nombre de directions positives p : M_n(R) : p = n(n+1)/2 ; M_m(H) : p = m(2m-1)  [dim 4 : M_2(R) -> p = 3 ; H -> p = 1].
Aucun caractere, aucun indicateur de Frobenius-Schur n'est utilise : c'est un controle croise de la preuve de l'article.
On traite : (S4, sgn) ; (S4, 1) ; (S3, sgn) ; et tous les couples (G, eta) pour les groupes d'ordre 24 construits ci-dessous
(S4, SL(2,3), A4 x C2, D12, C3 x D4, C4 x S3, C2 x C2 x S3 ... comme groupes de permutations) -- controle partiel de l'unicite.
Usage : python3 AF2_algebre_tordue_S4.py       (numpy)
"""
import itertools, numpy as np

def closure(gens):
    n = len(gens[0]); e = tuple(range(n)); G = {e}; front = [e]
    while front:
        new = []
        for g in front:
            for s in gens:
                h = tuple(g[s[i]] for i in range(n))          # (g o s)(i) = g[s[i]]
                if h not in G: G.add(h); new.append(h)
        front = new
    return sorted(G)

def mul(g, h): return tuple(g[h[i]] for i in range(len(g)))

def characters_pm1(G):
    """tous les homomorphismes G -> {+1,-1} : determines par un sous-groupe d'indice <= 2 ; recherche par les valeurs sur des generateurs"""
    idx = {g: i for i, g in enumerate(G)}; n = len(G)
    # G/G^2[G,G] : on calcule le sous-groupe N engendre par les carres et les commutateurs ; les caracteres sont ceux de G/N (2-groupe elementaire)
    N = set()
    gens = [mul(g, g) for g in G]
    N = set(closure(gens)) if gens else {G[0]}
    # N contient les carres ; G/N est abelien elementaire (les commutateurs sont des produits de carres)
    cosets = []; rep = {}
    for g in G:
        for c in cosets:
            if mul(g, inverse(c[0])) in N: rep[g] = c[0]; c.append(g); break
        else:
            cosets.append([g]); rep[g] = g
    reps = [c[0] for c in cosets]; r = len(reps).bit_length() - 1
    # base de G/N sur F2
    basis = []; span = {reps[0]: ()} if True else {}
    span = {rep[G[0]] if False else rep[tuple(range(len(G[0])))]: 0}
    vec = {rep[tuple(range(len(G[0])))]: 0}
    for g in reps:
        if g in vec: continue
        b = len(basis); basis.append(g)
        for h, v in list(vec.items()):
            vec[rep[mul(h, g)]] = v | (1 << b)
    chars = []
    for mask in range(1 << len(basis)):
        chars.append({g: (-1) ** bin(vec[rep[g]] & mask).count('1') for g in G})
    return chars

def inverse(g):
    inv = [0] * len(g)
    for i, x in enumerate(g): inv[x] = i
    return tuple(inv)

def twisted_regular(G, eta):
    """matrices L_g de la multiplication a gauche par u_g dans la base (u_h)"""
    idx = {g: i for i, g in enumerate(G)}; n = len(G)
    e = {g: (0 if eta[g] == 1 else 1) for g in G}
    L = np.zeros((n, n, n))
    for g in G:
        for h in G:
            L[idx[g], idx[mul(g, h)], idx[h]] = (-1) ** (e[g] * e[h])
    return L

def wedderburn(G, eta, tol=1e-8):
    n = len(G); L = twisted_regular(G, eta)
    # associativite (controle) sur un echantillon
    rng = np.random.default_rng(0)
    for _ in range(20):
        a, b, c = rng.integers(0, n, 3)
        assert np.allclose(L[a] @ L[b] @ np.eye(n)[:, c], (L[a] @ L[b])[:, c])
    # multiplication a droite R_g : R[g][:, h] = u_h u_g
    idx = {g: i for i, g in enumerate(G)}; e = {g: (0 if eta[g] == 1 else 1) for g in G}
    R = np.zeros((n, n, n))
    for g in G:
        for h in G:
            R[idx[g], idx[mul(h, g)], idx[h]] = (-1) ** (e[h] * e[g])
    # centre : x tel que L_g x = R_g x pour tout g
    M = np.vstack([L[i] - R[i] for i in range(n)])
    _, s, Vt = np.linalg.svd(M); Z = Vt[np.sum(s > tol):].T                  # base du centre, colonnes
    zc = Z.shape[1]
    # element central generique -> decomposition par les valeurs propres de L_z (commutant de l'algebre = multiplications a droite)
    z = Z @ rng.normal(size=zc); Lz = sum(z[i] * L[i] for i in range(n))
    w, U = np.linalg.eig(Lz)
    # regrouper les valeurs propres : une composante de centre R donne une valeur reelle, de centre C une paire conjuguee
    comps = []
    used = np.zeros(n, bool)
    for i in range(n):
        if used[i]: continue
        same = np.abs(w - w[i]) < 1e-6; conj = np.abs(w - np.conj(w[i])) < 1e-6
        grp = same | conj; used |= grp
        comps.append((int(grp.sum()), 2 if abs(w[i].imag) > 1e-6 else 1, grp))
    out = []
    for dim, cdim, grp in comps:
        # projecteur spectral reel sur la composante
        P = np.real(U[:, grp] @ np.linalg.pinv(U)[grp, :])
        B = np.linalg.svd(P)[0][:, :dim]                                     # base orthonormee de la composante (ideal bilatere)
        # forme trace restreinte : T(x, y) = tr(L_x L_y)
        Lb = [sum(B[i, a] * L[i] for i in range(n)) for a in range(dim)]
        T = np.array([[np.trace(Lb[a] @ Lb[b]) for b in range(dim)] for a in range(dim)])
        ev = np.linalg.eigvalsh((T + T.T) / 2); p = int((ev > tol).sum()); q = int((ev < -tol).sum())
        if cdim == 2:
            m = int(round((dim / 2) ** .5)); name = f'M_{m}(C)' if m > 1 else 'C'
        else:
            m = int(round(dim ** .5))
            if p == m * (m + 1) // 2: name = f'M_{m}(R)' if m > 1 else 'R'
            elif m % 2 == 0 and p == (m // 2) * (m - 1): name = f'M_{m//2}(H)' if m > 2 else 'H'
            else: name = f'?(dim {dim}, p {p}, q {q})'
        out.append((name, dim, cdim, (p, q)))
    return sorted(out, key=lambda t: (t[1], t[0])), zc

def perm(n, *cycles):
    p = list(range(n))
    for c in cycles:
        for i in range(len(c)): p[c[i]] = c[(i + 1) % len(c)]
    return tuple(p)

def sign(g):
    s = 1; seen = [False] * len(g)
    for i in range(len(g)):
        if not seen[i]:
            j = i; l = 0
            while not seen[j]: seen[j] = True; j = g[j]; l += 1
            if l % 2 == 0: s = -s
    return s

def nclasses(G):
    seen = set(); k = 0
    for g in G:
        if g in seen: continue
        k += 1
        for h in G: seen.add(mul(mul(h, g), inverse(h)))
    return k

if __name__ == '__main__':
    S4 = closure([perm(4, (0, 1)), perm(4, (0, 1, 2, 3))]); assert len(S4) == 24
    S3 = closure([perm(3, (0, 1)), perm(3, (0, 1, 2))])
    print('R^sgn[S4] :', wedderburn(S4, {g: sign(g) for g in S4}))
    print('R[S4]     :', wedderburn(S4, {g: 1 for g in S4}))
    print('R^sgn[S3] :', wedderburn(S3, {g: sign(g) for g in S3}))
    res = wedderburn(S4, {g: sign(g) for g in S4})[0]
    ok = sorted(n for n, *_ in res) == ['C', 'H', 'M_3(C)']
    print('  => R^sgn[S4] = C + H + M_3(C) :', 'CONFORME' if ok else 'ECART')
    # ---- unicite : les 15 groupes d'ordre 24 comme groupes de permutations ; nombre de classes et algebres tordues
    groups = {
        'S4': [perm(4, (0, 1)), perm(4, (0, 1, 2, 3))],
        'A4 x C2': [perm(6, (0, 1, 2)), perm(6, (0, 1), (2, 3)), perm(6, (4, 5))],
        'SL(2,3)': [perm(8, (0, 1, 2, 3), (4, 5, 6, 7)), perm(8, (0, 4, 2, 6), (1, 7, 3, 5)), perm(8, (1, 4, 5), (3, 6, 7))],
        'D12': [perm(12, tuple(range(12))), perm(12, (1, 11), (2, 10), (3, 9), (4, 8), (5, 7))],
        'C24': [perm(24, tuple(range(24)))],
        'C12 x C2': [perm(14, tuple(range(12))), perm(14, (12, 13))],
        'C6 x C2 x C2': [perm(10, tuple(range(6))), perm(10, (6, 7)), perm(10, (8, 9))],
        'C3 x D4': [perm(7, (0, 1, 2)), perm(7, (3, 4, 5, 6)), perm(7, (3, 5))],
        'C3 x Q8': [perm(11, (0, 1, 2)), perm(11, (3, 4, 5, 6), (7, 8, 9, 10)), perm(11, (3, 7, 5, 9), (4, 10, 6, 8))],
        'C4 x S3': [perm(7, (0, 1, 2, 3)), perm(7, (4, 5)), perm(7, (4, 5, 6))],
        'C2 x C2 x S3': [perm(7, (0, 1)), perm(7, (2, 3)), perm(7, (4, 5)), perm(7, (4, 5, 6))],
    }
    # Dic3 = <a, x | a^6 = 1, x^2 = a^3, x a x^-1 = a^-1> : representation reguliere ; Dic6, C3:C8, C3:D4 de meme par presentation
    def from_presentation(order, elems_mul):
        """elems_mul : fonction (i, j) -> k sur range(order) ; renvoie les permutations regulieres gauches"""
        return [tuple(elems_mul(i, j) for j in range(order)) for i in range(order)]
    def dicyclic(n):          # ordre 4n : elements a^k x^e, k mod 2n, e in {0,1}
        def m(i, j):
            k1, e1 = divmod(i, 2); k2, e2 = divmod(j, 2)
            k = (k1 + (k2 if e1 == 0 else -k2)) % (2 * n); e = e1 + e2
            if e == 2: k = (k + n) % (2 * n); e = 0
            return 2 * k + e
        return from_presentation(4 * n, m)
    def semidirect_cyc(n, m_, r):   # C_n : C_m, b a b^-1 = a^r ; elements a^k b^l
        def m(i, j):
            k1, l1 = divmod(i, m_); k2, l2 = divmod(j, m_)
            return ((k1 + k2 * pow(r, l1, n)) % n) * m_ + (l1 + l2) % m_
        return from_presentation(n * m_, m)
    def direct(A, B):
        na, nb = len(A[0]), len(B[0])
        return [tuple(list(a) + [na + x for x in b]) for a in A for b in B]
    Dic3 = dicyclic(3); Dic6 = dicyclic(6); C3C8 = semidirect_cyc(3, 8, 2)
    C2 = [(0, 1), (1, 0)]
    D4 = closure([perm(4, (0, 1, 2, 3)), perm(4, (0, 2))]); C3 = closure([perm(3, (0, 1, 2))])
    C3D4 = None
    # C3 : D4 (= (C6 x C2) : C2, SmallGroup(24,8)) : D4 = <r, s>, r agit trivialement sur C3, s par inversion
    def c3d4():
        # C3 : D4 = SmallGroup(24,8) : D4 agit sur C3 par un caractere de noyau le groupe de Klein {1, r^2, s, s r^2} (avec le noyau <r> on obtiendrait D12)
        els = [(k, g) for k in range(3) for g in D4]
        klein = set(closure([perm(4, (0, 2), (1, 3)), perm(4, (0, 2))]))
        chi = {g: (1 if g in klein else -1) for g in D4}
        index = {x: i for i, x in enumerate(els)}
        return [tuple(index[((k1 + chi[g1] * k2) % 3, mul(g1, g2))] for (k2, g2) in els) for (k1, g1) in els]
    groups_full = {k: closure(v) for k, v in groups.items()}
    groups_full['Dic3 x C2'] = closure(direct(Dic3, C2)); groups_full['Dic6'] = closure(Dic6)
    groups_full['C3 : C8'] = closure(C3C8); groups_full['C3 : D4'] = closure(c3d4())
    print('\nGroupes d\'ordre 24 (15 attendus) : nombre de classes de conjugaison et couples (G, eta) donnant C + H + M_3(C)')
    assert len(groups_full) == 15
    found = []
    for name, G in sorted(groups_full.items()):
        assert len(G) == 24, (name, len(G))
        k = nclasses(G); line = f'  {name:14s} classes = {k:2d}'
        if k == 5:
            for eta in characters_pm1(G):
                r = wedderburn(G, eta)[0]
                names = sorted(n for n, *_ in r)
                line += f' ; eta {"trivial" if all(v == 1 for v in eta.values()) else "non trivial"} -> {" + ".join(names)}'
                if names == ['C', 'H', 'M_3(C)']: found.append(name)
        print(line)
    def order(g):
        e = tuple(range(len(g))); k = 1; h = g
        while h != e: h = mul(h, g); k += 1
        return k
    sig = [(nclasses(G), tuple(sorted(order(g) for g in G))) for G in groups_full.values()]
    print('  les 15 groupes sont deux a deux non isomorphes (classes + statistique des ordres distinctes) :', len(set(sig)) == 15)
    print('  => seul couple donnant A_F :', found, '  ', 'CONFORME' if found == ['S4'] else 'ECART')
