#!/usr/bin/env python3
"""AUDIT FABLE -- A3, Proposition « The flux is selected by the fermions » (registre C4 ; note V205, jamais relue avant cet audit).
Verification independante des valeurs de l'article pour les fermions libres a demi-remplissage, saut uniforme |t| = 1 :
  cube des enfants : -4 sqrt3 = -6.928 (flux pi sur les 6 faces), -6.472 au mieux avec quatre faces a pi, -6 sans flux ;
  amas ouverts 4^3 et 6^3 : -66.00 et -234.71 (flux pi), -55.96 et -197.61 (flux nul) ;
  reseau periodique 8^3 : -1.190 par site (flux pi) contre -0.989 (flux nul).
(Les comparaisons a 200 configurations Z/2 aleatoires et la coquille de K26 ne sont pas refaites ici.)
Usage : python3 AF6_flux_pi_energie.py      (numpy ; ~1 min)
"""
import itertools, numpy as np

def energy(H):
    w = np.linalg.eigvalsh(H); return w[w < 0].sum()

def cube_all():
    """toutes les configurations de signes sur les 12 aretes du cube ; energie minimale par nombre de faces a flux pi"""
    V = list(itertools.product((0, 1), repeat=3)); ix = {v: i for i, v in enumerate(V)}
    E = [(ix[v], ix[tuple(v[j] ^ (j == d) for j in range(3))], d) for v in V for d in range(3) if v[d] == 0]
    faces = []
    for d in range(3):
        for c in (0, 1):
            vs = [v for v in V if v[d] == c]; faces.append([k for k, (a, b, _) in enumerate(E) if V[a] in vs and V[b] in vs])
    best = {}
    for signs in itertools.product((1, -1), repeat=12):
        nf = sum(1 for f in faces if np.prod([signs[k] for k in f]) < 0)
        H = np.zeros((8, 8))
        for k, (a, b, _) in enumerate(E): H[a, b] = H[b, a] = signs[k]
        e = energy(H); best[nf] = min(best.get(nf, 0), e)
    return best

def lattice(L, flux_pi, periodic):
    idx = lambda x, y, z: (x % L) * L * L + (y % L) * L + (z % L)
    H = np.zeros((L**3, L**3))
    for x, y, z in itertools.product(range(L), repeat=3):
        for d, (dx, dy, dz) in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1))):
            if not periodic and (x + dx >= L or y + dy >= L or z + dz >= L): continue
            s = 1 if not flux_pi else (1, (-1) ** x, (-1) ** (x + y))[d]
            H[idx(x, y, z), idx(x + dx, y + dy, z + dz)] += s; H[idx(x + dx, y + dy, z + dz), idx(x, y, z)] += s
    return H

if __name__ == '__main__':
    ok = True
    b = cube_all()
    print('Cube des enfants : energie minimale selon le nombre de faces a flux pi :', {k: round(v, 3) for k, v in sorted(b.items())})
    for k, target in ((6, -4 * np.sqrt(3)), (4, -6.472), (0, -6.0)):
        good = abs(b[k] - target) < 1e-3; ok &= good
        print(f'   {k} faces a pi : {b[k]:.3f} (article : {target:.3f})  {"ok" if good else "ECART"}')
    gmin = min(b.values()); print(f'   minimum global sur les 2^12 configurations : {gmin:.3f}, atteint pour {[k for k, v in b.items() if abs(v - gmin) < 1e-9]} faces a pi', '  ok' if abs(gmin + 4 * np.sqrt(3)) < 1e-9 else '  ECART')
    ok &= abs(gmin + 4 * np.sqrt(3)) < 1e-9
    for L, (tp, t0) in ((4, (-66.00, -55.96)), (6, (-234.71, -197.61))):
        ep, e0 = energy(lattice(L, True, False)), energy(lattice(L, False, False))
        good = abs(ep - tp) < 6e-3 and abs(e0 - t0) < 6e-3; ok &= good
        print(f'Amas ouvert {L}^3 : flux pi {ep:.2f} (article {tp}) ; flux nul {e0:.2f} (article {t0})  {"ok" if good else "ECART"}')
    ep, e0 = energy(lattice(8, True, True)) / 512, energy(lattice(8, False, True)) / 512
    good = abs(ep + 1.190) < 1e-3 and abs(e0 + 0.989) < 1e-3; ok &= good
    print(f'Reseau periodique 8^3 : flux pi {ep:.3f} par site (article -1.190) ; flux nul {e0:.3f} (article -0.989)  {"ok" if good else "ECART"}')
    print('VERDICT :', 'conforme' if ok else 'ECART')
