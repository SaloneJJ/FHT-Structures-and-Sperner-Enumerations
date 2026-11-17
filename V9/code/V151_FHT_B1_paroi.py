"""V151 B1 : fermions de paroi (Wilson-Dirac 4D) dans l'arbre augmente ferme en T^3, et controle plat.
Criteres et predictions : V151 par. 2 (ecrits avant calcul).
Usage : python3 V151_FHT_B1_paroi.py N m [tree|flat|flatk] [K]
  tree  : arbre (profondeur n : 2^n cellules par cote), reduction par les 8 secteurs (Z/2)^3
  flat  : meme code en espace reel, couches de 2^N cellules (controle de coherence)
  flatk : controle plat exact par impulsion laterale
JJ SALONE, FHT 2026"""
import sys, time, itertools
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
t0 = time.time()
np.set_printoptions(precision=4, suppress=True, linewidth=160)
N = int(sys.argv[1]); m = float(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'tree'
K = int(sys.argv[4]) if len(sys.argv) > 4 else 40
NPRINT = int(sys.argv[5]) if len(sys.argv) > 5 else 140
SECT = sys.argv[6].split(',') if len(sys.argv) > 6 else None
r = 1.0
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2, dtype=complex)
GAM = [np.kron(sx, sx), np.kron(sx, sy), np.kron(sx, sz), np.kron(sy, I2), np.kron(sz, I2)]   # Gamma_1..Gamma_5
CHI = np.kron(sx, I2)                                                                       # sigma_x (x) 1 = -i Gamma_4 Gamma_5
for a in range(5):
    for b in range(5):
        assert np.allclose(GAM[a] @ GAM[b] + GAM[b] @ GAM[a], 2 * (a == b) * np.eye(4))
assert np.allclose(CHI, -1j * GAM[3] @ GAM[4])
print(f'=== B1 : N = {N}, m = {m}, r = {r}, mode = {mode} ===', flush=True)

def geometry(flat):
    sides = [2 ** N if flat else 2 ** n for n in range(N + 1)]
    offs = np.concatenate([[0], np.cumsum([s ** 3 for s in sides])]).astype(np.int64)
    return sides, offs

def grid(s):
    I, J, Kk = np.meshgrid(np.arange(s), np.arange(s), np.arange(s), indexing='ij')
    return I.ravel(), J.ravel(), Kk.ravel()

def operators(flat):
    sides, offs = geometry(flat); Nc = int(offs[-1])
    rowsT = [[], [], []]; colsT = [[], [], []]; rS = []; cS = []; vS = []
    for n, s in enumerate(sides):
        I, J, Kk = grid(s); base = offs[n] + (I * s + J) * s + Kk
        for d, (a, b, c) in enumerate([((I + 1) % s, J, Kk), (I, (J + 1) % s, Kk), (I, J, (Kk + 1) % s)]):
            rowsT[d].append(base); colsT[d].append(offs[n] + (a * s + b) * s + c)
        if n >= 1:
            sp_ = sides[n - 1]
            if flat: pi, pj, pk, val = I, J, Kk, 1.0
            else:    pi, pj, pk, val = I // 2, J // 2, Kk // 2, 1 / np.sqrt(8)
            rS.append(base); cS.append(offs[n - 1] + (pi * sp_ + pj) * sp_ + pk); vS.append(np.full(len(base), val))
    T = [sp.coo_matrix((np.ones(Nc), (np.concatenate(rowsT[d]), np.concatenate(colsT[d]))), shape=(Nc, Nc)).tocsr() for d in range(3)]
    S = sp.coo_matrix((np.concatenate(vS), (np.concatenate(rS), np.concatenate(cS))), shape=(Nc, Nc)).tocsr()
    depth = np.concatenate([np.full(s ** 3, n) for n, s in enumerate(sides)])
    return sides, offs, Nc, T, S, depth

def hamiltonian(T, S, Nc, mass):
    kr = lambda A, g: sp.kron(A, sp.csr_matrix(g), format='csr')
    H = sp.csr_matrix((4 * Nc, 4 * Nc), dtype=complex)
    for d in range(3):
        H = H + kr((T[d] - T[d].T) / 2j, GAM[d]) - (r / 2) * kr(T[d] + T[d].T, GAM[4])
    H = H + kr((S - S.T) / 2j, GAM[3]) - (r / 2) * kr(S + S.T, GAM[4]) + (mass + 4 * r) * kr(sp.identity(Nc), GAM[4])
    return H.tocsr()

def half_shift(sides, offs, Nc, d, flat):
    rows, cols = [], []
    for n, s in enumerate(sides):
        I, J, Kk = grid(s); base = offs[n] + (I * s + J) * s + Kk
        h = s // 2 if s >= 2 else 0
        a = [I, J, Kk]; a[d] = (a[d] + h) % s
        rows.append(base); cols.append(offs[n] + (a[0] * s + a[1]) * s + a[2])
    return sp.coo_matrix((np.ones(Nc), (np.concatenate(rows), np.concatenate(cols))), shape=(Nc, Nc)).tocsr()

def sector_basis(sides, offs, Nc, eps):
    rows, cols, vals = [], [], []; col = 0
    for n, s in enumerate(sides):
        if s == 1:
            if eps == (0, 0, 0):
                rows.append([offs[n]]); cols.append([col]); vals.append([1.0]); col += 1
            continue
        h = s // 2; I, J, Kk = grid(h); nrep = len(I)
        for dl in itertools.product((0, 1), repeat=3):
            sgn = (-1) ** (eps[0] * dl[0] + eps[1] * dl[1] + eps[2] * dl[2]) / np.sqrt(8)
            rows.append(offs[n] + ((I + dl[0] * h) * s + (J + dl[1] * h)) * s + (Kk + dl[2] * h))
            cols.append(col + np.arange(nrep)); vals.append(np.full(nrep, sgn))
        col += nrep
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(Nc, col)).tocsr()

def analyse(psi, sides, offs, depth, Ltop):
    p = psi.reshape(-1, 4); w = np.sum(np.abs(p) ** 2, axis=1)
    wn = np.bincount(depth, weights=w, minlength=N + 1); nbar = float(np.dot(np.arange(N + 1), wn))
    chi_cell = np.real(np.einsum('ia,ab,ib->i', p.conj(), CHI, p))
    chig = float(np.sum(chi_cell))
    def chi_at(n):
        sl = slice(offs[n], offs[n + 1]); ww = w[sl].sum()
        return float(chi_cell[sl].sum() / ww) if ww > 1e-8 else np.nan
    top = p[offs[N]:offs[N + 1]].reshape(Ltop, Ltop, Ltop, 4)
    F = np.fft.fftn(top, axes=(0, 1, 2)); P = np.sum(np.abs(F) ** 2, axis=3)
    fr = np.fft.fftfreq(Ltop, d=1.0 / Ltop).astype(int)
    n2 = fr[:, None, None] ** 2 + fr[None, :, None] ** 2 + fr[None, None, :] ** 2
    tot = P.sum()
    if tot > 1e-12:
        sh = np.bincount(n2.ravel(), weights=P.ravel()); s_dom = int(np.argmax(sh)); frac = float(sh[s_dom] / tot)
    else:
        s_dom, frac = -1, 0.0
    return wn, nbar, chig, chi_at(N), chi_at(0), s_dom, frac

def classe(nbar):
    return 'bord' if nbar >= N - 1.5 else ('racine' if nbar <= 1.0 else 'interieur')

def report(rows, Ltop, label):
    rows.sort(key=lambda x: abs(x[0]))
    print(f'\n  {label} : {len(rows)} etats ; colonnes : E, L*E/2pi, secteur, n_moyen, w(N), w(0), chi_global, chi(N), chi(0), |n|^2 dominant (part), classe', flush=True)
    for (E, eps, wn, nbar, chig, chiN, chi0, sdom, frac) in rows[:NPRINT]:
        print(f'   {E:+.5f} {Ltop*E/(2*np.pi):+7.3f}  {"".join(map(str,eps))}  {nbar:5.2f}  {wn[-1]:.3f} {wn[0]:.3f}  {chig:+.3f}  {chiN:+.3f}  {chi0:+.3f}   {sdom:3d} ({frac:.2f})  {classe(nbar):9s} w(n)=[{" ".join(f"{x:.2f}" for x in wn)}]')
    return rows

def summary(rows, Ltop, Ewin):
    print(f'\n  Bilan sous la fenetre |E| < {Ewin:.4f} (L*E/2pi < {Ltop*Ewin/(2*np.pi):.3f}) :')
    cnt = {}
    for (E, eps, wn, nbar, chig, chiN, chi0, sdom, frac) in rows:
        if abs(E) >= Ewin: continue
        key = (classe(nbar), '+' if chig > 0.3 else ('-' if chig < -0.3 else '0'))
        cnt[key] = cnt.get(key, 0) + 1
    for k_, v in sorted(cnt.items()): print(f'    classe {k_[0]:9s} helicite {k_[1]} : {v}')
    return cnt

Ltop = 2 ** N
shells = {s: np.sin(2 * np.pi * np.sqrt(s) / Ltop) for s in (1, 2, 3, 4)}
Ewin = 1.1 * np.sin(2 * np.pi * 2 / Ltop)      # juste au-dessus de la couche |n|^2 = 4 (estimation reseau plat, axe)

if mode == 'flatk':
    nl = N + 1; S1 = np.zeros((nl, nl)); S1[np.arange(1, nl), np.arange(0, nl - 1)] = 1
    Hv = np.kron((S1 - S1.T) / 2j, GAM[3]) - (r / 2) * np.kron(S1 + S1.T, GAM[4])
    fr = np.fft.fftfreq(Ltop, d=1.0 / Ltop).astype(int)
    ns = np.array(list(itertools.product(fr, repeat=3))); k = 2 * np.pi * ns / Ltop
    on = (np.einsum('kd,dab->kab', np.sin(k), np.array(GAM[:3])) + ((m + 4 * r - r * np.cos(k).sum(1)))[:, None, None] * GAM[4][None])
    Hk = np.einsum('ij,kab->kiajb', np.eye(nl), on).reshape(-1, 4 * nl, 4 * nl) + Hv[None]
    E, V = np.linalg.eigh(Hk)
    rows = []
    for ik in range(len(ns)):
        for a in range(4 * nl):
            if abs(E[ik, a]) > 1.6 * Ewin: continue
            p = V[ik, :, a].reshape(nl, 4); w = np.sum(np.abs(p) ** 2, 1); nbar = float(np.dot(np.arange(nl), w))
            cc = np.real(np.einsum('ia,ab,ib->i', p.conj(), CHI, p))
            chiN = cc[N] / w[N] if w[N] > 1e-8 else np.nan; chi0 = cc[0] / w[0] if w[0] > 1e-8 else np.nan
            s2 = int(np.sum(ns[ik] ** 2))
            rows.append((E[ik, a], tuple(int(x) % 2 for x in ns[ik]), w, nbar, float(cc.sum()), chiN, chi0, s2, 1.0))
    report(rows, Ltop, 'controle plat (exact par impulsion)'); summary(rows, Ltop, Ewin)
    print(f'\n  energies de reference (reseau plat, sin) : ' + ', '.join(f'|n|^2={s}: {e:.4f}' for s, e in shells.items()))
    print(f'temps {time.time()-t0:.1f} s'); sys.exit()

flat = (mode == 'flat')
sides, offs, Nc, T, S, depth = operators(flat)
H = hamiltonian(T, S, Nc, m)
print(f'  cellules : {Nc} ; dimension : {4*Nc} ; hermiticite : {abs(H - H.getH()).max():.1e}', flush=True)
for d in range(3):
    tau = sp.kron(half_shift(sides, offs, Nc, d, flat), sp.identity(4), format='csr')
    print(f'  [tau_{d}, H] = {abs(tau @ H - H @ tau).max():.1e}', end=';', flush=True)
print()
rows = []; dims = 0
for eps in itertools.product((0, 1), repeat=3):
    if SECT and ''.join(map(str, eps)) not in SECT: continue
    B = sp.kron(sector_basis(sides, offs, Nc, eps), sp.identity(4), format='csr')
    He = (B.getH() @ H @ B).tocsc(); dim = He.shape[0]; dims += dim
    if dim <= 3000:
        E, V = np.linalg.eigh(He.toarray()); sel = np.argsort(np.abs(E))[:K]; E, V = E[sel], V[:, sel]
    else:
        E, V = sla.eigsh(He, k=K, sigma=1.234e-5, which='LM')
    for a in range(len(E)):
        psi = B @ V[:, a]
        wn, nbar, chig, chiN, chi0, sdom, frac = analyse(psi, sides, offs, depth, Ltop)
        rows.append((E[a], eps, wn, nbar, chig, chiN, chi0, sdom, frac))
    print(f'  secteur {eps} : dim {dim}, {len(E)} valeurs propres, |E| <= {np.abs(E).max():.4f}  ({time.time()-t0:.0f} s)', flush=True)
print(f'  somme des dimensions des secteurs : {dims} (= {4*Nc} ?)')
report(rows, Ltop, f'arbre ({"plat reel" if flat else "hierarchique"})'); summary(rows, Ltop, Ewin)
neg = sorted([x for x in rows if x[4] < -0.3], key=lambda x: abs(x[0]))
print(f'\n  tous les etats calcules d helicite negative (chi_global < -0.3) : {len(neg)}')
for (E, eps, wn, nbar, chig, chiN, chi0, sdom, frac) in neg:
    print(f'   {E:+.5f}  {"".join(map(str,eps))}  n_moyen {nbar:5.2f}  chi {chig:+.3f}  |n|^2 au bord {sdom:3d}  w(n)=[{" ".join(f"{x:.2f}" for x in wn)}]')
print(f'\n  energies de reference (reseau plat, sin) : ' + ', '.join(f'|n|^2={s}: {e:.4f}' for s, e in shells.items()))
print(f'temps {time.time()-t0:.1f} s')
