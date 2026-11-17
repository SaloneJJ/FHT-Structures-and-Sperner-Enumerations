"""V357 (relecture de V352) : analyse indépendante des spectres de Bloch produits par V352_noyau.py (N = 2, m = 8).
z(t) = Σ_k w_k Σ_i e^{−t E_i²} / 8^N (par cellule fine) ; n_eff = z (4πt)^{3/2} ; d_s = −2 dln z / dln t.
Correction de V352 §3.3 : retrait des 2 modes nuls de la racine à k0 = 0 (profil vérifié ici, indépendamment).
Insertion : rapport z_ins / (2 z_0) comparé à e^{−t ε_D²}. Contrôle plat exact : 8 [e^{−2t} I0(2t)]^3 sur la même grille finie."""
import sys, os, numpy as np
sys.path[:0] = [os.path.join(os.path.dirname(os.path.abspath(__file__)), d) for d in ('v352run', '../V352')] + ['.']
import V352_arbre as M
d = sys.argv[1] if len(sys.argv) > 1 else '.'
def load(f):
    z = np.load(f, allow_pickle=True); return list(z['E']), z['w'], int(z['N']), int(z['m']), [tuple(r) for r in z['reps']]
T = np.array([5, 10, 15, 20, 30, 40, 60, 100.])
E0, w, N, m, reps = load(f'{d}/V352_N2_m8_s1_mu1.0_e0.0AY.npz')
# profil des modes nuls à k0 = 0 (contrôle indépendant de « 2 sur la racine, 2 sur la couche fine »)
H = M.Hnat(N, (0., 0., 0.), +1); E, V = np.linalg.eigh(H); z0 = np.where(np.abs(E) < 1e-8)[0]
dep = M.depth_of(N); P = (np.abs(V[:, z0]) ** 2).reshape(-1, 8, len(z0)).sum(1)
print(f'k0 = 0 : {len(z0)} modes nuls ; poids par niveau (somme sur les modes) : {[round(P[dep == n].sum(), 3) for n in range(N + 1)]}')
def heat(Es, w, N, T): return np.array([sum(wk * np.exp(-t * E ** 2).sum() for wk, E in zip(w, Es)) for t in T]) / 8 ** N
z = heat(E0, w, N, T)
i0 = reps.index((0, 0, 0)); Es2 = [e.copy() for e in E0]; e0 = Es2[i0]; Es2[i0] = np.delete(e0, np.where(np.abs(e0) < 1e-8)[0][:2])
zc = heat(Es2, w, N, T)
Lf = 2 ** N * m; kk = 2 * np.pi * np.arange(Lf) / Lf
plat = np.array([np.mean(np.exp(-t * 4 * np.sin(kk / 2) ** 2)) ** 3 for t in T])   # cône plat à 1 composante sur la même grille, par site
print(' t     n_eff brut   n_eff corrigé   n_grille (corrigé / plat)   d_s corrigé')
lz = np.log(zc); ds = -2 * np.gradient(lz, np.log(T))
for i, t in enumerate(T): print(f'{t:5.0f}  {z[i]*(4*np.pi*t)**1.5:9.4f}  {zc[i]*(4*np.pi*t)**1.5:9.4f}  {zc[i]/plat[i]:9.4f}  {ds[i]:7.3f}')
A = np.vstack([np.ones(4), 1 / T[1:6:1][:4], T[1:5]]).T
sel = (T >= 10) & (T <= 40); A = np.vstack([np.ones(sel.sum()), 1 / T[sel], T[sel]]).T
c = np.linalg.lstsq(A, (zc * (4 * np.pi * T) ** 1.5)[sel], rcond=None)[0]
print(f'ajustement n_eff = α + β/t + γt sur [10, 40] (N = 2) : α = {c[0]:.3f}, β = {c[1]:.3f}, γ = {c[2]:.2e}, |γ|·40/α = {abs(c[2])*40/c[0]:.3f}')
for O in ('AY', 'G5'):
    Ei, wi, *_ = load(f'{d}/V352_N2_m8_s1_mu1.0_e0.1{O}.npz')
    zi = heat(Ei, wi, N, T)
    print(f'insertion ε_D = 0,1, O = {O} : z_ins/(2 z0) = {np.round(zi / (2 * z), 5).tolist()} ; e^(−tε²) = {np.round(np.exp(-T * 0.01), 5).tolist()}')
