"""V339 relecture, partie 3 : V337 S1 — la jauge alternée (t_c -> -t_c, poids vertical 1) est-elle équivalente à (jauge K, ω = σ(enfant)σ(parent)) ?
Transformation G = ⊕_n diag σ(x) ⊗ 1_8 ; comparaison des matrices (où est l'écart ?) et des spectres ; contrôle J (structure réelle) sur l'arbre à champ nul."""
import sys, os, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
here = os.path.dirname(os.path.abspath(__file__))
def model(N, ga, om, mt=0.0, mu=1.0, s=1, th=0.3):
    ns = {'__name__': 'rl2', '__file__': os.path.join(here, 'rl2_flot.py')}
    sys.argv = ['rl2', str(N), '8', ga, om, str(mt), str(mu), str(s)]
    exec(open(os.path.join(here, 'rl2_flot.py')).read(), ns)
    return ns
for N in (2, 3):
    a = model(N, 'alt', '1'); k = model(N, 'K', 'ScSp')
    Nc, sides, grid, depth = a['Nc'], a['sides'], a['grid'], a['depth']
    sg = np.concatenate([((-1.0) ** sum(grid(s))).ravel() for s in sides])
    G = sp.kron(sp.diags(sg), sp.identity(8), format='csr')
    for th in (0.3, 2.0):
        Ha = a['Hs'](th); Hk_ = k['Hs'](th)
        Dm = (G @ Ha @ G - Hk_).tocoo(); big = np.abs(Dm.data) > 1e-12
        rows = np.unique(Dm.row[big] // 8)
        print(f'N={N} θ={th} : |G H_alt G − H_K,σσ′| max = {np.abs(Dm.data).max() if Dm.nnz else 0:.2e} ; cellules concernées : {rows.tolist()} (profondeurs {sorted(set(depth[rows].tolist()))})')
        if N == 2:
            Ea = np.linalg.eigvalsh(Ha.toarray()); Ek = np.linalg.eigvalsh(Hk_.toarray())
            print(f'   spectres complets : écart max {np.abs(Ea - Ek).max():.3e} ; 6 plus petits |E| alt {np.round(np.sort(np.abs(Ea))[:6], 5)} ; K,σσ′ {np.round(np.sort(np.abs(Ek))[:6], 5)}')
        else:
            Ea = sla.eigsh(Ha, k=24, sigma=0, which='LM')[0]; Ek = sla.eigsh(Hk_, k=24, sigma=0, which='LM')[0]
            print(f'   24 valeurs propres près de 0 : écart max {np.abs(np.sort(Ea) - np.sort(Ek)).max():.3e}')
    # racine corrigée : même modèle K,σσ′ mais t_c -> -t_c à la racine seulement
# structure réelle sur l'arbre (champ nul : θ = 0 et pas de flux -> on prend q = 0 en reconstruisant avec q = 0)
ns = {'__name__': 'rl2', '__file__': os.path.join(here, 'rl2_flot.py')}; sys.argv = ['rl2', '2', '8', 'K', 'ScSp', '0', '1', '1']
src = open(os.path.join(here, 'rl2_flot.py')).read().replace("'1', 'none', '1', str(NTH)", "'0', 'none', '1', str(NTH)")
exec(src, ns)
Nc = ns['Nc']; U8 = ns['X5'] @ ns['A_'] @ ns['B_']
Ut = sp.kron(sp.identity(Nc), U8, format='csr')
Hp = ns['Hs'](0.0, 1); Hm = ns['Hs'](0.0, -1)
print('champ nul (q = 0, θ = 0), ω = σσ′ : |U conj(H_+) U⁻¹ + H_−| =', f"{abs(Ut @ Hp.conj() @ Ut.getH() + Hm).max():.1e}", '(J_lat envoie le secteur γ_F = +1 sur −(secteur −1) : relation KO-1 de V331 conservée sur l\'arbre)')
