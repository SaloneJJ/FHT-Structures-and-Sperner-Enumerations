# Complement (c) : ou vivent les modes zero ? n_eff/d_s en retirant seulement ceux de la racine.
import itertools, json, numpy as np
from verif_independante_A3 import D_tree, arbre_blocs, neff_ds
N, m = 2, 8; off = arbre_blocs(N); ks = 2*np.pi*np.arange(m)/m
out = {}
for s in (+1, -1):
    E2all = []; zeros = []
    for k0 in itertools.product(ks, repeat=3):
        E, U = np.linalg.eigh(D_tree(N, k0, s))
        E2all.append(E**2)
        for j in np.where(np.abs(E) < 1e-7)[0]:
            w = np.abs(U[:, j])**2
            poids = [float(w[8*off[n]:8*off[n+1]].sum()) for n in range(N+1)]
            zeros.append(dict(k0=[float(x) for x in k0], poids_couches=poids))
    E2 = np.concatenate(E2all); ncell = m**3 * 8**N
    nroot = sum(1 for zz in zeros if zz["poids_couches"][N] < 0.5)
    print(f"s={s:+d} modes zero :")
    for zz in zeros: print("   k0 =", zz["k0"], " poids par couche (0..N) =", [round(p,4) for p in zz["poids_couches"]])
    rows = {}
    for t in (10, 15, 20, 30, 40, 60):
        ne, ds = neff_ds(lambda tt: (np.sum(np.exp(-tt*E2)) - nroot)/ncell, t)
        rows[t] = (ne, ds)
        print(f"   t={t:>3} : n_eff (sans les {nroot} zeros de la racine) = {ne:.4f}   d_s = {ds:.3f}")
    out[str(s)] = dict(zeros=zeros, table={str(t): dict(n_eff=v[0], d_s=v[1]) for t, v in rows.items()})
json.dump(out, open("verif_c_modes_zero.json", "w"), indent=1)
