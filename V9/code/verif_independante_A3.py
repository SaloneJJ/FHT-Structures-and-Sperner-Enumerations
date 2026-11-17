# -*- coding: utf-8 -*-
"""Verification independante de trois resultats de A3 (point 1, 16/11).
Ecrit a partir des seules definitions de l'article (V9), sans importer ni lire
le code V352/Fable. numpy seul.
 (a) Thm P7:thm:bbh : reseau octupolaire construit en ESPACE REEL (sites, sauts,
     jauge de flux pi), spectre compare a E^2 = sum_d (tf^2+tc^2+2 tf tc cos k_d),
     degenerescence 4, masse m = sqrt3 |tf-tc| en K.
 (b) Prop A3:prop:heat, niveau plat seul : n_eff(40)=8.04, d_s(10)=3.04,
     analytiquement (Bessel) et numeriquement sur le reseau de (a).
 (c) Prop A3:prop:heat, arbre N=2 (m=8) : D_s construit depuis la formule de
     la sous-section Setting (T_d, V de poids 1/sqrt8, X_d, Y_d, Gamma_4,
     Gamma_5, F), moyenne de Bloch exacte sur k0 in (2pi Z/m)^3.
Usage : python3 verif_independante_A3.py [a|b|c|all]
"""
import sys, json, itertools, numpy as np

# ---------- generateurs de Cl(6) sur C^8 (Jordan-Wigner = jauge de flux pi) ----------
s0 = np.eye(2); sx = np.array([[0,1],[1,0]],complex); sy = np.array([[0,-1j],[1j,0]]); sz = np.diag([1.,-1.]).astype(complex)
def k3(a,b,c): return np.kron(a,np.kron(b,c))
X = [k3(sx,s0,s0), k3(sz,sx,s0), k3(sz,sz,sx)]     # saut intra-cellule le long de d, signe (-1)^{s_<d}
Y = [k3(sy,s0,s0), k3(sz,sy,s0), k3(sz,sz,sy)]
I8 = np.eye(8)

# ======================= (a) ================================================
def bbh_reel(L, tf, tc):
    """Reseau cubique (2L)^3 periodique, flux pi par plaquette :
    sauts x : +1 ; y : (-1)^x ; z : (-1)^(x+y). Amplitude tf si le lien est
    interne a la cellule 2x2x2 (coordonnee paire -> impaire), tc sinon."""
    n = 2*L; idx = lambda x,y,z: ((x%n)*n + (y%n))*n + (z%n)
    H = np.zeros((n**3, n**3))
    for x,y,z in itertools.product(range(n), repeat=3):
        for d,(dx,dy,dz) in enumerate([(1,0,0),(0,1,0),(0,0,1)]):
            c = (x,y,z)[d]; amp = tf if c%2==0 else tc
            sg = [1, (-1)**x, (-1)**(x+y)][d]
            i, j = idx(x,y,z), idx(x+dx,y+dy,z+dz)
            H[i,j] += sg*amp; H[j,i] += sg*amp
    return H

def test_a():
    out = {}
    for (tf,tc) in [(1.0,1.0),(1.3,0.7),(0.6,1.1)]:
        L = 4
        E = np.sort(np.linalg.eigvalsh(bbh_reel(L,tf,tc)))
        ks = 2*np.pi*np.arange(L)/L
        Ef = []
        for k in itertools.product(ks, repeat=3):
            e = np.sqrt(sum(tf**2+tc**2+2*tf*tc*np.cos(kd) for kd in k))
            Ef += [e]*4 + [-e]*4
        Ef = np.sort(np.array(Ef))
        err = float(np.max(np.abs(E-Ef)))
        # masse en K=(pi,pi,pi) : L pair contient k=pi
        mK = float(np.sqrt(3)*abs(tf-tc))
        EK = np.sqrt(sum(tf**2+tc**2+2*tf*tc*np.cos(np.pi) for _ in range(3)))
        out[f"tf={tf},tc={tc}"] = dict(ecart_max_spectre=err, m_K_formule=mK, E_K_reseau=float(EK))
        print(f"(a) tf={tf} tc={tc} : ecart max spectre reel / formule = {err:.2e} ; |E(K)| = {EK:.6f} ; sqrt3|tf-tc| = {mK:.6f}")
    return out

# ======================= (b) ================================================
def I0(x): return np.i0(x)
def z_plat_analytique(t):
    # H^2 = 4 sum sin^2(k_d/2) 1_8 ; <exp(-4t sin^2(k/2))>_k = e^{-2t} I0(2t)
    return 8.0*(np.exp(-2*t)*I0(2*t))**3
def neff_ds(zfun, t, h=1e-4):
    z = zfun(t); lp, lm = np.log(zfun(t*(1+h))), np.log(zfun(t*(1-h)))
    ds = -2*(lp-lm)/(np.log(1+h)-np.log(1-h))
    return float(z*(4*np.pi*t)**1.5), float(ds)

def test_b():
    out = {}
    for t in [10, 20, 40]:
        ne, ds = neff_ds(z_plat_analytique, t)
        out[f"analytique t={t}"] = dict(n_eff=ne, d_s=ds)
        print(f"(b) analytique  t={t:>3} : n_eff = {ne:.4f}  d_s = {ds:.4f}")
    # numerique : reseau de (a), tf=tc=1, jauge alternee (tc -> -tc), L=24 en Bloch
    L = 64; ks = 2*np.pi*np.arange(L)/L
    E2 = []
    for k in itertools.product(ks, repeat=3):
        E2.append(sum(2+2*np.cos(kd + np.pi) for kd in k))  # tc -> -tc : cos k -> -cos k
    E2 = np.array(E2)
    zf = lambda t: 8*np.mean(np.exp(-t*E2))
    for t in [10, 40]:
        ne, ds = neff_ds(zf, t)
        out[f"Bloch L={L} t={t}"] = dict(n_eff=ne, d_s=ds)
        print(f"(b) Bloch L={L} t={t:>3} : n_eff = {ne:.4f}  d_s = {ds:.4f}")
    return out

# ======================= (c) ================================================
G5 = (X[0]+X[1]+X[2])/np.sqrt(3)
G4 = G5 @ Y[0] @ Y[1] @ Y[2]
A_ = (X[0]-X[1])/np.sqrt(2); B_ = (X[0]+X[1]-2*X[2])/np.sqrt(6)
F = 1j*A_ @ B_
FG5 = F @ G5

def arbre_blocs(N):
    """indexation des cellules d'une racine : couche n a (2^n)^3 cellules"""
    off = [0]; 
    for n in range(N+1): off.append(off[-1] + 8**n)
    return off

def D_tree(N, k0, s, tf=1., tc=1., r=1., mu=1., alt=True):
    tcs = -tc if alt else tc
    off = arbre_blocs(N); nc = off[-1]
    T = [np.zeros((nc,nc),complex) for _ in range(3)]
    V = np.zeros((nc,nc))
    for n in range(N+1):
        S = 2**n
        cid = lambda x: off[n] + (x[0]*S + x[1])*S + x[2]
        for x in itertools.product(range(S), repeat=3):
            for d in range(3):
                y = list(x); y[d] += 1; ph = 1.0
                if y[d] == S: y[d] = 0; ph = np.exp(1j*k0[d])
                T[d][cid(x), cid(tuple(y))] += ph       # (T f)(x) = f(x+e_d)
            if n < N:
                S2 = 2*S
                for sig in itertools.product((0,1), repeat=3):
                    c = tuple(2*x[i]+sig[i] for i in range(3))
                    V[off[n+1] + (c[0]*S2 + c[1])*S2 + c[2], cid(x)] = 1/np.sqrt(8)  # parent -> enfant
    Inc = np.eye(nc)
    D = np.zeros((8*nc,8*nc),complex)
    for d in range(3):
        Td = T[d]
        D += tcs*np.kron((Td-Td.conj().T)/2j, Y[d])
        D += s*(tf*np.kron(Inc, X[d]) + tcs*np.kron((Td+Td.conj().T)/2, X[d]))
    D += np.kron((V-V.T)/2j, G4)
    D += s*r/2*np.kron(2*Inc - V - V.T, G5)
    D += mu*np.kron(Inc, FG5)
    return D

def test_c(N=2, m=8, ts=(10,15,20,30,40,60)):
    ks = 2*np.pi*np.arange(m)/m
    res = {}
    for s in (+1,-1):
        Es = []
        for k0 in itertools.product(ks, repeat=3):
            D = D_tree(N, k0, s)
            assert np.allclose(D, D.conj().T)
            Es.append(np.linalg.eigvalsh(D))
        Es = np.concatenate(Es); E2 = Es**2
        nzero = int(np.sum(np.abs(Es) < 1e-7))
        ncell = m**3 * 8**N
        def zf(t, drop=0):
            return (np.sum(np.exp(-t*E2)) - drop)/ncell
        res[s] = dict(n_zero=nzero, E_min_nonzero=float(np.min(np.abs(Es[np.abs(Es)>1e-7]))))
        print(f"(c) N={N} m={m} s={s:+d} : {len(Es)} valeurs propres, {nzero} zeros (|E|<1e-7)")
        for t in ts:
            ne_raw, ds_raw = neff_ds(lambda tt: zf(tt), t)
            ne_c, ds_c = neff_ds(lambda tt: zf(tt, drop=nzero), t)
            res[s][f"t={t}"] = dict(n_eff_brut=ne_raw, d_s_brut=ds_raw, n_eff_sans_zeros=ne_c, d_s_sans_zeros=ds_c)
            print(f"     t={t:>3} : n_eff brut {ne_raw:.4f}  sans zeros {ne_c:.4f}   d_s brut {ds_raw:.3f}  sans zeros {ds_c:.3f}")
    return res

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = {}
    if what in ("a","all"): out["a"] = test_a()
    if what in ("b","all"): out["b"] = test_b()
    if what in ("c","all"): out["c"] = {str(k):v for k,v in test_c().items()}
    with open(f"verif_{what}.json","w") as f: json.dump(out, f, indent=1)
