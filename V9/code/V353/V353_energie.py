"""V353 piste (d) : énergie du vide fermionique et action spectrale finie en fonction de φ. Critères : 01_synthese/V353 §1 (d1–d3), écrits avant.
(1) réseau périodique (mode « unif » de V336 : t_f = 1, t_c = −1, γ_F habille les termes « cos », μ FΓ5 nu) ⊗ fibre jouet (γ_F = σz, D_F = σx),
    insertion ε_D O_φ ⊗ σx ; grille k NK³.  (2) arbre N = 2 sans flux (Hnat de V336, θ = 0,3), diagonalisation complète.
Usage : python3 V353_energie.py [NK]     (l'arbre exige V336_FHT_flot.py, V162_FHT_R3_jauge.py, rvA_common.py dans le dossier courant)   JJ SALONE, FHT 2026"""
import sys, numpy as np, itertools as it
NK = int(sys.argv[1]) if len(sys.argv) > 1 else 16
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
X = [k3(sx, I2, I2), k3(sz, sx, I2), k3(sz, sz, sx)]; Y = [k3(sy, I2, I2), k3(sz, sy, I2), k3(sz, sz, sy)]
G5 = sum(X)/np.sqrt(3); A = (X[0]-X[1])/np.sqrt(2); B = (X[0]+X[1]-2*X[2])/np.sqrt(6); Y123 = Y[0]@Y[1]@Y[2]; FG = 1j*A@B@G5
def O(phi): return (np.cos(phi)*A + np.sin(phi)*B)@Y123
ks = 2*np.pi*(np.arange(NK) + 0.5)/NK                       # grille décalée (pas de point fixe particulier)
KK = np.array(list(it.product(ks, ks, ks)))
def Hs(s, mu):
    c = np.cos(KK); sn = np.sin(KK)
    H = np.zeros((len(KK), 8, 8), complex)
    for d in range(3):
        H += -sn[:, d, None, None]*Y[d] + s*(1 - c[:, d, None, None])*X[d]
    return H + mu*FG
def Hfull(phi, eps, mu):
    H = np.zeros((len(KK), 16, 16), complex)
    H[:, :8, :8] = Hs(+1, mu); H[:, 8:, 8:] = Hs(-1, mu); H[:, :8, 8:] = eps*O(phi); H[:, 8:, :8] = eps*O(phi)
    return H
def evac(phi, eps, mu):
    E = np.linalg.eigvalsh(Hfull(phi, eps, mu)); return E[E < 0].sum(), E
PH = np.deg2rad(np.arange(0, 180, 5))
def harm(vals):
    v = np.array(vals) - np.mean(vals); out = {}
    for n in (2, 4, 6, 12):
        out[n] = (2*np.mean(v*np.cos(n*PH)), 2*np.mean(v*np.sin(n*PH)))
    return out
SKIP = len(sys.argv) > 2 and sys.argv[2] == "arbre"
print(f"=== (1) réseau périodique, grille {NK}³ = {len(KK)} points k")
res = {}
for mu in (() if SKIP else (0.5, 1.0, 1.5)):
    for eps in (0.05, 0.1, 0.2, 0.4):
        vals = []; spec = {}
        for p in PH:
            e, E = evac(p, eps, mu); vals.append(e); spec[round(np.rad2deg(p))] = np.sort(E.ravel())
        h = harm(vals); res[(mu, eps)] = (vals, h)
        imin = int(np.argmin(vals)); amp6 = np.hypot(*h[6]); ph6 = np.rad2deg(np.arctan2(h[6][1], h[6][0]))
        print(f"  μ={mu} ε_D={eps:4.2f} : E_vac moyen {np.mean(vals):.6f} ; ΔE max−min {max(vals)-min(vals):.3e} ; |c2| {np.hypot(*h[2]):.1e} |c4| {np.hypot(*h[4]):.1e} "
              f"|c6| {amp6:.3e} (phase {ph6:+.2f}°) |c12| {np.hypot(*h[12]):.1e} ; minimum à φ = {np.rad2deg(PH[imin]):.0f}° (grille 5°)", flush=True)
        if eps == 0.2:
            d = lambda a, b: np.abs(spec[a] - spec[b]).max()
            print(f"      spectre complet : |S(φ)−S(φ+60°)| {d(10, 70):.1e} ; |S(φ)−S(−φ)| (10° vs 170°) {d(10, 170):.1e} ; |S(10°)−S(40°)| {d(10, 40):.1e}")
print("  exposant p (|c6| ∝ ε_D^p), par μ :")
for mu in (() if SKIP else (0.5, 1.0, 1.5)):
    e = np.array([0.05, 0.1, 0.2]); a = np.array([np.hypot(*res[(mu, x)][1][6]) for x in e])
    print(f"    μ={mu} : p(0,05→0,1) = {np.log(a[1]/a[0])/np.log(2):.2f} ; p(0,1→0,2) = {np.log(a[2]/a[1])/np.log(2):.2f} ; phase de c6 : "
          + " ".join(f"{np.rad2deg(np.arctan2(*res[(mu, x)][1][6][::-1])):+.2f}°" for x in (0.05, 0.1, 0.2, 0.4)))
print("=== (d3) action spectrale finie : Σ_k Tr H^{2n} en fonction de φ (μ = 1, ε_D = 0,2) : écart relatif max sur φ")
for n in (range(0) if SKIP else range(1, 8)):
    tr = []
    for p in PH[::3]:
        H = Hfull(p, 0.2, 1.0); M = np.linalg.matrix_power(H, 2*n); tr.append(np.real(np.trace(M, axis1=1, axis2=2)).sum())
    tr = np.array(tr); print(f"  2n = {2*n:2d} : Σ Tr H^2n = {tr.mean():.6e} ; (max−min)/moyenne = {(tr.max()-tr.min())/abs(tr.mean()):.1e}")
# (2) arbre N = 2 sans flux
try:
    sys.argv = ['x', '2', '0', '1', '1', '4', '0', 'unif', '1.0']
    src = open('V336_FHT_flot.py').read().split("t0 = time.time(); H0")[0]
    exec(src)
    # CORRECTION (13/11, 14 h 25) : V162 redéfinit B (champ magnétique scalaire) ; on reprend A, B de l'octet (Am, Bm de V336)
    def O(phi): return (np.cos(phi)*Am + np.sin(phi)*Bm)@(Yg[0]@Yg[1]@Yg[2])
    import scipy.sparse as sp
    print(f"=== (2) arbre N = 2 (Nc = {Nc}), sans flux, θ = 0,3, μ = 1, R-unif ; dimension {16*Nc}")
    Hp, Hm = Hnat(0.3, s=+1).toarray(), Hnat(0.3, s=-1).toarray(); Id = np.eye(Nc)
    for eps in (0.05, 0.1, 0.2, 0.4):
        vals = []; SP = {}
        for p in PH:
            Hx = eps*np.kron(Id, O(p)); H = np.block([[Hp, Hx], [Hx, Hm]]); E = np.linalg.eigvalsh(H); vals.append(E[E < 0].sum()); SP[round(np.rad2deg(p))] = E
        h = harm(vals); imin = int(np.argmin(vals))
        print(f"  ε_D={eps:4.2f} : ΔE max−min {max(vals)-min(vals):.3e} ; |c2| {np.hypot(*h[2]):.1e} |c4| {np.hypot(*h[4]):.1e} |c6| {np.hypot(*h[6]):.3e} "
              f"(phase {np.rad2deg(np.arctan2(h[6][1], h[6][0])):+.2f}°) ; minimum à φ = {np.rad2deg(PH[imin]):.0f}°", flush=True)
        print(f"      |S(10°)−S(130°)| {np.abs(SP[10]-SP[130]).max():.1e} ; |S(10°)−S(170°)| {np.abs(SP[10]-SP[170]).max():.1e} ; |S(10°)−S(40°)| {np.abs(SP[10]-SP[40]).max():.1e}")
except FileNotFoundError as e:
    print("  arbre : fichiers de V336 absents :", e)
