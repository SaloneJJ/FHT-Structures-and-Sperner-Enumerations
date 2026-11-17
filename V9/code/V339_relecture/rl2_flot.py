"""V339 relecture, partie 2 : modèle natif réécrit + flot spectral par couche.
Usage : python3 rl2_flot.py N NTH jauge(K|alt) omega(1|Sc|Sp|ScSp) mtot mu s [conv(cos|v327)] [tf]
  jauge K : BBH dans sa jauge propre (cône en K) ; alt : t_c -> -t_c à chaque niveau (mode bbhG/unif de V335-V338).
  omega : poids du saut vertical enfant -> parent, σ = (-1)^(x+y+z) sur la grille de chaque niveau (Sc : enfant, Sp : parent).
  mtot : masse au cône (coefficient de Γ5 ajouté à √3(t_f - t_c)) ; mu : terme nu μ FΓ5 ; s = γ_F.
  conv : cos = toute la partie paire habillée (D-C3) ; v327 = seule la composante Γ5 de la partie paire latérale habillée (contrôle).
Construction : (Hψ)_x = H0 ψ_x + Σ_d M_d U ψ_{x+e_d} + h.c. ; partie paire/impaire = (T ± T^†) ⊗ (M ± M^T)/2 (parité en q autour du cône).
Géométrie, champ de fond (un quantum à travers (x,y), torsion θ selon z) et transport (b, r0) : V162 (exec). Flot : Hellmann–Feynman,
fenêtre carrée de demi-largeur ε = 0,25·2π/L (comme V166) ET fenêtres ε/2, 2ε et gaussienne (contrôle de robustesse)."""
import sys, time, os
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rl_common import *
a = sys.argv
N_, NTH, GA, OM, MT, MU, SG = int(a[1]), int(a[2]), a[3], a[4], float(a[5]), float(a[6]), int(a[7])
CONV = a[8] if len(a) > 8 else 'cos'; TF = float(a[9]) if len(a) > 9 else 1.0
sys.argv = ['x', str(N_), '-1', 'b', '1', 'none', '1', str(NTH), 'r0', '0/1']
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'V162_FHT_R3_jauge.py')).read().split("print(f'=== V162")[0])
H0o, Mo = octet_mats(TF, 1.0)
sg = -1.0 if GA == 'alt' else 1.0
Md = [sg * Mo[d] for d in range(3)]
X5 = sum(Mo[d] + Mo[d].T for d in range(3)) / np.sqrt(3)          # Γ5 (convention V335 : ΣX/√3)
Yd = [1j * (Mo[d] - Mo[d].T) for d in range(3)]
G4 = X5 @ Yd[0] @ Yd[1] @ Yd[2]                                  # orientation contrôlée dans rl1 : ε = −1 comme l'emprunté
Xd = [Mo[d] + Mo[d].T for d in range(3)]
A_ = (Xd[0] - Xd[1]) / np.sqrt(2); B_ = (Xd[0] + Xd[1] - 2 * Xd[2]) / np.sqrt(6); FG = 1j * A_ @ B_ @ X5
csr = sp.csr_matrix
def sig(I, J, K): return ((-1.0) ** (I + J + K)).ravel()
def tree_ops(th):
    U, TIN = link_fields(th, 0.0)
    T = []
    for d, sh in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
        R, C, V = [], [], []
        for n in range(N + 1):
            s = sides[n]; I, J, K = grid(s)
            R.append(lin(n, I, J, K).ravel()); C.append(lin(n, I + sh[0], J + sh[1], K + sh[2]).ravel()); V.append(np.conj(U[n][d]).ravel())
        T.append(csr((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(Nc, Nc)))
    R, C, V = [], [], []
    for n in range(1, N + 1):
        s = sides[n]; I, J, K = grid(s)
        u = np.zeros((s, s, s), complex)
        for c, P in TIN[n - 1].items(): u[c[0]::2, c[1]::2, c[2]::2] = P
        w = {'1': 1.0, 'Sc': sig(I, J, K), 'Sp': sig(I // 2, J // 2, K // 2), 'ScSp': sig(I, J, K) * sig(I // 2, J // 2, K // 2)}[OM]
        R.append(lin(n, I, J, K).ravel()); C.append(lin(n - 1, I // 2, J // 2, K // 2).ravel()); V.append(u.ravel() * w / np.sqrt(8))
    S = csr((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(Nc, Nc))
    return T, S
kr = lambda X, g: sp.kron(X, csr(g), format='csr')
def parts(th):
    T, S = tree_ops(th); Id = sp.identity(Nc, format='csr')
    odd = sum(kr(T[d] - T[d].getH(), (Md[d] - Md[d].T) / 2) for d in range(3)) + kr((S - S.getH()) / 2j, G4)
    ev_lat = kr(Id, H0o) + sum(kr(T[d] + T[d].getH(), (Md[d] + Md[d].T) / 2) for d in range(3))
    ev_v = kr(Id, MT * X5) + 0.5 * kr(2 * Id - S - S.getH(), X5)            # masse au cône + Wilson vertical (r = 1)
    return odd, ev_lat, ev_v
def Hs(th, s=SG):
    odd, ev_lat, ev_v = parts(th)
    if CONV == 'cos': H = odd + s * (ev_lat + ev_v)
    else:                                                               # v327 : seule la projection sur Γ5 (par cellule et par lien) habillée
        P5 = lambda Z: Z                                                # (contrôle grossier : on habille Γ5-part = (Z + Γ5 Z Γ5)/2 bloc par bloc)
        G = kr(sp.identity(Nc, format='csr'), X5)
        evp = (ev_lat + G @ ev_lat @ G) / 2; evr = ev_lat - evp
        H = odd + evr + s * (evp + ev_v)
    return (H + MU * kr(sp.identity(Nc, format='csr'), FG)).tocsr()
EPS = 0.25 * 2 * np.pi / L
def prof(v): return np.bincount(depth, weights=np.sum(np.abs(v.reshape(Nc, 8)) ** 2, axis=1), minlength=N + 1)
if __name__ == '__main__':
    t0 = time.time(); H = Hs(0.3)
    print(f'=== V339/rl2 : N={N} jauge={GA} ω={OM} m_tot={MT} μ={MU} γ_F={SG} conv={CONV} t_f={TF} pas={NTH} ; dim {H.shape[0]} ; herm {abs(H - H.getH()).max():.1e}', flush=True)
    ths = np.linspace(0, 2 * np.pi, NTH, endpoint=False); dth = ths[1] - ths[0]; h = 1e-6
    WIN = {'box': (EPS, 'box'), 'box/2': (EPS / 2, 'box'), 'box*2': (2 * EPS, 'box'), 'gauss': (EPS / 2, 'gauss')}
    flow = {k: np.zeros(N + 1) for k in WIN}; gross = np.zeros(N + 1); kmax_used = 0; nwin_max = 0
    KE = 40 if N <= 2 else 64
    for th in ths:
        H = Hs(th); D = (Hs(th + h) - Hs(th - h)) / (2 * h)
        if H.shape[0] <= 700: E, V = np.linalg.eigh(H.toarray())
        else:
            lu = sla.splu(H.tocsc()); Op = sla.LinearOperator(H.shape, matvec=lu.solve, dtype=complex)
            E, V = sla.eigsh(H, k=KE, sigma=0.0, OPinv=Op, which='LM')
            assert np.max(np.abs(E)) > 2.2 * EPS, f'fenêtre saturée à θ={th:.3f} : max|E| = {np.max(np.abs(E)):.3f}'
        nwin_max = max(nwin_max, int(np.sum(np.abs(E) < EPS)))
        for i in np.where(np.abs(E) < 2.2 * EPS)[0]:
            v = V[:, i]; sl = np.real(np.vdot(v, D @ v)); w = prof(v)
            for k, (e, typ) in WIN.items():
                if typ == 'box': wt = (abs(E[i]) < e) / (2 * e)
                else: wt = np.exp(-E[i] ** 2 / (2 * e ** 2)) / np.sqrt(2 * np.pi) / e
                flow[k] += dth * sl * w * wt
            if abs(E[i]) < EPS: gross += dth / (2 * EPS) * abs(sl) * w
    for k in WIN: print(f'  [{k:6s}] flot net par couche : [{" ".join(f"{x:+.3f}" for x in flow[k])}] ; total {flow[k].sum():+.3f} ; couche fine {flow[k][N]:+.4f}')
    print(f'  brut (box) : [{" ".join(f"{x:.3f}" for x in gross)}] ; max états dans la fenêtre {nwin_max} ; temps {time.time() - t0:.0f} s', flush=True)
