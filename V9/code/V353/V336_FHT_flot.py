"""V335 : C3' -- piste B native (octet BBH, Cl(5) natif, convention « cos ») + terme de goût nu mu*FΓ5, flot spectral sur l'arbre augmenté.
Critères : 01_synthese/V335 (écrits avant calcul). Géométrie, flux, transport : V162 (b, r0) ; flot : Hellmann–Feynman comme V166.
Usage : python3 V335_FHT_C3prime_flot.py N mtot mu s NTH [q]     (s = secteur γ_F = ±1)       JJ SALONE, FHT 2026"""
import sys, time, itertools
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
N_, MT, MU, SG, NTH_ = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
Q_ = sys.argv[6] if len(sys.argv) > 6 else '1'
MODE = sys.argv[7] if len(sys.argv) > 7 else 'bbh'
TFU = float(sys.argv[8]) if len(sys.argv) > 8 else 1.0   # V336 : mode unif = jauge alternée, t_f = TFU, t_c = 1, AUCUN décalage de masse (m_K = √3(t_f − t_c))     # bbh : BBH natif (R(q) anisotrope) ; iso : contrôle, Γ natifs + Wilson isotrope de V151 (deux copies exactes du modèle emprunté)
sys.argv = ['x', N_, '-1', 'b', Q_, 'none', '1', str(NTH_), 'r0', '0/1']
exec(open('V162_FHT_R3_jauge.py').read().split("print(f'=== V162")[0])
from rvA_common import clifford, H as Hbbh, KPT, TF, TC
sxm = np.array([[0, 1], [1, 0]], complex); sym = np.array([[0, -1j], [1j, 0]]); szm = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
Xg = [k3(sxm, I2, I2), k3(szm, sxm, I2), k3(szm, szm, sxm)]; Yg = [k3(sym, I2, I2), k3(szm, sym, I2), k3(szm, szm, sym)]
# Clifford natif exact (identique à V331 §0 au signe près) : Γ5 ∝ ΣX, Am, Bm = composantes E, α ∝ Yg, Γ4 = Γ5 α1 α2 α3
G5 = sum(Xg)/np.sqrt(3); Am = (Xg[0] - Xg[1])/np.sqrt(2); Bm = (Xg[0] + Xg[1] - 2*Xg[2])/np.sqrt(6); G4 = G5@Yg[0]@Yg[1]@Yg[2]; FG = 1j*Am@Bm@G5
Cl = clifford(); print('  accord avec clifford() de V331 (|tr G5 G5ref|/8, |tr FG FGref|/8) :', round(abs(np.trace(G5@Cl['G5']))/8, 6), round(abs(np.trace(FG@(1j*Cl['A']@Cl['B']@Cl['G5'])))/8, 6))
for nm, Z in (('G5', G5), ('G4', G4), ('FG', FG)): assert np.allclose(Z@Z, np.eye(8)) and np.allclose(Z, Z.conj().T), nm
mBBH = np.real(np.trace(Hbbh(KPT)@G5))/8
# contrôles algébriques : BBH réécrit = H(k) de rvA_common ; G4 anticommute avec tout ; FΓ5 commute avec Xg, anticommute avec Yg, G4
kk = np.array([0.3, -1.1, 2.0]); Hk = sum((TF + TC*np.cos(kk[d]))*Xg[d] + TC*np.sin(kk[d])*Yg[d] for d in range(3))
assert np.allclose(Hk, Hbbh(kk)), 'BBH réel != BBH en k'
assert all(np.allclose(G4@Z, -Z@G4) for Z in Yg + [G5]) and all(np.allclose(G4@Z, Z@G4) for Z in (Am, Bm)) and all(np.allclose(FG@Z, Z@FG) for Z in Xg) and all(np.allclose(FG@Z, -Z@FG) for Z in Yg + [G4])
def ops(th):
    U, TIN = link_fields(th, 0.0)
    R, C, Vv = [[], [], []], [[], [], []], [[], [], []]
    for n in range(N + 1):
        s = sides[n]; I, J, K = grid(s); i = lin(n, I, J, K).ravel()
        for jd, sh in enumerate([(1, 0, 0), (0, 1, 0), (0, 0, 1)]):
            nb = lin(n, I + sh[0], J + sh[1], K + sh[2]).ravel()
            R[jd].append(i); C[jd].append(nb); Vv[jd].append(np.conj(U[n][jd]).ravel())
    T = [sp.csr_matrix((np.concatenate(Vv[jd]), (np.concatenate(R[jd]), np.concatenate(C[jd]))), shape=(Nc, Nc)) for jd in range(3)]
    rs, cs, vs = [], [], []
    for n in range(1, N + 1):
        s = sides[n]; I, J, K = grid(s); ch = lin(n, I, J, K).ravel(); pa = lin(n - 1, I // 2, J // 2, K // 2).ravel()
        u = np.zeros((s, s, s), complex)
        for c, P in TIN[n - 1].items(): u[c[0]::2, c[1]::2, c[2]::2] = P
        rs.append(ch); cs.append(pa); vs.append(u.ravel() / np.sqrt(8))
    S = sp.csr_matrix((np.concatenate(vs), (np.concatenate(rs), np.concatenate(cs))), shape=(Nc, Nc))
    return T, S
kr = lambda A_, g: sp.kron(A_, sp.csr_matrix(g), format='csr')
def Hnat(th, mtot=MT, mu=MU, s=SG):
    T, S = ops(th); Id = sp.identity(Nc, format='csr')
    H = sp.csr_matrix((8*Nc, 8*Nc), dtype=complex)
    if MODE == 'iso':
        for d in range(3):
            H = H + kr((T[d] - T[d].getH())/2j, Yg[d]) - s*(r/2)*kr(T[d] + T[d].getH(), G5)
        H = H + s*(mtot + 3*r)*kr(Id, G5)
    else:
      TCs = -TC if MODE in ('bbhG', 'unif') else TC
      TFs = TFU if MODE == 'unif' else TF      # bbhG : jauge alternée ψ -> (-1)^(x+y+z) ψ à chaque niveau (cône en Γ, couplé au parent)
      for d in range(3):
        H = H + TCs*kr((T[d] - T[d].getH())/2j, Yg[d])                                          # impair en q : nu
        H = H + s*(TFs*kr(Id, Xg[d]) + TCs*kr((T[d] + T[d].getH())/2, Xg[d]))                     # pair en q : habillé (« cos »)
      if MODE != 'unif': H = H + s*(mtot - mBBH)*kr(Id, G5)                                                         # masse en K = mtot (habillée)
    H = H + kr((S - S.getH())/2j, G4) + s*(r/2)*kr(2*Id - S - S.getH(), G5)                  # vertical : Γ4 nu, Wilson habillé
    H = H + mu*kr(Id, FG)                                                                      # terme de goût NU
    return H.tocsr()
t0 = time.time(); H0 = Hnat(0.3); print(f'=== V335 [{MODE}] : N = {N}, m_tot = {MT}, μ = {MU}, γ_F = {SG}, q = {q}, {NTH_} pas ; dim {H0.shape[0]} ; herm {abs(H0 - H0.getH()).max():.1e}', flush=True)
EPS = 0.25*2*np.pi/L
def prof(v): return np.bincount(depth, weights=np.sum(np.abs(v.reshape(Nc, 8))**2, axis=1), minlength=N + 1)
ths = np.linspace(0, 2*np.pi, NTH_, endpoint=False); dth = ths[1] - ths[0]; dd = 1e-5
flow = np.zeros(N + 1); gross = np.zeros(N + 1); sat = 0
for th in ths:
    H = Hnat(th); Dm = (Hnat(th + dd) - Hnat(th - dd))/(2*dd)
    if H.shape[0] <= 1500: E, V = np.linalg.eigh(H.toarray())
    else: E, V = sla.eigsh(H, k=48, sigma=1e-4, which='LM')
    sl_ = np.abs(E) < EPS
    if H.shape[0] > 1500 and sl_.sum() >= 46: sat += 1
    for i in np.where(sl_)[0]:
        v = V[:, i]; s_ = np.real(np.vdot(v, Dm@v)); w = prof(v)
        flow += dth/(2*EPS)*s_*w; gross += dth/(2*EPS)*abs(s_)*w
print(f'  fenêtre {EPS:.3f} ; flot net par couche : [{" ".join(f"{x:+.3f}" for x in flow)}] ; total {flow.sum():+.3f} ; saturées {sat}')
print(f'  flot brut par couche : [{" ".join(f"{x:.3f}" for x in gross)}] ; temps {time.time() - t0:.0f} s', flush=True)
