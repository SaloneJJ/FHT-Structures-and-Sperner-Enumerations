# Controle V146 §3 / V144 §5.3 : "la torsion par gamma_F est forcee". Alternative Gamma' = gamma_F . (+1 sur I, -1 sur II).
import io, contextlib, itertools, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    import V146_FHT_produit_bord as P
dF, dI, GF, UF, DF, DS, I2, iSy = P.dF, P.dI, P.GF, P.UF, P.DF, P.DS, P.I2, P.iSy
S3 = np.diag([1.0] * dI + [-1.0] * dI)
Gp = GF @ S3
def test(DFm, nom):
    k = np.array([0.4, -0.3, 0.8])
    D1 = np.kron(DFm, I2) + np.kron(Gp, DS(k)); D1m = np.kron(DFm, I2) + np.kron(Gp, DS(-k))
    pyth = np.abs(D1 @ D1 - np.kron(DFm @ DFm, I2) - (k @ k) * np.eye(2 * dF)).max()
    out = [f"{nom} : |D'^2 - D_F^2 - k^2| = {pyth:.1e}"]
    for Jnom, U in (("J_F g_F (x) J_S", np.kron(UF @ GF, iSy).real), ("J_F (x) J_S", np.kron(UF, iSy).real)):
        J2 = np.round((U @ U.conj())[0, 0].real, 6)
        a = np.abs(U @ D1.conj() @ U.T + D1m).max(); b = np.abs(U @ D1.conj() @ U.T - D1m).max()
        out.append(f"{Jnom} : J^2 = {J2:+.0f}, JD+DJ {a:.1e}, JD-DJ {b:.1e}")
    # helicite des etats d'energie positive du secteur II, doublet d1, pour D_F = 0
    H = np.kron(Gp, DS(k)); w, Vv = np.linalg.eigh(H); pos = Vv[:, w > 1e-9]
    proj = np.zeros(dF); proj[dI + P.iI(1, 0, 0)] = 1
    Pm = np.kron(np.diag(proj), I2); hel = np.kron(np.eye(dF), DS(k / np.linalg.norm(k)))
    h = np.real(np.trace(pos.conj().T @ Pm @ hel @ Pm @ pos) / np.trace(pos.conj().T @ Pm @ pos))
    out.append(f"helicite (D_F=0) du doublet d1 de II : {h:+.0f}")
    print(" ; ".join(out))
DF0 = DF.copy(); DF0[dI:, :dI] = 0; DF0[:dI, dI:] = 0      # sans terme de Majorana
test(DF0, "sans Majorana"); test(DF, "avec Majorana")
