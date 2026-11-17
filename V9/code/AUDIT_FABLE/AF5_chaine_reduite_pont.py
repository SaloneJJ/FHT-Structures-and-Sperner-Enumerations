#!/usr/bin/env python3
"""AUDIT FABLE -- A3, sections 8-9 (registre C14, M1, M2, M4) : chaine reduite ecrite independamment.
Modele : octet au cone (q = 0) x une direction de profondeur ouverte (n = 0 : couche fine ... n = L-1 : racine), sauts verticaux de Kaplan, r = 1,
secteurs gamma_F = +-1 (fibre jouet gamma_F = sigma_z, D_F = sigma_x) :
   H = K (x) Gamma4 (x) 1 + (W + m_K) (x) Gamma5 (x) sigma_z + mu 1 (x) F Gamma5 (x) 1 + eps_D P_h (x) O (x) sigma_x ,
   K = (V - V^+)/(2i), W = (r/2)(2 - V - V^+), V = decalage d'un cran, m_K = sqrt3 (eps - 1) (dimerisation eps = t_f / t_c), mu = 1.
On verifie :
  (a) Proposition « Two branches » : nombre et gout des modes de bord, decroissance lambda = sqrt3 (1 - eps) (survivant) et |2 - sqrt3 (1 - eps)| (jumeau),
      transition a eps = 1 - 1/sqrt3 ;
  (b) Proposition « No mass through Gamma5 x D_F » : O = Gamma5 laisse le survivant sans masse ;
  (c) Proposition « Insertion O_phi » (iii) : m = eps_D (1 - lambda^2) ;
  (d) Proposition « Bridge » : m(h) = eps_D (1 - lambda^2) lambda^(2h).
Usage : python3 AF5_chaine_reduite_pont.py       (numpy)
"""
import numpy as np
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex); I2 = np.eye(2)
k3 = lambda a, b, c: np.kron(np.kron(a, b), c)
X = [k3(sx, I2, I2), k3(sz, sx, I2), k3(sz, sz, sx)]; Y = [k3(sy, I2, I2), k3(sz, sy, I2), k3(sz, sz, sy)]
G5 = sum(X) / np.sqrt(3); A = (X[0] - X[1]) / np.sqrt(2); B = (X[0] + X[1] - 2 * X[2]) / np.sqrt(6)
F = 1j * A @ B; G4 = G5 @ Y[0] @ Y[1] @ Y[2]; FG = F @ G5; O0 = A @ Y[0] @ Y[1] @ Y[2]

def chain(L, eps, mu=1.0, r=1.0, epsD=0.0, O=None, h=None, orient=+1):
    V = np.diag(np.ones(L - 1), 1).astype(complex)                 # V|n+1> = |n>
    K = orient * (V - V.conj().T) / 2j; W = (r / 2) * (2 * np.eye(L) - V - V.conj().T); mK = np.sqrt(3) * (eps - 1)
    H = np.kron(np.kron(K, G4), I2) + np.kron(np.kron(W + mK * np.eye(L), G5), sz) + mu * np.kron(np.kron(np.eye(L), FG), I2)
    if epsD:
        P = np.eye(L) if h is None else np.diag([1.0 if n == h else 0.0 for n in range(L)])
        H = H + epsD * np.kron(np.kron(P, O), sx)
    return H

def boundary_modes(L, eps, orient=+1, tol=1e-6):
    """modes de |E| < tol ; renvoie pour chaque extremite (fine n < L/2, racine) : nombre, <F>, <Gamma_chir>, decroissance"""
    H = chain(L, eps, orient=orient); w, v = np.linalg.eigh(H)
    sel = np.abs(w) < tol; out = {}
    if not sel.any(): return out, H
    Vz = v[:, sel]
    Pn = np.kron(np.kron(np.diag([1.0 if n < L // 2 else 0.0 for n in range(L)]), np.eye(8)), I2)
    # separer les deux extremites : diagonaliser le projecteur « moitie fine » dans l'espace des modes nuls
    pw, pu = np.linalg.eigh(Vz.conj().T @ Pn @ Vz); Vz = Vz @ pu
    for name, mask in (('fine', pw > 0.5), ('racine', pw < 0.5)):
        if not mask.any(): continue
        Vm = Vz[:, mask]
        Fop = np.kron(np.kron(np.eye(L), F), I2); gF = np.kron(np.kron(np.eye(L), np.eye(8)), sz)
        prof = (np.abs(Vm) ** 2).sum(1).reshape(L, 16).sum(1)
        n0 = 0 if name == 'fine' else L - 1; n1 = 1 if name == 'fine' else L - 2; n2 = 2 if name == 'fine' else L - 3
        lam = np.sqrt(prof[n2] / prof[n1]) if prof[n1] > 1e-200 and prof[n2] > 0 else 0.0
        if prof[n1] < 1e-25: lam = 0.0
        fl = np.real(np.trace(Vm.conj().T @ Fop @ gF @ Vm)) / Vm.shape[1]       # <F gamma_F> : -1 pour le survivant (f = -gamma_F), +1 pour le jumeau
        out[name] = (Vm.shape[1], round(fl, 3), lam)
    return out, H

if __name__ == '__main__':
    ok = True
    print('(a) Deux branches (chaine de 160 couches, mu = r = 1) ; attendu : survivant pour eps > 1 - 1/sqrt3 = 0.4226, jumeau pour 1 - sqrt3 < eps < 0.4226')
    print('    eps   | modes nuls cote fin | <F gamma_F> | lambda mesure | lambda predit')
    for eps in (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1, 0.0, -0.5):
        out, _ = boundary_modes(160, eps, tol=1e-7)
        Ms, Mt = np.sqrt(3) * (eps - 1) - 1, np.sqrt(3) * (eps - 1) + 1
        surv, twin = -2 < Ms < 0, -2 < Mt < 0
        pred = abs(1 + Ms) if surv else (abs(1 + Mt) if twin else None)
        got = out.get('fine') or out.get('racine')
        sides = {k: v for k, v in out.items()}
        exp_f = -1.0 if surv else (1.0 if twin else None)
        line_ok = (got is not None) == (pred is not None) and (got is None or (abs(got[1] - exp_f) < 1e-3 and abs(got[2] - pred) < 2e-3))
        ok &= line_ok
        print(f'    {eps:5.2f} | {sides} | attendu <F gamma_F> = {exp_f}, lambda = {None if pred is None else round(pred, 4)}  {"ok" if line_ok else "ECART"}')
    print('    (chaque extremite porte 2 modes nuls : un par secteur gamma_F = +-1 ; l\'extremite opposee porte la contrepartie)')

    print('(b) Insertion Gamma5 x D_F uniforme, eps = 1 (lambda = 0), eps_D = 0.3 : plus petites |E|')
    w = np.sort(np.abs(np.linalg.eigvalsh(chain(24, 1.0, epsD=0.3, O=G5))))[:8]
    print('    ', np.array2string(w, precision=2), '-> les 8 modes de bord (2 extremites x 2 secteurs x 2 spins) restent a E = 0 (a l\'hybridation de taille finie pres, < 1e-6) :', 'ok' if w[7] < 1e-6 else 'ECART'); ok &= w[7] < 1e-6
    Hm = chain(24, 1.0, mu=0.0, epsD=0.3, O=G5) + np.kron(np.kron(np.eye(24), -1.0 * G5), sz)    # controle : mu = 0, m_tot = -1 (deux copies de la piste B)
    w = np.sort(np.abs(np.linalg.eigvalsh(Hm)))[:8]
    print('     controle mu = 0, m_tot = -1 :', np.round(w, 4), '-> masse eps_D = 0.3 (les deux copies ont alors le meme gout dans les deux secteurs)'); ok &= abs(w[0] - 0.3) < 1e-9
    print('(c) Insertion O_0 x D_F : sur TOUTES les couches (modele de la Proposition « Insertion O_phi ») et sur la seule couche fine (h = 0)')
    print('     au cone, O_0 x sigma_x anticommute avec tous les termes de H_0 : H^2 = H_0^2 + eps_D^2 exactement pour l\'insertion uniforme')
    H0 = chain(12, 0.8); Ins = np.kron(np.kron(np.eye(12), O0), sx)
    ac = np.abs(H0 @ Ins + Ins @ H0).max(); print(f'     |{{H_0, O_0 x sigma_x}}| = {ac:.1e}'); ok &= ac < 1e-12
    for eps in (1.0, 0.9, 0.8, 0.6, 0.2):
        Ms = np.sqrt(3) * (eps - 1) - 1; Mt = np.sqrt(3) * (eps - 1) + 1
        lam = abs(1 + Ms) if -2 < Ms < 0 else abs(1 + Mt)
        mu_ = np.sort(np.abs(np.linalg.eigvalsh(chain(40, eps, epsD=1e-3, O=O0))))[0]
        mf = np.sort(np.abs(np.linalg.eigvalsh(chain(40, eps, epsD=1e-3, O=O0, h=0))))[4]
        print(f'     eps = {eps:.1f}, lambda = {lam:.4f} : uniforme m/eps_D = {mu_ / 1e-3:.6f} (attendu 1) ; couche fine seule m/eps_D = {mf / 1e-3:.6f} (attendu 1 - lambda^2 = {1 - lam**2:.6f})')
        ok &= abs(mu_ / 1e-3 - 1) < 1e-6 and abs(mf / 1e-3 - (1 - lam**2)) < 2e-3 * (1 - lam**2) + 1e-9
    print('     => « m = eps_D (1 - lambda^2) » (new_sec9, Prop. Insertion (iii)) vaut pour l\'insertion sur la couche fine seule ; pour l\'insertion uniforme m = eps_D.')
    print('(d) Pont : insertion sur la seule couche h ; m(h) / [eps_D (1 - lambda^2) lambda^(2h)]')
    for eps in (0.9, 0.8, 0.7, 0.6, 0.2, 0.1):
        Ms = np.sqrt(3) * (eps - 1) - 1; Mt = np.sqrt(3) * (eps - 1) + 1
        lam = abs(1 + Ms) if -2 < Ms < 0 else abs(1 + Mt)
        out, _ = boundary_modes(40, eps); side = 'fine' if 'fine' in out and out['fine'][0] else 'racine'
        ratios = []
        for h in range(6):
            w = np.sort(np.abs(np.linalg.eigvalsh(chain(40, eps, epsD=1e-3, O=O0, h=h))))
            # les 4 modes de l'extremite racine restent a E ~ 0 : la masse des modes du cote fin est la 5e plus petite |E|
            m = w[4]; ratios.append(m / (1e-3 * (1 - lam**2) * lam**(2 * h)))
        good = all(abs(x - 1) < 5e-3 for x in ratios); ok &= good
        print(f'     eps = {eps:.1f} (lambda = {lam:.4f}) : h = 0..5 -> {np.round(ratios, 4)}  {"ok" if good else "ECART"}')
    print('VERDICT :', 'conforme' if ok else 'ECART')
