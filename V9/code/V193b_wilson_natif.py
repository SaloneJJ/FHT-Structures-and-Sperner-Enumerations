# V193b — décomposition de H(K+q) : partie Γ5 (masse + Wilson natif) et partie hors Γ5 (brisure de goût) à l'ordre q².
import numpy as np
exec(open("V193_FHT_dimerisation_et_coins.py").read().split('print("=== Q1')[0])
exec("def Hk" + open("V193_FHT_dimerisation_et_coins.py").read().split("def Hk")[1].split("K = np.array")[0])
K = np.array([np.pi]*3); Hc = cube(True); G5 = Hc/np.sqrt(3)
for tf_, tc_ in [(1.0, 0.8), (0.6, 1.0)]:
    for qv in [np.array([0.05, 0, 0]), np.array([0.05, 0.05, 0.05])]:
        H = Hk(K + qv, tf_, tc_); Hs = (H + G5 @ H @ G5)/2          # partie qui commute avec Γ5
        c5 = np.real(np.trace(G5 @ H))/8                              # coefficient de Γ5
        rest = Hs - c5*G5
        m = np.sqrt(3)*(tf_ - tc_); rW = tc_/np.sqrt(3)
        pred = m + rW*np.sum(1 - np.cos(qv))
        print(" t_f=%.1f t_c=%.1f q=%s : coeff Γ5 = %.6f ; prédit m + r_W Σ(1 − cos q) = %.6f (m = √3(t_f − t_c), r_W = t_c/√3) ; reste commutant (hors Γ5, brisure de goût) ‖·‖/8 = %.2e"
              % (tf_, tc_, qv.tolist(), c5, pred, np.linalg.norm(rest)/8))
