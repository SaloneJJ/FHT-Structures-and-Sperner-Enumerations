"""Integrale de la fenetre carree sans erreur de discontinuite : entre deux points du balayage, chaque etat est apparie
(prediction E + pente x h), E(theta) est interpole lineairement, et on integre exactement l indicatrice |E| < eps
fois (pente x poids) interpole. Usage : python3 rv3_an_exact.py q tags(sep ,)"""
import sys, numpy as np
q = float(sys.argv[1]); P = 2 * np.pi / q; fac = 2 * int(round(2 * np.pi / P))
TH, EE, SS, WW = [], [], [], []
for t in sys.argv[2].split(','):
    z = np.load(f'/home/claude/fht/relecture3/scan_{t}.npz'); ok = ~np.isnan(z['E'][:, 0]); eps0 = float(z['eps'])
    TH.append(z['ths'][ok]); EE.append(z['E'][ok]); SS.append(z['S'][ok]); WW.append(z['W'][ok])
TH = np.concatenate(TH); EE = np.concatenate(EE); SS = np.concatenate(SS); WW = np.concatenate(WW)
o = np.argsort(TH, kind='stable'); TH, EE, SS, WW = TH[o], EE[o], SS[o], WW[o]
_, iu = np.unique(np.round(TH, 10), return_index=True); TH, EE, SS, WW = TH[iu], EE[iu], SS[iu], WW[iu]
def frac_inside(Ea, Eb, eps):
    """fraction de [0,1] ou |Ea + (Eb - Ea) t| < eps"""
    if Ea == Eb: return float(abs(Ea) < eps)
    t1, t2 = sorted([(-eps - Ea) / (Eb - Ea), (eps - Ea) / (Eb - Ea)])
    return max(0.0, min(1.0, t2) - max(0.0, t1))
for eps in (eps0, eps0 / 2):
    F = np.zeros(WW.shape[1]); B = np.zeros(WW.shape[1]); nbad = 0
    for i in range(len(TH) - 1):
        h = TH[i + 1] - TH[i]; used = set()
        for k in np.where(np.abs(EE[i]) < 2 * eps)[0]:
            pred = EE[i, k] + SS[i, k] * h; d = np.abs(EE[i + 1] - pred); d[list(used)] = np.inf; k2 = int(np.argmin(d)); used.add(k2)
            if d[k2] > 0.2 * eps: nbad += 1
            fr = frac_inside(EE[i, k], EE[i + 1, k2], eps)
            sw = 0.5 * (SS[i, k] * WW[i, :, k] + SS[i + 1, k2] * WW[i + 1, :, k2]); swa = 0.5 * (abs(SS[i, k]) * WW[i, :, k] + abs(SS[i + 1, k2]) * WW[i + 1, :, k2])
            F += fr * h * sw / (2 * eps); B += fr * h * swa / (2 * eps)
        # etats entrant dans la fenetre depuis l exterieur de 2 eps a gauche : negligeable si les pas sont petits
    F *= fac; B *= fac
    print(f'q = {q}, fenetre {eps:.4f} : flot [' + ' '.join(f'{x:+.3f}' for x in F) + f'] total {F.sum():+.3f} ; brut fine {B[-1]:.3f} ; appariements douteux {nbad}')
