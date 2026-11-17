"""V357 (relecture de V352 §3.7) : rapports de traces d'une génération et couplages MS à une boucle (normalisation GUT), recalcul indépendant."""
import numpy as np
from fractions import Fraction as Fr
# (multiplicité couleur × isospin, Y) pour les 16 Weyl d'une génération (Q = T3 + Y)
reps = [('Q_L', 3, 2, Fr(1, 6)), ('u_R', 3, 1, Fr(2, 3)), ('d_R', 3, 1, Fr(-1, 3)), ('L', 1, 2, Fr(-1, 2)), ('e_R', 1, 1, Fr(-1)), ('nu_R', 1, 1, Fr(0))]
TrY2 = sum(c * w * Y ** 2 for _, c, w, Y in reps); TrT2_2 = sum(Fr(1, 2) * c for _, c, w, Y in reps if w == 2); TrT2_3 = sum(Fr(1, 2) * w for _, c, w, Y in reps if c == 3)
print('Tr Y² =', TrY2, ' Tr T²(SU2) =', TrT2_2, ' Tr T²(SU3) =', TrT2_3, ' rapport Tr Y²/Tr T² =', TrY2 / TrT2_2)
aem, s2w, as_, MZ = 1 / 127.95, 0.2312, 0.118, 91.1876
a1i = 0.6 * (1 - s2w) / aem; a2i = s2w / aem; a3i = 1 / as_; b = np.array([41 / 10, -19 / 6, -7])
for L in (1e16, 1e17, 2.4e18, 1.22e19):
    a = np.array([a1i, a2i, a3i]) - b / (2 * np.pi) * np.log(L / MZ)
    print(f'Λ = {L:.2e} GeV : α⁻¹ = {np.round(a, 1).tolist()} ; écart α1/α2 = {100 * (a[1] / a[0] - 1):+.0f} % (α2⁻¹/α1⁻¹ − 1) ; {100 * (a[0] / a[1] - 1):+.0f} % (α1⁻¹/α2⁻¹ − 1) ; α3/α2 : {100 * (a[1] / a[2] - 1):+.0f} %')
