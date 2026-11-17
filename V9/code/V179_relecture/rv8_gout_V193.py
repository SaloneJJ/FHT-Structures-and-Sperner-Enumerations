# Relecture V193 : le terme γ⊥ brise-t-il le goût ? Commutant de {H(k)} pour tous k (t_f ≠ t_c).
import numpy as np, itertools as it
exec(open("rv4_dimer_V193.py").read().split("# (a) H(K)")[0])   # fonction bloch du relecteur (rv4), pas du code relu
rng = np.random.default_rng(8)
Hs = [bloch(rng.uniform(-np.pi, np.pi, 3), 0.6, 1.0) for _ in range(12)]
M = np.vstack([np.kron(h, np.eye(8)) - np.kron(np.eye(8), h.T) for h in Hs])
s = np.linalg.svd(M, compute_uv=False)
print("dimension du commutant de {H(k)} (12 k aléatoires, t_f=0.6, t_c=1) :", int(np.sum(s < 1e-9)))
# Γ4 unique (V183/V193) : comparaison avec i Γ5 α1 α2 α3, et commutation avec γ⊥
exec("G5" + open("rv4_dimer_V193.py").read().split("G5", 1)[1].split("worst = 0")[0])
G4c = G5 @ alpha[0] @ alpha[1] @ alpha[2]
print("Γ4c = Γ5α1α2α3 hermitien : %s, carré 1 : %s, anticommute avec Γ5 et α_d : %s"
      % (np.allclose(G4c, G4c.conj().T), np.allclose(G4c @ G4c, np.eye(8)),
         all(np.allclose(G4c @ Y + Y @ G4c, 0) for Y in [G5] + alpha)))
from rv_common import ROTS, lift_diag
Us = [lift_diag(Hcube, pos, R)[0] for R in ROTS]
print("Γ4c invariant (A1) :", all(np.allclose(U @ G4c @ U.conj().T, G4c) for U in Us))
print("[Γ4c, γ⊥_d] = 0 :", [bool(np.allclose(G4c @ g - g @ G4c, 0)) for g in gperp], "; {Γ4c, γ⊥_d} = 0 :", [bool(np.allclose(G4c @ g + g @ G4c, 0)) for g in gperp])
