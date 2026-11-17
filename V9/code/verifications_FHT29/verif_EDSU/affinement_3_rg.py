# Affinement 3 : course RG à une boucle (critères : CRITERES_AFFINEMENT_3_RG.md)
import numpy as np, math
s13 = 0.02248; dm21 = 7.530e-5; dm31 = 2.511e-3
ytau, ymu = 0.0102, 6.1e-4; L = math.log(91.19/1e16); C0 = -1.5
def U_tm1(s13, s23=0.5, dl=-90.0):
    s12 = 1 - 2/(3*(1-s13)); t12, t13, t23 = (math.asin(math.sqrt(x)) for x in (s12, s13, s23))
    c12, c13, c23, S12, S13, S23 = math.cos(t12), math.cos(t13), math.cos(t23), math.sin(t12), math.sin(t13), math.sin(t23)
    e = np.exp(1j*math.radians(dl))
    return np.array([[c12*c13, S12*c13, S13/e],
                     [-S12*c23 - c12*S23*S13*e, c12*c23 - S12*S23*S13*e, S23*c13],
                     [S12*S23 - c12*c23*S13*e, -c12*S23 - S12*c23*S13*e, c23*c13]])
def angles(U):
    s13 = abs(U[0,2])**2; s12 = abs(U[0,1])**2/(1-s13); s23 = abs(U[1,2])**2/(1-s13)
    J = np.imag(U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0]))
    c12, c23, c13 = math.sqrt(1-s12), math.sqrt(1-s23), math.sqrt(1-s13)
    sd = J/(c12*math.sqrt(s12)*c23*math.sqrt(s23)*c13**2*math.sqrt(s13))
    S12, S23, S13 = math.sqrt(s12), math.sqrt(s23), math.sqrt(s13)
    cd = (abs(U[1,0])**2 - (S12*c23)**2 - (c12*S23*S13)**2)/(2*S12*c23*c12*S23*S13)
    return s12, s23, math.degrees(math.atan2(max(-1, min(1, sd)), max(-1, min(1, cd)))), s13
def run(m1, eta, fac, U0):
    m = np.array([m1, math.sqrt(m1**2+dm21), math.sqrt(m1**2+dm31)])
    k = U0.conj() @ np.diag(eta*m) @ U0.conj().T
    eps = np.array([0.0, fac*C0*ymu**2*L/(16*math.pi**2), fac*C0*ytau**2*L/(16*math.pi**2)])
    D = np.diag(1+eps); kl = D @ k @ D
    w, V = np.linalg.eigh(kl.conj().T @ kl); o = np.argsort(w); return angles(V[:, o]), eps
U0 = U_tm1(s13); a0 = angles(U0)
print(f'à Λ : s12={a0[0]:.5f} s23={a0[1]:.5f} δ={a0[2]:+.3f}° s13={a0[3]:.5f}')
for fac in (1, 2, -2, 10, -10):
    worst = [0, 0, 0]
    for m1 in np.linspace(0, 4.25e-3, 41):
        for e2 in (1, -1):
            for e3 in (1, -1):
                (s12, s23, dl, _), eps = run(m1, np.array([1, e2, e3]), fac, U0)
                d = [abs(s12-a0[0]), abs(s23-a0[1]), abs(((dl-a0[2]+180) % 360)-180)]
                worst = [max(x, y) for x, y in zip(worst, d)]
    print(f'facteur {fac:+d} (ε_τ = {eps[2]:+.2e}) : max |Δs12| = {worst[0]:.2e}, |Δs23| = {worst[1]:.2e}, |Δδ| = {worst[2]:.3f}°')
print('--- grandeur complémentaire : m1 (quasi dégénéré) pour |Δs23| = 0.035, facteur 1')
for m1 in [0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
    best = max(abs(run(m1, np.array([1, e2, e3]), 1, U0)[0][1]-0.5) for e2 in (1, -1) for e3 in (1, -1))
    print(f'  m1 = {m1:5.2f} eV : max |Δs23| = {best:.2e} ; Σmν ≈ {3*m1:.2f} eV')
