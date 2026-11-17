# Affinement 4 (critères : CRITERES_AFFINEMENT_4_TYPEI.md)
import numpy as np, math
rng = np.random.default_rng(11)
Bz = np.array([[1, -1, 0, 0], [1, 1, -2, 0], [1, 1, 1, -3]], float); Bz = (Bz.T / np.linalg.norm(Bz, axis=1)).T
def rho(g):
    P = np.zeros((4, 4))
    for i in range(4): P[g[i], i] = 1
    return Bz @ P @ Bz.T
C = rho((0, 2, 3, 1)); Z = rho((1, 0, 2, 3))
w, V = np.linalg.eig(C)
v1 = np.real(V[:, np.argmin(np.abs(w-1))]); v1 /= np.linalg.norm(v1)
vo = V[:, np.argmin(np.abs(w-np.exp(2j*np.pi/3)))]; vo /= np.linalg.norm(vo)
Q = np.column_stack([v1, vo, vo.conj()])
Pax = np.outer(v1, v1); Jr = (C - C.T)/np.sqrt(3)
def z2tex(spec):
    wz, Uz = np.linalg.eigh(Z); u = Uz[:, np.argmin(wz)] if np.sum(wz < 0) == 1 else Uz[:, np.argmax(wz)]
    Qz = np.linalg.qr(np.column_stack([u, rng.normal(size=(3, 2))]))[0]
    th = rng.uniform(0, np.pi); c, s = np.cos(th), np.sin(th); Qz = Qz @ np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return Qz @ np.diag(np.array(spec)[rng.permutation(3)]) @ Qz.T
def angles(U):
    s13 = abs(U[0,2])**2; s12 = abs(U[0,1])**2/(1-s13); s23 = abs(U[1,2])**2/(1-s13)
    J = np.imag(U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0]))
    c12, c23, c13 = math.sqrt(1-s12), math.sqrt(1-s23), math.sqrt(1-s13)
    sd = J/(c12*math.sqrt(s12)*c23*math.sqrt(s23)*c13**2*math.sqrt(s13))
    S12, S23, S13 = math.sqrt(s12), math.sqrt(s23), math.sqrt(s13)
    cd = (abs(U[1,0])**2 - (S12*c23)**2 - (c12*S23*S13)**2)/(2*S12*c23*c12*S23*S13)
    return s12, s23, math.degrees(math.atan2(max(-1, min(1, sd)), max(-1, min(1, cd)))), s13
def lepton_angles(m):
    # m réelle symétrique dans la base géométrique ; base des leptons chargés : Q (e = v1)
    mf = Q.conj().T @ m @ Q.conj()          # m dans la base de saveur (κ_f = Q† m Q*)
    wv, Vv = np.linalg.eigh(mf.conj().T @ mf); o = np.argsort(wv)
    return angles(Vv[:, o])
worst23 = worstd = 0.0; ds12 = []
for eps in (0.0, 0.05, 0.2):
    w23 = wd = 0.0
    for spec in ((1, 1, 1), (0.5, 1, 1), (0.1, 0.5, 1)):
        for R in (0.3, 1, 3, 10, 30, 100):
            for _ in range(300):
                al, be, ga = rng.uniform(0.5, 1.5, 2).tolist() + [rng.uniform(-0.5, 0.5)]
                Y = al*Pax + be*(np.eye(3)-Pax) + ga*Jr
                MR = z2tex(spec); YD = z2tex(rng.permutation([0.02, 0.3, 1.0]))
                mI = Y.T @ np.linalg.inv(MR) @ Y; mII = YD/np.linalg.norm(YD)*np.linalg.norm(mI)*R
                m = mI + mII
                if eps:
                    S = Y.T @ Y; S /= np.linalg.norm(S); A = np.eye(3) + eps*S; m = A @ m @ A
                s12, s23, dl, s13 = lepton_angles(m)
                w23 = max(w23, abs(s23-0.5)); wd = max(wd, abs(abs(dl)-90))
                if eps == 0: ds12.append((R, s12, s13))
    print(f'ε = {eps}: max |s23 − 0.5| = {w23:.2e} ; max ||δ| − 90°| = {wd:.2e}')
    worst23 = max(worst23, w23); worstd = max(worstd, wd)
print('critère :', 'PASSE' if worst23 < 1e-8 and worstd < 1e-6 else 'ÉCHOUE')
ds12 = np.array(ds12)
for R in (0.3, 1, 3, 10, 30, 100):
    x = ds12[ds12[:, 0] == R]; tm1 = 1 - 2/(3*(1-x[:, 2]))
    print(f'  R = {R:5.1f} : médiane |Δ sin²θ12| (contre TM1 au même s13) = {np.median(np.abs(x[:,1]-tm1)):.4f}')
