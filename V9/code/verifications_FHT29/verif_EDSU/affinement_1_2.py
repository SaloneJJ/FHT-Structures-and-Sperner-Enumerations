# Affinements 1 et 2 (critères : CRITERES_AFFINEMENT_1_2.md, écrits avant ce script)
import sys, os, math, itertools as it
import numpy as np
class chi2:
    @staticmethod
    def sf(x, k):
        if x <= 0: return 1.0
        if k == 1: return math.erfc(math.sqrt(x/2))
        if k == 2: return math.exp(-x/2)
        if k == 3: return math.erfc(math.sqrt(x/2)) + math.sqrt(2*x/math.pi)*math.exp(-x/2)
        raise ValueError
class norm:
    @staticmethod
    def isf(p):
        lo, hi = 0.0, 40.0
        for _ in range(200):
            m = (lo+hi)/2
            if 0.5*math.erfc(m/math.sqrt(2)) > p: lo = m
            else: hi = m
        return lo
D = sys.argv[1]
FILES = {('yes','NO'):'506694.txt', ('yes','IO'):'460581.txt', ('off','NO'):'144289.txt', ('off','IO'):'196693.txt'}
WANT = ['T23/DCP', 'T12/DCP', 'T12', 'T23']
def read(path):
    secs, cur = {}, None
    for line in open(path):
        if line.startswith('###'): continue
        if line.startswith('#'):
            tag = line[1:].split('projection')[0].strip(); cur = tag if tag in WANT else None
            if cur: secs[cur] = []
            continue
        if cur and line.strip(): secs[cur].append([float(x) for x in line.split()])
    return {k: np.array(v) for k, v in secs.items()}
def grid2(tab):
    xs = np.unique(tab[:, 0]); ys = np.unique(tab[:, 1]); Z = np.full((len(xs), len(ys)), np.nan)
    ix = {v: k for k, v in enumerate(xs)}; iy = {v: k for k, v in enumerate(ys)}
    for x, y, z in tab[:, :3]: Z[ix[x], iy[y]] = z
    return xs, ys, Z
def i2(g, x, y):
    xs, ys, Z = g
    y = ((y + 180) % 360) - 180
    i = min(max(np.searchsorted(xs, x) - 1, 0), len(xs)-2); j = min(max(np.searchsorted(ys, y) - 1, 0), len(ys)-2)
    u = (x-xs[i])/(xs[i+1]-xs[i]); v = (y-ys[j])/(ys[j+1]-ys[j])
    return float((1-u)*(1-v)*Z[i,j] + u*(1-v)*Z[i+1,j] + (1-u)*v*Z[i,j+1] + u*v*Z[i+1,j+1])
def i1(t, x): return float(np.interp(x, t[:, 0], t[:, 1]))
def pstr(c, dof):
    p = chi2.sf(c, dof); z = norm.isf(p/2)
    v = 'compatible' if p > 0.05 else ('en tension' if p > 0.0027 else 'exclu')
    return f'χ²={c:.2f} ({dof} ddl) p={p:.3g} ≈{z:.1f}σ {v}'
# ---------- courbe GJ (V358_B_bande, recopiée) ----------
S4 = list(it.permutations(range(4)))
X = np.array([[1, 0, 0, -1], [0, 1, 0, -1], [0, 0, 1, -1]], float); Qr, _ = np.linalg.qr(X.T); Bs = Qr.T
def rho(g):
    P = np.zeros((4, 4))
    for i in range(4): P[g[i], i] = 1
    return Bs @ P @ Bs.T
c3 = rho((0, 2, 3, 1)); w3, P3 = np.linalg.eig(c3)
v1 = np.real(P3[:, np.argmin(np.abs(w3 - 1))]); v1 /= np.linalg.norm(v1)
vo = P3[:, np.argmin(np.abs(w3 - np.exp(2j * np.pi / 3)))]; vo /= np.linalg.norm(vo)
Q = np.column_stack([v1, vo, vo.conj()])
r01 = rho((1, 0, 2, 3)); w, P = np.linalg.eigh(r01); u = P[:, np.argmin(w)]; p1, p2 = P[:, w > 0].T
y = dict(d=0.65e-5, s=1.31e-4, b=0.607e-2)
F0 = np.sqrt(y['d'] * y['s']); E0 = y['s'] - y['d']; D0 = y['b']
def Ylep(alpha, kappa=1.0, r=-3):
    M = np.array([[0, F0*np.sqrt(kappa)*np.exp(1j*alpha), 0], [F0/np.sqrt(kappa)*np.exp(-1j*alpha), r*E0, 0], [0, 0, D0]])
    return Q @ M @ Q.conj().T
def Ue_of(Yl):
    UL, s, _ = np.linalg.svd(Yl); o = np.argsort(s); return UL[:, o]
def angles(U):
    s13 = abs(U[0, 2])**2; s12 = abs(U[0, 1])**2/(1-s13); s23 = abs(U[1, 2])**2/(1-s13)
    J = np.imag(U[0, 0]*U[1, 1]*np.conj(U[0, 1])*np.conj(U[1, 0]))
    c12, c23, c13 = np.sqrt(1-s12), np.sqrt(1-s23), np.sqrt(1-s13)
    sd = J/(c12*np.sqrt(s12)*c23*np.sqrt(s23)*c13**2*np.sqrt(s13))
    S12, S23, S13 = np.sqrt(s12), np.sqrt(s23), np.sqrt(s13)
    cd = (abs(U[1, 0])**2 - (S12*c23)**2 - (c12*S23*S13)**2)/(2*S12*c23*c12*S23*S13)
    return s13, s12, s23, np.degrees(np.arctan2(np.clip(sd, -1, 1), np.clip(cd, -1, 1)))
def Ufull(Ue, th):
    n2 = np.cos(th)*p1 + np.sin(th)*p2; n3 = -np.sin(th)*p1 + np.cos(th)*p2
    return Ue.conj().T @ np.column_stack([u, n2, n3])
def solve(Ue, s13t):
    ths = np.linspace(0, np.pi, 2001); f = np.array([angles(Ufull(Ue, t))[0] for t in ths]) - s13t; out = []
    for i in np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]:
        a, b = ths[i], ths[i+1]
        for _ in range(50):
            m = (a+b)/2
            if np.sign(angles(Ufull(Ue, m))[0]-s13t) == np.sign(angles(Ufull(Ue, a))[0]-s13t): a = m
            else: b = m
        Uf = Ufull(Ue, (a+b)/2); out.append(Uf)
    return out
def gj_curve(s13t, step=0.5):
    pts = []
    for ad in np.arange(0, 360, step):
        for k, U in enumerate(solve(Ue_of(Ylep(np.radians(ad))), s13t)):
            a = angles(U); pts.append((ad, k, 0, a[1], a[2], a[3]))
            Us = U[[0, 2, 1], :]; b = angles(Us); pts.append((ad, k, 1, b[1], b[2], b[3]))  # orientation échangée
    return np.array(pts)   # α, branche, orientation, s12, s23, δ
print('=== contrôle de non-régression GJ (s13 = 0.0220, orientation 0, pas 5°)')
C = gj_curve(0.0220, step=5.0); C0 = C[C[:, 2] == 0]
print(f'  s12 ∈ [{C0[:,3].min():.4f}, {C0[:,3].max():.4f}] (A2 : [0.269, 0.374]) ; s23 ∈ [{C0[:,4].min():.4f}, {C0[:,4].max():.4f}] (A2 : [0.493, 0.510]) ; |δ| ∈ [{np.abs(C0[:,5]).min():.1f}, {np.abs(C0[:,5]).max():.1f}] (A2 : [68.5, 111.6])')
ok = abs(C0[:,3].min()-0.269) < 0.002 and abs(C0[:,3].max()-0.374) < 0.002 and abs(C0[:,4].min()-0.493) < 0.002 and abs(C0[:,4].max()-0.510) < 0.002 and abs(np.abs(C0[:,5]).min()-68.5) < 0.5 and abs(np.abs(C0[:,5]).max()-111.6) < 0.5
print('  non-régression :', 'OK' if ok else 'ÉCHEC'); 
if not ok: sys.exit(1)
C1 = C[C[:, 2] == 1]
print(f'  orientation échangée : s23 ∈ [{C1[:,4].min():.4f}, {C1[:,4].max():.4f}], δ ∈ [{C1[:,5].min():.1f}, {C1[:,5].max():.1f}]')
# diagnostic CP-conjugué : (s12, s23, −δ) atteint par le balayage (orientation 0) ?
d = [np.min(np.hypot((C0[:,3]-a)*1000, np.hypot((C0[:,4]-b)*1000, C0[:,5]+c))) for a, b, c in C0[:, 3:6]]
print(f'  diagnostic CP-conjugué dans orientation 0 : distance max = {max(d):.2f} (unités : 1e-3 en s², degrés en δ)')
# ---------- évaluation ----------
for var in ('yes', 'off'):
    for o in ('NO', 'IO'):
        s = read(os.path.join(D, FILES[(var, o)]))
        m = s['T23'][:, 1].min(); g = grid2(s['T23/DCP']); g12 = grid2(s['T12/DCP'])
        s13 = 0.02248 if o == 'NO' else 0.02262
        def chiA(s12, s23, dl): return i2(g, s23, dl) + i1(s['T12'], s12) - 2*m
        def chiB(s12, s23, dl): return i2(g12, s12, dl) + i1(s['T23'], s23) - 2*m
        print(f'\n######## {var} {o} (minimum de l\'ordre m = {m:.3f} au-dessus du minimum global)')
        # TM1 exact
        t12 = 1 - 2/(3*(1-s13))
        ex = min((chiA(t12, 0.5, sg*90), sg) for sg in (1, -1)); exB = min(chiB(t12, 0.5, sg*90) for sg in (1, -1))
        print(f'  TM1 exact (s12={t12:.4f}, s23=0.5, δ={ex[1]*90:+d}°) : {pstr(ex[0], 3)} [contrôle B : {exB:.2f}] ; rel. global {ex[0]+m:.2f}')
        # TM1 seul
        best = None; pts = []
        for s23 in np.arange(0.35, 0.65+1e-9, 0.0005):
            c2 = math.cos(2*math.asin(math.sqrt(s23))); s2 = math.sin(2*math.asin(math.sqrt(s23)))
            cd = -(1-5*s13)*(c2/s2)/(2*math.sqrt(2)*math.sqrt(s13)*math.sqrt(1-3*s13))
            if abs(cd) > 1: continue
            for sg in (1, -1):
                dl = sg*math.degrees(math.acos(cd)); c = chiA(t12, s23, dl); pts.append((c, s23, dl))
        pts = np.array(pts); k = pts[:, 0].argmin(); cmin, s23b, dlb = pts[k]
        print(f'  TM1 seul : min en s23={s23b:.4f}, δ={dlb:+.1f}° : {pstr(cmin, 2)} [contrôle B : {chiB(t12, s23b, dlb):.2f}]')
        print(f'    coût de (M3′) = {ex[0]-cmin:.2f} (1 ddl, p={chi2.sf(ex[0]-cmin,1):.3g}) ; rapport χ²(TM1 seul)/χ²(TM1 exact) = {cmin/ex[0]:.2f}')
        for dc in (1, 4):
            sel = pts[pts[:, 0] <= cmin + dc]
            print(f'    prédiction Δχ² ≤ min+{dc} : s23 ∈ [{sel[:,1].min():.3f}, {sel[:,1].max():.3f}], δ ∈ {sorted(set(np.sign(sel[:,2]).astype(int)))} signe(s), |δ| ∈ [{np.abs(sel[:,2]).min():.1f}, {np.abs(sel[:,2]).max():.1f}]°')
        for ds in (-0.00057, 0.00057):
            s13p = s13 + ds; t12p = 1 - 2/(3*(1-s13p)); cc = []
            for s23 in np.arange(0.35, 0.65+1e-9, 0.001):
                c2 = math.cos(2*math.asin(math.sqrt(s23))); s2 = math.sin(2*math.asin(math.sqrt(s23)))
                cd = -(1-5*s13p)*(c2/s2)/(2*math.sqrt(2)*math.sqrt(s13p)*math.sqrt(1-3*s13p))
                if abs(cd) <= 1: cc += [chiA(t12p, s23, sg*math.degrees(math.acos(cd))) for sg in (1, -1)]
            print(f'    sensibilité s13 = {s13p:.5f} : χ²min(TM1 seul) = {min(cc):.2f} ; TM1 exact = {min(chiA(t12p, 0.5, sg*90) for sg in (1,-1)):.2f}')
        # GJ
        G = gj_curve(s13, step=0.5)
        cg = np.array([chiA(a, b, c) for a, b, c in G[:, 3:6]]); kk = cg.argmin()
        al, br, ori, a12, a23, adl = G[kk]
        print(f'  GJ : min en α={al:.1f}°, branche {int(br)}, orientation {int(ori)} : s12={a12:.4f} s23={a23:.4f} δ={adl:+.1f}° : {pstr(cg[kk], 2)} [contrôle B : {chiB(a12, a23, adl):.2f}]')
        for dc in (1, 4):
            sel = G[cg <= cg[kk] + dc]
            print(f'    prédiction Δχ² ≤ min+{dc} : α ∈ [{sel[:,0].min():.0f}, {sel[:,0].max():.0f}]°, s12 ∈ [{sel[:,3].min():.4f}, {sel[:,3].max():.4f}], s23 ∈ [{sel[:,4].min():.4f}, {sel[:,4].max():.4f}], δ ∈ [{sel[:,5].min():.1f}, {sel[:,5].max():.1f}]°')
