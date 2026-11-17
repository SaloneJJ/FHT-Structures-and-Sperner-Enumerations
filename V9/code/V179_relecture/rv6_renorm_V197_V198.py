# Relecture V197–V198 : renormalisation. Tout le modèle est séparable (rv5) : on travaille sur la chaîne 1D, puis on contrôle en 3D.
import numpy as np, itertools as it
from scipy.optimize import least_squares

def depth(x):
    h = 1
    while x % 2 == 1: h += 1; x //= 2
    return h
def chain(N, t):           # t : dict profondeur -> amplitude
    h = np.zeros((N, N))
    for x in range(N - 1): h[x, x + 1] = h[x + 1, x] = t[depth(x)]
    return h
def ks(d, r): return [1, (-1) ** r[0], (-1) ** (r[0] + r[1])][d]
def cluster3d(N, t):
    sites = list(it.product(range(N), repeat=3)); idx = {s: i for i, s in enumerate(sites)}
    H = np.zeros((N ** 3,) * 2)
    for s in sites:
        for d in range(3):
            if s[d] + 1 < N:
                u = list(s); u[d] += 1; H[idx[s], idx[tuple(u)]] = H[idx[tuple(u)], idx[s]] = t[depth(s[d])] * ks(d, s)
    return H

# --- V197 : amplitude uniforme entre octets -> pas de groupe de 64 états séparés ---
print("=== V197 (uniforme) : chaîne 8 sites ε,1,ε,1,ε,1,ε ; 3D 8×8×8 ===")
for eps in [0.2, 0.5]:
    e1 = np.sort(np.abs(np.linalg.eigvalsh(chain(8, {1: eps, 2: 1, 3: 1}))))
    a = np.sort(np.abs(np.linalg.eigvalsh(cluster3d(8, {1: eps, 2: 1, 3: 1}))))
    print(" ε=%.1f : |e| 1D = %s ; 3D : 8e=%.4f 9e=%.4f 64e=%.4f 65e=%.4f" % (eps, np.round(e1[::2], 4).tolist(), a[7], a[8], a[63], a[64]))

# --- V198 : t(1)=ε, t(2)=1, t(3)=τ ---
def fit_eff(target4):
    """chaîne effective à 4 sites (x, y, x) : |e| = √(x²+y²/4) ± y/2 ; solution EXACTE x = √(ab), y = a − b."""
    a, b = max(target4), min(target4); return np.sqrt(a * b), a - b
print("=== V198 : 3D (512 sites) contre 1D (8 sites) ; ajustement 3D par moindres carrés indépendant ===")
def eff3d(p):
    return np.sort(np.linalg.eigvalsh(cluster3d(4, {1: p[0], 2: p[1]})))
print(" ε    τ    | 65e/64e | t_f' (1D exact)  t_c' (1D exact) | t_f', t_c' (ajust. 3D) résidu | ε'=t_f'/t_c' | δ/ε·(1+δ²/ε²) (τ*, 1er ordre)")
for eps in [0.2, 0.3, 0.5]:
    for tau in [0.05, 0.1, 0.3]:
        t = {1: eps, 2: 1.0, 3: tau}
        e = np.sort(np.abs(np.linalg.eigvalsh(chain(8, t))))
        x, y = fit_eff([e[0], e[2]])
        H = cluster3d(8, t); w = np.linalg.eigvalsh(H); o = np.argsort(np.abs(w)); low = np.sort(w[o[:64]]); aa = np.sort(np.abs(w))
        r = least_squares(lambda p: eff3d(p) - low, [x * 1.3, y * 0.7], bounds=([0, 0], [5, 5]))
        res = np.linalg.norm(eff3d(r.x) - low) / np.linalg.norm(low)
        delta = np.sqrt(eps ** 2 + 0.25) - 0.5; u = delta / eps
        print(" %.1f  %.2f |  %5.2f  |  %.6f        %.6f        | %.6f %.6f  %.1e | %.4f       | %.4f"
              % (eps, tau, aa[64] / aa[63], x, y, r.x[0], r.x[1], res, x / y, u * (1 + u * u)))
print(" (toute paire de niveaux 1D ±a, ±b est reproduite EXACTEMENT par la chaîne (x, y, x) : le test C2 « résidu < 5·10⁻⁴ » ne peut pas échouer)")

# point fixe τ*(ε) exact (tous ordres en τ) : résoudre x(τ)/y(τ) = ε
from scipy.optimize import brentq
def eps_prime(eps, tau):
    e = np.sort(np.abs(np.linalg.eigvalsh(chain(8, {1: eps, 2: 1.0, 3: tau})))); x, y = fit_eff([e[0], e[2]]); return x / y
print(" τ*(ε) exact :", {eps: round(brentq(lambda tau: eps_prime(eps, tau) - eps, 1e-4, 0.99), 4) for eps in [0.1, 0.2, 0.3, 0.5, 0.7]})

# --- deux pas : chaîne de 16 sites, t(1)=ε, t(2)=1, t(3)=τ, t(4)=T4 ; ce que donne vraiment l'autosimilarité ---
print("=== Deux pas (chaîne de 16 sites = amas 16×16×16 par séparabilité) ===")
def step(N, t):
    """Un pas : sites effectifs = bords des blocs de 4 ; amplitudes effectives (t'(1), t'(2), t'(3), ...) obtenues
    par ajustement d'une chaîne hiérarchique de N/2 sites sur les N/2 plus petits |e| (paramètres libres par profondeur)."""
    e = np.sort(np.linalg.eigvalsh(chain(N, t))); o = np.argsort(np.abs(e)); low = np.sort(e[o[:N // 2]])
    hmax = max(depth(x) for x in range(N // 2 - 1))
    def r(p): return np.sort(np.linalg.eigvalsh(chain(N // 2, {h: p[h - 1] for h in range(1, hmax + 1)}))) - low
    guess = [abs(low).min() * 1.7] + [t[h + 1] * 0.9 for h in range(1, hmax)]
    s = least_squares(r, guess, bounds=(0, 5))
    return {h: s.x[h - 1] for h in range(1, hmax + 1)}, np.linalg.norm(s.fun) / np.linalg.norm(low)
eps = 0.3
tau = brentq(lambda tt: eps_prime(eps, tt) - eps, 1e-4, 0.99)
for T4, lab in [(tau ** 2, "T4 = τ²  (décroissance géométrique)"), (1.0, "T4 = 1   (profil « alterné » de V198)"), (tau, "T4 = τ")]:
    t = {1: eps, 2: 1.0, 3: tau, 4: T4}
    tp, res = step(16, t)
    print(" ε=%.1f τ=%.4f %s : t'(1)/t'(2) = %.4f, t'(3)/t'(2) = %.4f (autosimilaire si = ε = %.4f et τ = %.4f) ; résidu %.1e"
          % (eps, tau, lab, tp[1] / tp[2], tp[3] / tp[2], eps, tau, res))
# profil autosimilaire cherché sur la profondeur 4 : t'(3)/t'(2) = τ
T4s = brentq(lambda T: (lambda tp: tp[3] / tp[2] - tau)(step(16, {1: eps, 2: 1.0, 3: tau, 4: T})[0]), 1e-4, 2.0)
print(" T4 autosimilaire = %.5f ; τ² = %.5f ; profil t(h) = %s" % (T4s, tau ** 2, [round(v, 4) for v in [eps, 1.0, tau, T4s]]))
