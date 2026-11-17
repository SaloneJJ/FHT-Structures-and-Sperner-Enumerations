"""V357 (relecture) : chaîne de Wilson–Kaplan en profondeur, code indépendant du relecteur (numpy seul).
Modèle réduit au point du cône (impulsion latérale nulle) : niveaux n = 0 (racine) .. N (bord fin), saut parent -> enfant d'amplitude 1
(combinaison symétrique des 8 enfants), H_s,f = (S - S†)/2i ⊗ g4 + [(r/2)(2 - S - S†) + m_K + s f μ η_n] ⊗ g5 (secteur γ_F = s, goût f),
g4 = σ_y, g5 = σ_x, chiralité c = σ_z (anticommute avec g4, g5 ; commute avec les α latéraux).
Sorties : modes nuls (noyau de Q, Q† du bloc chiral), leur chiralité, leur goût, leur profil en profondeur.
Usages : V340 (bit par niveau), V341 (borne ε), V345 (pont : m(h) = ε_D (1 − λ′²) λ′^(2h))."""
import numpy as np, itertools
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex)

def chain(N, mK, mu, eta, s=+1, f=+1, r=1.0):
    L = N + 1; S = np.zeros((L, L)); 
    for n in range(N): S[n + 1, n] = 1.0                       # enfant <- parent
    M = (r / 2) * (2 * np.eye(L) - S - S.T) + np.diag(mK + s * f * mu * np.asarray(eta, float))   # s·M = s(W + m_K) + f μ η (terme de goût nu)
    return np.kron((S - S.T) / 2j, sy) + s * np.kron(M, sx)          # secteur s : partie paire habillée par s (convention « cos »)

def zero_modes(N, mK, mu, eta, s=+1, tol=1e-6):
    out = []
    for f in (+1, -1):
        H = chain(N, mK, mu, eta, s, f)
        E, V = np.linalg.eigh(H)
        idx = np.where(np.abs(E) < tol)[0]
        if len(idx):
            Vz = V[:, idx]; C = Vz.conj().T @ np.kron(np.eye(N + 1), sz) @ Vz   # diagonaliser la chiralité dans le noyau
            c, U = np.linalg.eigh(C); Vz = Vz @ U
            for j in range(len(idx)):
                prof = (np.abs(Vz[:, j]) ** 2).reshape(N + 1, 2).sum(1)
                out.append((f, int(round(c[j])), prof))
    return out

def gaps(N, mK, mu, eta, s=+1):
    return {f: np.sort(np.abs(np.linalg.eigvalsh(chain(N, mK, mu, eta, s, f))))[:2] for f in (+1, -1)}

if __name__ == '__main__':
    print('=== R340 : motifs η par niveau, μ = r = 1, m_K = 0 (R-unif), secteur s = +1')
    for N in (3, 10):
        motifs = [tuple([1] + list(m)) for m in itertools.product((1, -1), repeat=N)] if N == 3 else \
                 [tuple([1] * N + [1]), tuple([1] * N + [-1]), tuple([1] * (N - 1) + [-1, 1]), tuple([1, -1] * (N // 2) + [1] * (N + 1 - 2 * (N // 2))),
                  tuple([1] * 5 + [-1] * (N - 4))]
        for eta in motifs:
            zm = zero_modes(N, 0.0, 1.0, eta)
            desc = []; net_bord = {k: 0 for k in range(1, N + 1)}
            for f, c, p in zm:
                desc.append(f'f={f:+d} c={c:+d} prof={np.round(p, 2).tolist() if N == 3 else int(np.argmax(p))}')
                for k in range(1, N + 1): net_bord[k] += c * p[k:].sum()
            g = gaps(N, 0.0, 1.0, eta)
            print(f'N={N} η={"".join("+" if e > 0 else "-" for e in eta)} : {len(zm)} modes nuls ; ' + ' | '.join(desc))
            print(f'      chiralité nette de la région [k, N] pour k = 1..N : {[round(net_bord[k], 3) for k in range(1, N + 1)]} ; 2 plus petites |E| par goût : '
                  f'{ {f: np.round(v, 3).tolist() for f, v in g.items()} }')
    print('\n=== R341 : chaîne N = 40, μ = r = t_c = 1, m_K = √3(ε − 1) ; mode de bord (région [30, 40]) : goût, chiralité, λ')
    for eps in (1.0, 0.8, 0.6, 0.5, 0.45, 0.40, 0.3, 0.2, 0.1, 0.0, -0.5, -0.8):
        mK = np.sqrt(3) * (eps - 1); N = 40; eta = np.ones(N + 1)
        res = []
        for f in (+1, -1):
            E, V = np.linalg.eigh(chain(N, mK, 1.0, eta, +1, f)); i = np.argsort(np.abs(E))[:2]
            for j in i:
                p = (np.abs(V[:, j]) ** 2).reshape(N + 1, 2).sum(1)
                if p[30:].sum() > 0.5 and abs(E[j]) < 1e-3:
                    c = np.real(np.vdot(V[:, j], np.kron(np.eye(N + 1), sz) @ V[:, j]))
                    lam = np.sqrt(p[N - 1] / p[N]) if p[N] > 1e-300 else np.nan
                    res.append(f'f={f:+d} (|E|={abs(E[j]):.1e}, chiralité ⟨c⟩={c:+.3f}, λ mesuré={lam:.3f})')
        ls, lt = abs(1 + (mK - 1)), abs(1 + (mK + 1))
        print(f'ε={eps:+.2f} m_K={mK:+.3f} λ_surv={ls:.3f} λ_jum={lt:.3f} : ' + ('; '.join(res) if res else 'aucun mode de bord'))
