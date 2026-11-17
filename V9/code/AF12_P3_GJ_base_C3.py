#!/usr/bin/env python3
"""PROPOSITION P3 (05_audit/AUDIT_FABLE_PROPOSITIONS.md, criteres ecrits le 13/11/2026 a 20 h 30) : l'affectation de Georgi-Jarlskog
« une sur vingt » recomptee dans la base propre (v1, v_omega, v_omegabar) de C3 sous (P0), avec e = v1.
Classe minimale (V177/V178) : matrices symetriques 3x3 a exactement 3 entrees non nulles (parmi 11, 22, 33, 12, 13, 23), chaque entree portant un
operateur de coin pur : S (sur place, rapport lepton/quark 1) ou H (saut, rapport -3).  Les entrees (reelles, signes compris) sont ajustees
sur les masses des quarks d a 10^16 GeV (AHS25, Table 1, valeurs de V174/V178) ; les leptons charges s'en deduisent sans parametre.
Base C3 : indice 0 = v1 (axe fixe de C3), 1 = v_omega, 2 = v_omegabar ; e = v1 : l'etat propre le plus leger des leptons charges a son plus grand
poids sur v1 (> 1/2).  Affectations comptees a l'echange v_omega <-> v_omegabar pres (seule renumerotation qui fixe v1).
Reussite (critere P3) : X2 dans [2,6 ; 3,1] et X3 dans [0,030 ; 0,045].  Prediction ecrite : exactement une (celle de l'article).
Note : entrees reelles (comme V178) ; la phase invariante d'une texture « triangle » (12, 13, 23) n'est pas balayee (signale dans le journal).
Usage : python3 AF12_P3_GJ_base_C3.py        (numpy seul ; n'ecrit aucun fichier)
"""
import itertools as it
import numpy as np
rng = np.random.default_rng(20261114)
md = np.array([0.65e-5, 1.31e-4, 0.607e-2]); ml = np.array([2.6952e-6, 5.6779e-4, 0.9644e-2])   # V178 (AHS25 Table 1, 10^16 GeV)
def X(m, q): return m[2] / q[2], (m[1] / m[2]) / (q[1] / q[2]), (m[0] * m[2] / m[1] ** 2) / (q[0] * q[2] / q[1] ** 2)
X1o, X2o, X3o = X(ml, md)
POS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
def mat(P, x):
    M = np.zeros((3, 3))
    for (i, j), v in zip(P, x): M[i, j] = M[j, i] = v
    return M
def sv(M): return np.sort(np.abs(np.linalg.eigvalsh(M)))
def resid(P, x): return np.log(np.maximum(sv(mat(P, x)), 1e-300)) - np.log(md)
def fit(P, ntry=400):
    sols = []
    for _ in range(ntry):
        x = rng.choice([-1, 1], 3) * np.exp(rng.uniform(np.log(1e-6), np.log(1e-1), 3)); lam = 1e-3
        r = resid(P, x); c = r @ r
        for _ in range(300):
            J = np.empty((3, 3))
            for k in range(3):
                h = 1e-7 * max(abs(x[k]), 1e-12); xp = x.copy(); xp[k] += h; J[:, k] = (resid(P, xp) - r) / h
            A = J.T @ J; g = J.T @ r
            try: dx = -np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-30), g)
            except np.linalg.LinAlgError: break
            xn = x + dx; rn = resid(P, xn); cn = rn @ rn
            if cn < c: x, r, c, lam = xn, rn, cn, lam / 3
            else: lam *= 5
            if c < 1e-24 or lam > 1e12: break
        if np.max(np.abs(r)) < 1e-7 and not any(np.allclose(x, s, rtol=1e-5) for s in sols): sols.append(x)
    return sols
def canon(P, ops):
    best = None
    for p in ((0, 1, 2), (0, 2, 1)):
        key = tuple(sorted((tuple(sorted((p[i], p[j]))), o) for (i, j), o in zip(P, ops)))
        best = key if best is None or key < best else best
    return best
names = {0: 'v1', 1: 'vw', 2: 'vwb'}
fmt = lambda k: ', '.join(f'{names[a]}{names[b]}:{o}' for (a, b), o in k)
classes = {}; tri = 0
for P in it.combinations(POS, 3):
    sols = fit(list(P))
    if set(P) == {(0, 1), (0, 2), (1, 2)}: tri = len(sols)
    for ops in it.product('SH', repeat=3):
        k = canon(P, ops); rr = np.array([1. if o == 'S' else -3. for o in ops])
        rec = classes.setdefault(k, dict(fit=False, ev1=False, pass_P3=False, pass_V178=False, best=None))
        for x in sols:
            rec['fit'] = True
            Lm = mat(P, rr * x); w, V = np.linalg.eigh(Lm); i0 = np.argmin(np.abs(w))
            if abs(V[0, i0]) ** 2 <= 0.5: continue                         # e = v1 exige
            rec['ev1'] = True
            x1, x2, x3 = X(np.sort(np.abs(w)), sv(mat(P, x)))
            p3 = 2.6 <= x2 <= 3.1 and 0.030 <= x3 <= 0.045
            p178 = abs(x2 / X2o - 1) <= 0.3 and abs(x3 / X3o - 1) <= 0.3
            rec['pass_P3'] |= p3; rec['pass_V178'] |= p178
            if p3 or (rec['best'] is None): rec['best'] = (round(x1, 3), round(x2, 3), round(x3, 4), round(abs(V[0, i0]) ** 2, 3))
print(f'observe (10^16 GeV) : X1 = {X1o:.3f}, X2 = {X2o:.3f}, X3 = {X3o:.4f}')
nfit = sum(r['fit'] for r in classes.values()); nev1 = sum(r['ev1'] for r in classes.values())
P3 = [k for k, r in classes.items() if r['pass_P3']]; P178 = [k for k, r in classes.items() if r['pass_V178']]
print(f'affectations (positions x operateurs, a v_w <-> v_wb pres) : {len(classes)} ; ajustables sur les quarks d : {nfit} ; avec e = v1 : {nev1}')
print(f'reussissent (critere P3 : X2 in [2.6, 3.1], X3 in [0.030, 0.045]) : {len(P3)}')
for k in P3: print('   ', fmt(k), ' (X1, X2, X3, poids de e sur v1) =', classes[k]['best'])
print(f'pour comparaison, tolerance de V178 (+-30 %) : {len(P178)}')
for k in P178: print('   ', fmt(k), ' ', classes[k]['best'])
gj = canon([(0, 1), (1, 1), (2, 2)], ('S', 'H', 'S'))
print('affectation de l\'article (v1vw : S, vwvw : H, vwbvwb : S) :', 'reussit' if gj in P3 else 'ECHOUE', classes.get(gj, {}).get('best'))
print(f'texture triangle (12, 13, 23) : {tri} solutions reelles ; phase invariante non balayee')
verdict = ('confirme : exactement une affectation, celle de l\'article' if P3 == [gj] else
           'aucune : la phrase est a retirer' if not P3 else f'{len(P3)} affectations : le texte « exactly one out of twenty » est faux en base C3')
print('VERDICT P3 :', verdict)
