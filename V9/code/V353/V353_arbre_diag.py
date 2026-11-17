"""V353 (d2) diagnostic, ajouté APRÈS le résultat de l'arbre N = 2 (forte harmonique cos 2φ, absente sur le réseau périodique) :
l'arbre de V336 est-il invariant sous la permutation cyclique des axes (C3 de T) ? Test terme à terme : [H, P_g ⊗ U_g] = 0 ?"""
import sys, numpy as np
sys.argv = ['x', '2', '0', '1', '1', '4', '0', 'unif', '1.0']
exec(open('V336_FHT_flot.py').read().split("t0 = time.time(); H0")[0])
import rvA_common as R
Y123 = Yg[0]@Yg[1]@Yg[2]
def O(p): return (np.cos(p)*Am + np.sin(p)*Bm)@Y123
for name, g in (('C3 (x,y,z)->(y,z,x)', np.eye(3)[[1, 2, 0]]), ('C3 inverse', np.eye(3)[[2, 0, 1]]), ('C2z', np.diag([-1., -1., 1.]))):
    U, ok = R.lift(g)
    # permutation des cellules : x -> g x (mod taille du niveau)
    P = np.zeros((Nc, Nc))
    for n in range(N + 1):
        s = sides[n]
        for I in range(s):
            for J in range(s):
                for K in range(s):
                    v = g@np.array([I, J, K]); P[lin(n, *[int(t) % s for t in v]), lin(n, I, J, K)] = 1
    W = np.kron(P, U)
    T, S = ops(0.3); Id = sp.identity(Nc, format='csr')
    terms = {}
    for d in range(3):
        terms[f'sin Y{d+1}'] = kr((T[d] - T[d].getH())/2j, Yg[d]).toarray()
        terms[f'cos X{d+1}'] = (kr(Id, Xg[d]) - kr((T[d] + T[d].getH())/2, Xg[d])).toarray()
    terms['vert Γ4'] = kr((S - S.getH())/2j, G4).toarray(); terms['vert Wilson Γ5'] = kr(2*Id - S - S.getH(), G5).toarray(); terms['μ FΓ5'] = kr(Id, FG).toarray()
    Hs = {k: v for k, v in terms.items()}
    Hsum_h = sum(v for k, v in terms.items() if k.startswith('sin')) ; Hsum_c = sum(v for k, v in terms.items() if k.startswith('cos'))
    out = [f"{k}: {np.abs(W@v@W.conj().T - v).max():.1e}" for k, v in (('Σ sin', Hsum_h), ('Σ cos', Hsum_c), ('Γ4 vert', terms['vert Γ4']), ('Wilson vert', terms['vert Wilson Γ5']), ('μFΓ5', terms['μ FΓ5']))]
    print(name, '|', ' ; '.join(out))
    # action de U sur le doublet O
    M = [[np.real(np.trace(U@O(a)@U.conj().T@O(b)))/8 for a in (0, np.pi/2)] for b in (0, np.pi/2)]
    print('   U O U† dans (O_A, O_B) :', np.round(M, 3))
# vertical : décomposer S en enfants c ; S est-il symétrique sous permutation des axes ?
T, S = ops(0.3); Sd = S.toarray()
print('S (vertical) : valeurs distinctes non nulles :', np.unique(np.round(Sd[np.abs(Sd) > 1e-12], 6)))
