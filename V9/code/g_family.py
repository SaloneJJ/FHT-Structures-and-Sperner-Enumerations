"""H1.3 / G1, G2 : famille K^(n) = (A^L * B^R * {n points}) u (B^L * A^R * squelette_(n-3)(simplexe a n sommets)).
n = 2 : squelette_(-1) = complexe vide (G1 = K25 sans pendant) ; n = 4 : K26 ; n = 6 : G2."""
import sys, json, itertools, time
sys.path.insert(0, '.')
from fht_tools import *
d = json.load(open('V25_K3_pieces.json'))
A_ = [frozenset(f) for f in maximal(d['A'])]
B_ = [frozenset(f) for f in maximal([b for b in d['B'] if not any(v.startswith('y') for v in b)])]
def Kn(n):
    G = ['g%d' % i for i in range(n)]
    I = [frozenset(['L_' + v for v in a] + ['R_' + v for v in b] + [g]) for a in A_ for b in B_ for g in G]
    sk = [frozenset(c) for c in itertools.combinations(G, n - 2)] if n >= 3 else [frozenset()]
    II = [frozenset(['L_' + v for v in b] + ['R_' + v for v in a]) | s for b in B_ for a in A_ for s in sk]
    K = maximal(I + II)
    return K, sorted(set(v for f in K for v in f))
if __name__ == '__main__':
    for n in [int(a) for a in sys.argv[1:]] or [2, 4, 6]:
        t0 = time.time()
        K, V = Kn(n)
        print(f"n = {n} : feuilles {len(V)}, facettes {len(K)}, Betti Z/2 {betti_z2(K)}, |Aut| = {aut_order(V, K)}  ({time.time()-t0:.0f} s)", flush=True)
