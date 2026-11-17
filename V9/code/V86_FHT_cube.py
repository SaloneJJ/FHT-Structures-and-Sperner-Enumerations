"""V86 : verification du cube des huit x (V85). JJ SALONE, FHT 2026"""
import itertools, time, sys
from collections import defaultdict, Counter
from V86_FHT_k26 import build_K, faces, betti_from_faces, maximal
K, I, II = build_K()
V = sorted(set().union(*K))
FK = faces(K)

def xn(S, c): return f'{S}_x{c}'
# regles de recollement : liste de (direction, [(v, w)]) : (t, v) ~ (t+e_dir, w)
def rules_cube(diag=3, pi=(0,1,2,3)):
    # carre sur {0,1,2,3} avec diagonale P_diag ; cotes = les deux autres appariements
    # V4 en Z2^2 : a1=1,a2=2,a3=3 ; H_k={0,a_k}
    sides = [k for k in (1,2,3) if k != diag]
    R = {0: [('R_z','L_z')] + [(xn('R',c), xn('L',pi[c])) for c in range(4)]}
    for j,k in enumerate(sides):
        other = sides[1-j]
        s = other           # translation par a_other : preserve H_other, echange les blocs de H_k
        plus = [c for c in range(4) if c not in (0,k)]   # coset de H_k
        R[j+1] = [(xn(S,c), xn(S,c ^ s)) for S in 'LR' for c in plus]
    return R

def glue(tiles, R):
    par = {}
    def f(a):
        while par.get(a, a) != a:
            par[a] = par.get(par[a], par[a]); a = par[a]
        return a
    def u(a,b):
        a,b = f(a),f(b)
        if a!=b: par[max(a,b)] = min(a,b)
    T = set(tiles)
    for t in tiles:
        for d, pairs in R.items():
            t2 = tuple(t[i] + (1 if i==d else 0) for i in range(len(t)))
            if t2 in T:
                for v,w in pairs: u((t,v),(t2,w))
    return f

def analyse(tiles, R, nom, homologie=True):
    print(f'--- {nom} : {len(tiles)} tuiles', flush=True)
    f = glue(tiles, R)
    lab = {(t,v): f((t,v)) for t in tiles for v in V}
    # 1. pas d'identification interne a une tuile
    ok = all(len(set(lab[(t,v)] for v in V)) == len(V) for t in tiles)
    cls = defaultdict(list)
    for (t,v),c in lab.items(): cls[c].append((t,v))
    sizes = Counter((v_name(cl), len(cl)) for cl in cls.values())
    print('  injectivite par tuile :', ok, '| classes (type, taille) :', dict(sizes))
    # 2. facettes, Sperner, connexite
    tileF = {t: [frozenset(lab[(t,v)] for v in F_) for F_ in K] for t in tiles}
    allF = [F_ for t in tiles for F_ in tileF[t]]
    print('  Sperner (aucune facette contenue dans une autre) :', len(maximal(allF)) == len(set(allF)) == len(allF))
    # 3. intersections de tuiles (complexes) : acycliques ; nerf
    tileFaces = {t: set(frozenset(lab[(t,v)] for v in s) for s in FK) for t in tiles}
    nerve = []; bad = 0; types = Counter()
    for r in range(2, len(tiles)+1):
        found = False
        for sub in itertools.combinations(tiles, r):
            inter = set.intersection(*(tileFaces[t] for t in sub))
            if not inter: continue
            found = True; nerve.append(frozenset(sub))
            b = betti_from_faces(set(tuple(sorted(s)) for s in inter))
            b = {k:v for k,v in b.items()}
            if b: bad += 1
            types[(r, max(map(len, inter)))] += 1
        if not found: break
    print('  intersections non vides (nb tuiles, dim max+1) :', dict(types), '| non acycliques :', bad)
    nb = betti_from_faces(set(tuple(sorted(map(str,s))) for s in faces(maximal(nerve + [frozenset([t]) for t in tiles]))))
    print('  nerf du recouvrement, Betti reduits :', nb)
    if homologie:
        t0 = time.time()
        fs = set()
        for t in tiles: fs |= set(tuple(sorted(map(str, s))) for s in tileFaces[t])
        b = betti_from_faces(fs)
        print(f'  Betti reduits du recollement : {b}  (attendu {{3: {12*len(tiles)}, 4: {36*len(tiles)}, 5: {36*len(tiles)}, 6: {12*len(tiles)}}})  [{len(fs)} faces, {time.time()-t0:.0f} s]', flush=True)

def v_name(cl):
    return '/'.join(sorted(set(v.split('_')[-1] if '_' in v else v for (_,v) in cl)))

if __name__ == '__main__':
    R = rules_cube()
    analyse([(0,0,0),(1,0,0)], R, 'chaine L-R, 2 tuiles')
    analyse([(0,0,0),(0,1,0)], R, 'chaine e1, 2 tuiles')
    analyse([(0,0,0),(1,0,0),(0,1,0),(1,1,0)], R, 'plaquette e0-e1')
    analyse([(0,0,0),(0,1,0),(0,0,1),(0,1,1)], R, 'plaquette e1-e2')
    box = [t for t in itertools.product(range(2), repeat=3)]
    analyse(box, R, 'bloc 2x2x2')
