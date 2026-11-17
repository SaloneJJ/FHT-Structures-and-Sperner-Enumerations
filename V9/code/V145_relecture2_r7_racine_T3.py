# V146 §0 / V145 : "fermeture T3 de la tuile racine par la meme regle de faces" : la racine recollee a elle-meme
# identifie ses 8 coins en un seul point. Que devient l'homologie de la copie de K26 portee par la racine ?
import time
from V86_FHT_k26 import build_K, faces, betti_from_faces
from V87_FHT_stella import rules_stella
K, I, II = build_K()
par = {}
def f(a):
    while par.get(a, a) != a: a = par[a]
    return a
def u(a, b):
    a, b = f(a), f(b)
    if a != b: par[max(a, b)] = min(a, b)
for k, pairs in rules_stella().items():
    for v, w in pairs: u(v, w)
classes = {}
for v in sorted(set().union(*K)): classes.setdefault(f(v), []).append(v)
print('classes non triviales :', [c for c in classes.values() if len(c) > 1])
Kq = set(frozenset(f(v) for v in F) for F in K)
t = time.time(); fs = faces(Kq)
print('Betti reduits (Z/2) de K26 recolle a lui-meme (fermeture T3 de la racine) :', betti_from_faces(fs), f'({time.time()-t:.0f} s)')
print('rappel K26 :', {3: 12, 4: 36, 5: 36, 6: 12})
