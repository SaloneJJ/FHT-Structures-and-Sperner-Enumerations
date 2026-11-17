from twisted_alg import *
def s4(off,n): # generateurs de S4 sur {off..off+3} dans S_n
    t=list(range(n)); t[off],t[off+1]=t[off+1],t[off]
    c=list(range(n)); c[off:off+4]=[off+1,off+2,off+3,off]
    return [tuple(t),tuple(c)]
G=closure(s4(0,4)); decompose(G, sign, 'S4, sgn (controle A_F)')
G5=closure([(1,0,2,3,4),(1,2,3,4,0)]); decompose(G5, sign, 'S5, sgn'); decompose(G5, lambda p:1, 'S5, trivial')
# couronne S4 wr C2 sur 8 points
sw=tuple([4,5,6,7,0,1,2,3])
W=closure(s4(0,8)+s4(4,8)+[sw]); print(len(W))
sL=lambda p: sign(tuple(p[i] for i in range(4)))
def base_sign(p):  # signe de la partie base (gauche*droite) : signe de la permutation sur 8 points, corrige du swap
    return sign(p)
def swap_part(p): return -1 if p[0]>=4 else 1
for nm,eps in [('S4 wr C2, sgn8',sign),('S4 wr C2, sgn8*swap',lambda p:sign(p)*swap_part(p)),('S4 wr C2, swap',swap_part),('S4 wr C2, trivial',lambda p:1)]:
    decompose(W,eps,nm)
