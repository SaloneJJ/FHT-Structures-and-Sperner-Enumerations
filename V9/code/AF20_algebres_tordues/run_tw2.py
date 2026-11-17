from twisted_alg import *
def sgn_on(p, pts): return sign(tuple(sorted(pts).index(p[i]) for i in sorted(pts)))
# S4 x S5 sur 4+5 points
g=[]
t=list(range(9)); t[0],t[1]=t[1],t[0]; g.append(tuple(t)); c=list(range(9)); c[0:4]=[1,2,3,0]; g.append(tuple(c))
t=list(range(9)); t[4],t[5]=t[5],t[4]; g.append(tuple(t)); c=list(range(9)); c[4:9]=[5,6,7,8,4]; g.append(tuple(c))
G=closure(g); print(len(G))
s4=lambda p: sgn_on(p,range(4)); s5=lambda p: sgn_on(p,range(4,9))
for nm,eps in [('S4xS5, sgn x 1',s4),('S4xS5, sgn x sgn',lambda p:s4(p)*s5(p)),('S4xS5, 1 x sgn',s5)]: decompose(G,eps,nm)
G6=closure([(1,0,2,3,4,5),(1,2,3,4,5,0)]); decompose(G6, sign, 'S6, sgn')
