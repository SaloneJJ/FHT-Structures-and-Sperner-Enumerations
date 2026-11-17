# Extension de (c2) : meme formule, mais |U_e1|^2 = c12^2 c13^2 laisse a sa valeur mesuree (niveau "realite seule")
from itertools import product
from math import sqrt
D21=(6.92e-5,8.05e-5); IO=dict(d3=(2.421e-3,2.547e-3),s13=(0.02060,0.02409)); NO=dict(d3=(2.451e-3,2.578e-3),s13=(0.02030,0.02388))
def io(m3,d21,d3,s13,s12,e2,e3):
    m2=sqrt(m3*m3+d3); m1=sqrt(m2*m2-d21); c13=1-s13
    return abs((1-s12)*c13*m1+e2*s12*c13*m2+e3*s13*m3)*1e3
def no(m1,d21,d3,s13,s12,e2,e3):
    m2=sqrt(m1*m1+d21); m3=sqrt(m1*m1+d3); c13=1-s13
    return abs((1-s12)*c13*m1+e2*s12*c13*m2+e3*s13*m3)*1e3
for name,rng in (("NuFIT 3s [0.275,0.345]",(0.275,0.345)),("JUNO 3s 0.3092+-0.0261",(0.2831,0.3353)),("JUNO 1s",(0.3005,0.3179)),("TM1 exact",None)):
    out=[]
    for f,par,m0 in ((io,IO,0.0),(io,IO,0.010),(no,NO,0.0)):
        for e2 in (1,-1):
            vals=[]
            for d21,d3,s13,e3 in product(D21,par['d3'],par['s13'],(1,-1)):
                ss=[1-(2/3)/(1-s13)] if rng is None else rng
                for s12 in ss: vals.append(f(m0,d21,d3,s13,s12,e2,e3))
            out.append(f"{min(vals):5.1f}-{max(vals):5.1f}")
    print(f"{name:26s} IO m3=0: +:{out[0]} -:{out[1]} | IO m3=10: +:{out[2]} -:{out[3]} | NO m1=0: +:{out[4]} -:{out[5]}")
