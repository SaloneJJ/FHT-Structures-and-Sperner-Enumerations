import numpy as np, math
def cosd(s13sq, s23sq):
    s = math.sqrt(s13sq); t = math.asin(math.sqrt(s23sq))
    return -(1-5*s13sq)*(math.cos(2*t)/math.sin(2*t))/(2*math.sqrt(2)*s*math.sqrt(1-3*s13sq))
for nom, (a, b) in (('NO', (0.435, 0.585)), ('IO', (0.440, 0.584))):
    vals = [math.degrees(math.acos(max(-1, min(1, cosd(x, y))))) for x in np.linspace(0.02077, 0.02419, 21) for y in np.linspace(a, b, 301)]
    print(f'{nom}: |δ| ∈ [{min(vals):.1f}°, {max(vals):.1f}°] sur s13² ∈ [0.02077, 0.02419], s23² ∈ [{a}, {b}]')
for x in (0.02077, 0.02248, 0.02419):
    k = 2*math.sqrt(2)*math.sqrt(x)*math.sqrt(1-3*x)/(1-5*x)   # |cot 2θ23| pour |cos δ| = 1
    t1 = 0.5*math.atan(1/k); t2 = math.pi/2 - t1
    print(f's13²={x}: δ = 180° ⇔ sin²θ23 = {math.sin(t1)**2:.4f} ; δ = 0 ⇔ sin²θ23 = {math.sin(t2)**2:.4f}')
# contrôle croisé sur matrices explicites
Utbm = np.array([[2/math.sqrt(6), 1/math.sqrt(3), 0], [-1/math.sqrt(6), 1/math.sqrt(3), -1/math.sqrt(2)], [-1/math.sqrt(6), 1/math.sqrt(3), 1/math.sqrt(2)]])
pts = []
for th in np.linspace(0, math.pi, 721):
    for ph in np.linspace(-math.pi, math.pi, 361):
        R = np.array([[1, 0, 0], [0, math.cos(th), math.sin(th)*np.exp(-1j*ph)], [0, -math.sin(th)*np.exp(1j*ph), math.cos(th)]])
        U = Utbm @ R
        s13 = abs(U[0, 2])**2
        if not (0.02077 <= s13 <= 0.02419): continue
        s12 = abs(U[0, 1])**2/(1-s13); s23 = abs(U[1, 2])**2/(1-s13)
        if not (0.435 <= s23 <= 0.585): continue
        J = np.imag(U[0, 0]*U[1, 1]*np.conj(U[0, 1])*np.conj(U[1, 0]))
        S12, S23, S13 = math.sqrt(s12), math.sqrt(s23), math.sqrt(s13); c12, c23, c13 = math.sqrt(1-s12), math.sqrt(1-s23), math.sqrt(1-s13)
        cd = (abs(U[1, 0])**2 - (S12*c23)**2 - (c12*S23*S13)**2)/(2*S12*c23*c12*S23*S13)
        sd = J/(c12*S12*c23*S23*c13**2*S13)
        pts.append(abs(math.degrees(math.atan2(sd, cd))))
print(f'contrôle matrices explicites ({len(pts)} points, NO) : |δ| ∈ [{min(pts):.1f}°, {max(pts):.1f}°]')
