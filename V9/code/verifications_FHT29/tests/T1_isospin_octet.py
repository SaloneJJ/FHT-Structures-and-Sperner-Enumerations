# T1 (Fable, 16/11/2026) -- L'indice de saveur E de l'octet est-il l'indice d'isospin ?
# Partie algebrique, numpy seul, quelques secondes. Conventions de V335/V340 (Xg, Yg, G5, Am, Bm, G4, FG).
# Criteres : CRITERES_T1_T2.md (ecrits avant ce calcul ; la lecture du journal V335 est declaree post hoc).
import numpy as np
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.diag([1.,-1]).astype(complex); I2=np.eye(2)
k3=lambda a,b,c:np.kron(np.kron(a,b),c)
Xg=[k3(sx,I2,I2),k3(sz,sx,I2),k3(sz,sz,sx)]; Yg=[k3(sy,I2,I2),k3(sz,sy,I2),k3(sz,sz,sy)]
G5=sum(Xg)/np.sqrt(3); A=(Xg[0]-Xg[1])/np.sqrt(2); B=(Xg[0]+Xg[1]-2*Xg[2])/np.sqrt(6); G4=G5@Yg[0]@Yg[1]@Yg[2]; F=1j*A@B; I8=np.eye(8)
comm=lambda X,Y:np.abs(X@Y-Y@X).max(); acomm=lambda X,Y:np.abs(X@Y+Y@X).max()
TF,TC=0.4,1.0
def H(k,tf=TF,tc=TC): return sum((tf+tc*np.cos(k[d]))*Xg[d]+tc*np.sin(k[d])*Yg[d] for d in range(3))
K=np.array([np.pi]*3)
print("== (a) algebre de saveur au point K")
rows=[np.kron(I8,X)-np.kron(X.T,I8) for X in [G5]+Yg]; s=np.linalg.svd(np.vstack(rows),compute_uv=False)
print("   dim_C du commutant de {G5, alpha_d} :",int((s<1e-9).sum()),"(attendu 4)")
quat=[I8,1j*G4@A,1j*G4@B,A@B]
print("   {1, iG4A, iG4B, AB} : carres des trois unites = -1 :",[bool(np.allclose(q@q,-I8)) for q in quat[1:]],
      "| anticommutent :",max(acomm(quat[i],quat[j]) for i in (1,2,3) for j in (1,2,3) if i<j)<1e-12,
      "| commutent a G5, alpha_d, G4 :",max(comm(q,X) for q in quat for X in [G5,G4]+Yg)<1e-12)
print("== (b) structure reelle de l'octet : kappa (fonctions reelles sur les huit sites)")
Sg=G5@A@B@Yg[0]@Yg[1]@Yg[2]
print("   H(k) reel en position (conj H(-k) = H(k)) :",all(np.abs(H(-k).conj()-H(k)).max()<1e-12 for k in np.random.default_rng(1).uniform(-3,3,(5,3))),
      "| les quatre quaternions sont des matrices reelles :",[bool(np.abs(q.imag).max()<1e-12) for q in quat])
print("   donc H_saveur = {1, iG4A, iG4B, AB} est la forme reelle fixee par kappa ; F = iAB est imaginaire pur :",bool(np.abs(F.real).max()<1e-12))
print("== (c) chiralite et algebre gauche-droite native")
Gch=1j*Yg[0]@Yg[1]@Yg[2]
print("   gamma = i Y1Y2Y3 : carre 1 :",bool(np.allclose(Gch@Gch,I8)),"| commute aux alpha_d :",max(comm(Gch,Y) for Y in Yg)<1e-12,"| anticommute a G5 et G4 :",max(acomm(Gch,X) for X in (G5,G4))<1e-12,
      "| commute a H_saveur :",max(comm(Gch,q) for q in quat)<1e-12,"| i G4 G5 = -gamma :",bool(np.allclose(1j*G4@G5,-Gch)))
PL,PR=(I8+Gch)/2,(I8-Gch)/2
rk=lambda gens:np.linalg.matrix_rank(np.array([np.concatenate([g.real.ravel(),g.imag.ravel()]) for g in gens]))
print("   dim_R de P_L H + P_R H :",rk([PL@q for q in quat]+[PR@q for q in quat]),"(H + H = 8) ; dim_R de P_L H + P_R {1, AB} :",rk([PL@q for q in quat]+[PR,PR@(A@B)]),"(H + C = 6)")
print("   sur P_R : (P_R AB)^2 = -P_R :",bool(np.allclose((PR@A@B)@(PR@A@B),-PR)),"; F = iAB y vaut +1 et -1 avec les multiplicites",[int(round(np.trace(PR@(I8+e*F)/2).real)) for e in (1,-1)])
print("== (d) le terme (m + mu F) G5 est-il une masse haut/bas ?")
for m,mu in ((0.3,0.0),(0.3,0.1),(0.0,0.2)):
    ev=np.sort(np.abs(np.linalg.eigvalsh(H(K,1.0,1.0)+(m*I8+mu*F)@G5))); print(f"   t_f=t_c, k=K, m={m}, mu={mu} : |E| =",np.round(ev,6),"-> m+mu et |m-mu|, quatre fois chacun")
print("== (e) conservation loin du cone (t_f = t_c) : ordre q^2 attendu")
for q in (0.0,0.05,0.1,0.2,0.4):
    k=K+np.array([q,0.7*q,-0.4*q]); Hk=H(k,1.0,1.0)
    print(f"   |q|~{q:4.2f} : ||[H,F]|| = {comm(Hk,F):.2e} ; ||[H, P_L iG4A]|| = {comm(Hk,PL@quat[1]):.2e} ; ||[H, P_R AB]|| = {comm(Hk,PR@A@B):.2e}")
print("== (f) comptage (arithmetique)")
print("   lecture actuelle de A3 : 96 etats de F x 2 (E) = 192 branches -> Majorana 96 -> masse de saveur 48")
print("   lecture testee : (2 chiralites) x (1+3 lepton-couleur) x (3 generations) = 24 secteurs x 2 (E = isospin) = 48 fermions de Weyl, sans masse de saveur")
