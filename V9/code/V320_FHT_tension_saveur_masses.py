# V320 (FHT_25) — Tension saveur (V68 : leptons chargés en base C3 complexe, e = direction S3-invariante)
# contre masses (V177-V178 : base de chaîne réelle S3 > C2 > 1, v1 = S3-invariante = la plus lourde).
# Colonne PMNS fixée par le résiduel neutrino (0 j) dans chacune des deux bases. Rep 3 de S4 = {x in R^4, somme nulle}. numpy seul.
import numpy as np
w=np.exp(2j*np.pi/3)
v1=np.array([3,-1,-1,-1])/np.sqrt(12)           # fixe par S3 = Stab(g0)
v2=np.array([0,1,1,-2])/np.sqrt(6)              # fixe par C2 = <(12)>
v3=np.array([0,1,-1,0])/np.sqrt(2)              # impair sous (12)
c1=np.array([0,1,w,w*w])/np.sqrt(3); c2=c1.conj()   # base propre complexe de C3 = <(123)>
for j in (1,2,3):
    wv=np.zeros(4); wv[0]=1; wv[j]=-1; wv/=np.sqrt(2)   # vecteur impair de la transposition (0 j)
    C3col=[abs(np.vdot(b,wv))**2 for b in (v1,c1,c2)]
    CHcol=[abs(np.vdot(b,wv))**2 for b in (v1,v2,v3)]
    print(f"(0 {j}) : base C3 (v1, w, w2) -> {np.round(C3col,4)} ; base de chaîne (v1, v2, v3) -> {np.round(CHcol,4)}")
print("V68 : e = v1 (2/3 = |U_e1|^2, TM1). V178 : v1 = la plus lourde = tau ⇒ |U_tau1|^2 = 2/3 (données : ~0,1–0,25).")
# T1 : Y_e réel, Y_e Y_e^T commute avec C3 = <(123)> dans la rep 3 ⇒ valeurs singulières (mu, tau) égales ?
P=np.array([[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]],float)   # (123) sur g1,g2,g3, g0 fixe
B=np.array([v1,v2,v3]).T                                          # base orthonormée réelle de la rep 3
R=B.T@P@B                                                        # rotation d'ordre 3 dans la rep 3
rng=np.random.default_rng(320); worst=0
for _ in range(2000):
    Y=rng.normal(size=(3,3))
    Yc=sum(np.linalg.matrix_power(R,k)@Y@np.linalg.matrix_power(R,-k) for k in range(3))/3  # moyenne : [R, Yc]=0
    s=np.sort(np.linalg.svd(Yc,compute_uv=False))
    S=Yc@Yc.T
    # valeurs singulières sur le plan perpendiculaire à v1 (deux d'entre elles) :
    ev=np.linalg.eigvalsh(S); worst=max(worst, min(abs(ev[0]-ev[1]),abs(ev[1]-ev[2]),abs(ev[0]-ev[2])))
print("T1 : Y réel commutant avec C3 : plus petit écart entre deux valeurs propres de Y Y^T (max sur 2000 tirages) =",worst)
# contrôle : Y complexe C3-invariant ⇒ non dégénéré
Y=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
Yc=sum(np.linalg.matrix_power(R,k)@Y@np.linalg.matrix_power(R,-k) for k in range(3))/3
print("contrôle, Y complexe C3-invariant : valeurs singulières",np.round(np.linalg.svd(Yc,compute_uv=False),4))
