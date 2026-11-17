"""Point 2 : lemme de V158 §1 (C^8 = 3 qubits, generateurs de Jordan-Wigner)."""
import numpy as np, itertools
X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1.,-1]).astype(complex); I=np.eye(2)
def op(lst):
    out=np.eye(1)
    for o in lst: out=np.kron(out,o)
    return out
G=[]
for j in range(3):
    for P in (X,Y):
        G.append(op([Z]*j+[P]+[I]*(2-j)))
ac=max(np.abs(G[a]@G[b]+G[b]@G[a]-2*(a==b)*np.eye(8)).max() for a in range(6) for b in range(6))
print('6 generateurs JW : ecart a {G_a,G_b}=2 delta =',ac)
G5=G[:5]; Om=G5[0]@G5[1]@G5[2]@G5[3]@G5[4]
print('Omega = G1..G5 : Omega^2 =',np.round((Om@Om)[0,0],6),'; max|[Omega,G_j]| =',max(np.abs(Om@g-g@Om).max() for g in G5))
print('G6 commute-t-il avec G1..G5 ? max|[G6,G_j]| =',max(np.abs(G[5]@g-g@G[5]).max() for g in G5),' ; max|{G6,G_j}| =',max(np.abs(G[5]@g+g@G[5]).max() for g in G5))
print('Omega proportionnel a G6 ? ecart min sur la phase :',min(np.abs(Om-c*G[5]).max() for c in (1,-1,1j,-1j)))
w=np.linalg.eigvals(Om/np.sqrt((Om@Om)[0,0])); print('valeurs propres de Omega normalise :',np.round(sorted(w.real),3))
# helicite du mode de paroi dans chaque secteur : sur le sous-espace -i G4 G5 = +1 (projecteur de paroi), G1G2G3 = ?
Pw=(np.eye(8)+(-1j*G5[3]@G5[4]))/2
for s in (+1,-1):
    Ps=(np.eye(8)+s*Om/np.sqrt((Om@Om)[0,0]))/2; P=Ps@Pw
    v=np.linalg.svd(P)[0][:,:np.linalg.matrix_rank(P)]
    h=v.conj().T@(G5[0]@G5[1]@G5[2])@v
    print(f'secteur Omega = {s:+d} : dim mode de paroi {v.shape[1]}, G1G2G3 restreint = ',np.round(np.diag(h),6))
