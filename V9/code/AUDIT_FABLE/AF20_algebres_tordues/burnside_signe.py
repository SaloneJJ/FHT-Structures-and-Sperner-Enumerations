"""Condition necessaire (formule de Hopf) pour realiser un QG-module gradue M comme homologie rationnelle d'un G-complexe fini :
la classe d'Euler sum (-1)^k [M_k] (+ correction du simplexe vide) est dans le Z-module engendre par les caracteres Ind_H^G(eps),
H <= G, eps : H -> {+-1} (modules de chaines orientees). Calcul de ce reseau et test d'appartenance."""
from twisted_alg import *
from fractions import Fraction
import itertools, numpy as np
def subgroups(G):
    idx={g:i for i,g in enumerate(G)}; subs={}
    for a in G:
        for b in G:
            H=frozenset(closure([a,b])); subs[H]=1
    return list(subs)
def homs_pm(H):
    H=sorted(H); e=tuple(range(len(H[0])))
    # generateurs minimaux par fermeture gloutonne
    gens=[]; cur={e}
    for h in H:
        if h not in cur: gens.append(h); cur=set(closure(gens))
    out=[]
    for signs in itertools.product([1,-1],repeat=len(gens)):
        val={e:1}; fr=[e]; ok=True
        while fr and ok:
            x=fr.pop()
            for g,s in zip(gens,signs):
                y=mul(x,g); v=val[x]*s
                if y in val:
                    if val[y]!=v: ok=False; break
                else: val[y]=v; fr.append(y)
        if ok and len(val)==len(H): out.append(val)
    return out
def induced(G, H, eps, reps):
    Hs=set(H); res=[]
    for r in reps:
        tot=0
        for x in G:
            y=mul(mul(x,r),inv(x))
            if y in Hs: tot+=eps[y]
        res.append(Fraction(tot,len(H)))
    return res
def lattice_basis(vectors):
    """base d'Hermite (lignes) du Z-module engendre (entiers)"""
    M=[list(map(int,v)) for v in vectors]; n=len(M[0]); B=[]; row=0
    import math
    for col in range(n):
        piv=[i for i in range(len(M)) if M[i][col]!=0]
        while len(piv)>1:
            piv.sort(key=lambda i:abs(M[i][col])); p=piv[0]
            for i in piv[1:]:
                q=M[i][col]//M[p][col]; M[i]=[a-q*b for a,b in zip(M[i],M[p])]
            piv=[i for i in range(len(M)) if M[i][col]!=0]
        if piv:
            p=piv[0]; B.append(M[p]); M.pop(p)
    return B
def in_lattice(B, v):
    v=list(map(int,v)); n=len(v)
    for b in B:
        c=next(i for i in range(n) if b[i]!=0)
        if v[c]%b[c]!=0: return False
        q=v[c]//b[c]; v=[a-q*x for a,x in zip(v,b)]
    return all(x==0 for x in v)
def analyse(G, name, targets=None):
    C,cl,idx,reps,X=char_table(G)
    sizes=[len(c) for c in C]; k=len(C); e=tuple(range(len(G[0])))
    e0=[l for l in range(k) if reps[l]==e][0]
    gens=[]
    for H in subgroups(G):
        for eps in homs_pm(H): gens.append(induced(G,sorted(H),eps,reps))
    B=lattice_basis([[int(x) for x in v] for v in gens])
    # caracteres rationnels irreductibles : sommes des conjugues de Galois (approx : chi + conj(chi) si non reel)
    irr=[]; used=set()
    for a in range(k):
        if a in used: continue
        chi=X[a]
        if np.allclose(chi.imag,0,atol=1e-6): v=chi.real; used.add(a)
        else:
            b=[bb for bb in range(k) if np.allclose(X[bb],np.conj(chi),atol=1e-6)][0]; v=(chi+X[b]).real; used|={a,b}
        irr.append((int(round(v[e0])), np.round(v).astype(int)))
    print(f'{name} : |G|={len(G)}, {len(gens)} generateurs, rang du reseau {len(B)} / {k}')
    for d,v in irr:
        print(f'   irreductible rationnel de degre {d:2d} : dans le reseau (Moore realisable, condition necessaire) ? {in_lattice(B,v)}')
    return B, irr, reps
if __name__=='__main__':
    analyse(closure([(1,0,2,3),(1,2,3,0)]), 'S4')
    vecs=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
    def mp(M): return tuple(vecs.index(((M[0][0]*a+M[0][1]*b)%3,(M[1][0]*a+M[1][1]*b)%3)) for a,b in vecs)
    analyse(closure([mp([[1,1],[0,1]]),mp([[0,1],[2,0]]),mp([[2,0],[0,1]])]), 'GL(2,3)')
    Q8=closure([mp([[0,1],[2,0]]),mp([[1,1],[1,2]])]); analyse(Q8,'Q8 (controle classique)')
    analyse(closure([(1,0,2,3,4),(1,2,3,4,0)]), 'S5')
