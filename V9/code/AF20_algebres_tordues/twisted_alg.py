"""Algebre de groupe tordue R^eps[G] (eps : G -> +-1) : decomposition en facteurs simples via la table de caracteres (algorithme de Burnside)."""
import itertools, numpy as np, sys
def closure(gens):
    n=len(gens[0]); e=tuple(range(n)); G={e}; fr=[e]
    while fr:
        x=fr.pop()
        for g in gens:
            y=tuple(x[g[i]] for i in range(n))
            if y not in G: G.add(y); fr.append(y)
    return sorted(G)
def mul(a,b): return tuple(a[b[i]] for i in range(len(a)))
def inv(a):
    r=[0]*len(a)
    for i,x in enumerate(a): r[x]=i
    return tuple(r)
def classes(G):
    idx={g:i for i,g in enumerate(G)}; cl=[-1]*len(G); C=[]
    for g in G:
        if cl[idx[g]]>=0: continue
        c={mul(mul(h,g),inv(h)) for h in G}
        for x in c: cl[idx[x]]=len(C)
        C.append(sorted(c))
    return C, cl, idx
def char_table(G):
    C,cl,idx=classes(G); k=len(C); N=len(G)
    reps=[c[0] for c in C]
    # a[i][j][l] = #{x in C_i : x^-1 z_l in C_j}
    M=np.zeros((k,k,k))
    for i in range(k):
        for x in C[i]:
            xi=inv(x)
            for l,z in enumerate(reps): M[i][cl[idx[mul(xi,z)]]][l]+=1
    rng=np.random.default_rng(1); T=sum(rng.normal()*M[i] for i in range(k))
    w,V=np.linalg.eig(T)   # vecteurs propres a droite : M_i omega = omega_i omega
    chars=[]
    sizes=np.array([len(c) for c in C])
    for j in range(k):
        v=V[:,j]                       # v_l proportionnel a |C_l| chi(g_l)/chi(1) ?
        # normalisation : omega_l = |C_l| chi_l / chi_1 ; on prend v/v[classe identite]
        e0=[l for l in range(k) if reps[l]==tuple(range(len(G[0])))][0]
        om=v/v[e0]
        chi_over=om/sizes             # chi_l/chi_1
        d2=N/np.sum(np.abs(chi_over)**2*sizes)
        d=np.sqrt(d2.real); chars.append(d*chi_over)
    return C,cl,idx,reps,np.array(chars)
def decompose(G, eps, name=''):
    C,cl,idx,reps,X=char_table(G); N=len(G); k=len(C)
    sizes=np.array([len(c) for c in C]); epsc=np.array([eps(r) for r in reps])
    sq=[cl[idx[mul(r,r)]] for r in reps]
    out=[]; used=set()
    for a in range(k):
        if a in used: continue
        chi=X[a]; d=int(round(chi[[l for l in range(k) if reps[l]==tuple(range(len(G[0])))][0]].real))
        tw=np.conj(chi)*epsc
        b=[bb for bb in range(k) if np.allclose(X[bb],tw,atol=1e-6)]
        if not b: raise RuntimeError('table incoherente')
        b=b[0]
        if b!=a:
            used|={a,b}; out.append(f'M{d}(C)')
        else:
            used.add(a)
            nu=sum(sizes[l]*epsc[l]*chi[sq[l]] for l in range(k)).real/N
            if round(nu)==1: out.append(f'M{d}(R)')
            elif round(nu)==-1: out.append(f'M{d//2}(H)' if d>2 else 'H')
            else: out.append(f'?nu={nu:.2f}')
    dimR=sum({'R':1,'C':2,'H':4}[s.split('(')[1][0]]*int(s[1:s.index('(')])**2 if s.startswith('M') else 4 for s in out)
    print(f'{name:28s} |G|={N:6d} classes={k:3d} : '+' + '.join(sorted(out, key=lambda s:(len(s),s)))+f'   (dim {dimR})', flush=True)
    return out
def sign(p):
    s=1; seen=set()
    for i in range(len(p)):
        if i in seen: continue
        j=i; c=0
        while j not in seen: seen.add(j); j=p[j]; c+=1
        s*= -1 if c%2==0 else 1
    return s
