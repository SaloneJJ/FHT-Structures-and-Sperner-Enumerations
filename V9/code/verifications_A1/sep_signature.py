"""Separation of the 260 FHTs with at most five vertices (Article I, Section "How much do these invariants distinguish?").
Reads the tables of Appendices C, D, E (file A1_appendix_tables.tex) and counts the distinct values of the signatures.
The last lines ("noroot-dist 5, 6, 7") are the columns (v), (vi), (vii) of the article: purely numerical signatures,
root host line distinguished, no isomorphism type of an FHT or of a host.  Earlier lines are the variants discussed in
the text, in particular the variant keeping the flat types F_f of the hosts (230/259 classes at (v)), which identifies a flat
FHT by its own type.  Usage: python3 sep_signature.py A1_appendix_tables.tex
"""
import re,sys,collections
s=open(sys.argv[1],encoding='utf-8').read()
def sect(lab):
    i=s.index('\\label{%s}'%lab); j=s.find('\\section',i+10); return s[i:j if j>0 else len(s)]
C=sect('sec:FHT_list'); D=sect('sec:annex'); E=sect('sec:annexB')
info={}
for m in re.finditer(r'^\$\\F_\{(\d+)\}\$ & ([^&]*) & (\d+) & (\d+) & (\d+) \\\\',C,re.M):
    info[int(m.group(1))]=dict(v=int(m.group(3)),N=int(m.group(4)),E=int(m.group(5)))
cells=lambda line:[c.strip() for c in line.rstrip().rstrip('\\').split('&')]
glob={};hosts=collections.defaultdict(list);cur=None
for line in D.split('\n'):
    m=re.match(r'\\textbf\{\$\\F_\{(\d+)\}\$\}',line)
    if m: cur=int(m.group(1)); glob[cur]=cells(line)[1:]; continue
    if line.startswith('\\quad') and cur:
        c=cells(line); lab=c[0]
        mu=re.search(r'\{\\times\}(\d+)',lab); mu=int(mu.group(1)) if mu else 1
        hosts[cur].append((lab,mu,c[1:]))
eo={}
for m in re.finditer(r'^\$\\F_\{(\d+)\}\$ & .*$',E,re.M):
    c=cells(m.group(0)); eo[int(m.group(1))]=c[-1]
print(len(info),len(glob),len(eo))
def sig(n,lev,ff=True,mult=True,labels=False):
    g=glob[n]; i=info[n]
    t=[(i['v'],i['N'],i['E'])]
    if lev>=2: t.append(tuple(g[0:5]))
    if lev>=3: t.append(g[5])
    if lev>=4: t.append(tuple(g[6:15]))
    if lev>=5:
        L=[]
        for lab,mu,c in hosts[n]:
            c=list(c)
            if not ff: c[5]=re.sub(r'\\F_\{\d+\}','',c[5])
            key=tuple(c)+((lab,) if labels else ())
            L+= [key]*(mu if mult else 1)
        t.append(tuple(sorted(L)))
    if lev>=6: t.append(eo.get(n,'-'))
    return tuple(t)
def table(**kw):
    out=[]
    for rng in [[3],[4],[5],[1,2,3,4,5]]:
        ns=[n for n in info if info[n]['v'] in rng]
        out.append([len(ns)]+[len(set(sig(n,l,**kw) for n in ns)) for l in range(1,7)])
    return out
for kw in [dict(),dict(ff=False),dict(mult=False),dict(ff=False,mult=False),dict(labels=True)]:
    print(kw,table(**kw))
print('--- flat root line added (F_f = own type)')
for n in info:
    if info[n]['N']==1 and not hosts[n]:
        g=glob[n]; hosts[n].append(('root',1,g[0:5]+['%d|%d|\\F_{%d}'%(info[n]['v'],info[n]['E'],n)]+g[6:15]))
for kw in [dict(),dict(ff=False)]:
    print(kw,table(**kw))
for lev in (5,6):
    d=collections.defaultdict(list)
    for n in info: d[sig(n,lev)].append(n)
    print(lev,[v for v in d.values() if len(v)>1])
def sig2(n,lev,ff=True,dist=False,flatline=True):
    g=glob[n]; i=info[n]
    t=[(i['v'],i['N'],i['E'])]
    if lev>=2: t.append(tuple(g[0:5]))
    if lev>=3: t.append(g[5])
    if lev>=4: t.append(tuple(g[6:15]))
    if lev>=5:
        root=None;L=[]
        for lab,mu,c in hosts[n]:
            if lab=='root' and not flatline and i['N']==1: continue
            c=list(c)
            if not ff: c[5]=re.sub(r'\\F_\{\d+\}','',c[5])
            if lab in ('root','\\quad\\textit{root}') or 'textit{root}' in lab: root=tuple(c); continue
            dd=re.search(r'\^\{\((\d+)\)\}',lab); dd=int(dd.group(1)) if dd else 1
            L+=[tuple(c)+((dd,) if dist else ())]*mu
        t.append((root,tuple(sorted(L))))
    if lev>=6: t.append(eo.get(n,'-'))
    return tuple(t)
for kw in [dict(),dict(dist=True),dict(ff=False),dict(flatline=False),dict(ff=False,flatline=False),dict(ff=False,dist=True)]:
    out=[]
    for rng in [[3],[4],[5],[1,2,3,4,5]]:
        ns=[n for n in info if info[n]['v'] in rng]
        out.append([len(set(sig2(n,l,**kw) for n in ns)) for l in (5,6)])
    print(kw,out)
for lev in (5,6):
    d=collections.defaultdict(list)
    for n in info: d[sig2(n,lev,ff=False)].append(n)
    print('noff',lev,[(v,[info[k]['N'] for k in v]) for v in d.values() if len(v)>1])
Erow={}
for m in re.finditer(r'^\$\\F_\{(\d+)\}\$ & .*$',E,re.M):
    Erow[int(m.group(1))]=cells(m.group(0))[1:]
def sig3(n,variant):
    base=sig2(n,6,ff=False)
    r=Erow.get(n,['-']*7)
    if variant=='num': extra=(r[1],r[2],r[3],r[4])
    if variant=='num+prof': extra=(re.sub(r'\\F_\{\d+\}','',r[0]),r[1],r[2],r[3],r[4],r[5])
    if variant=='full': extra=tuple(r)
    return base+(extra,)
for var in ['num','num+prof','full']:
    out=[]
    for rng in [[3],[4],[5],[1,2,3,4,5]]:
        ns=[n for n in info if info[n]['v'] in rng]
        out.append(len(set(sig3(n,var) for n in ns)))
    d=collections.defaultdict(list)
    for n in info: d[sig3(n,var)].append(n)
    print(var,out,[v for v in d.values() if len(v)>1])
def sig4(n): return sig2(n,6,ff=False)+(Erow.get(n,['-']*7)[1],)
out=[]
for rng in [[3],[4],[5],[1,2,3,4,5]]:
    ns=[n for n in info if info[n]['v'] in rng]; out.append(len(set(sig4(n) for n in ns)))
print('sizes only',out)
# does root-distinguish matter without ff?
def sig5(n,lev):
    t=list(sig(n,min(lev,4)))
    if lev>=5:
        L=[]
        for lab,mu,c in hosts[n]:
            c=list(c); c[5]=re.sub(r'\\F_\{\d+\}','',c[5]); L+=[tuple(c)]*mu
        t.append(tuple(sorted(L)))
    if lev>=6: t.append(eo.get(n,'-'))
    if lev>=7: t.append(Erow.get(n,['-']*7)[1])
    return tuple(t)
hosts2={n:[h for h in hosts[n] if not (info[n]['N']==1)] for n in hosts}
for lev in (5,6,7):
    out=[]
    for rng in [[3],[4],[5],[1,2,3,4,5]]:
        ns=[n for n in info if info[n]['v'] in rng]; out.append(len(set(sig5(n,lev) for n in ns)))
    print('noroot-dist',lev,out)
roots={}
for n in info:
    if n==1: continue
    if info[n]['N']==1: roots[n]=glob[n][6:15]
    else:
        r=[c for lab,mu,c in hosts[n] if 'root' in lab]; roots[n]=r[0][6:15]
print('nroots',len(roots), collections.Counter(int(v[7]) - 2*int(v[4]) if False else v[7] for v in roots.values()))
import itertools
A=B=0
for (a,x),(b,y) in itertools.combinations(roots.items(),2):
    if x[0:2]==y[0:2] and x[3:5]==y[3:5] and x[6]!=y[6]: A+=1
    if x[0:6]==y[0:6] and x[8]!=y[8]: B+=1
print('pairs',A,B)
# dia(IG) vs dia(LG)
print(collections.Counter((int(v[7])-2*int(v[4])) for v in roots.values()))
print('diaIG-2max',collections.Counter(int(v[7])-2*max(int(v[1]),int(v[4])) for v in roots.values()))
# distinct root hosts (as flat types) 
