import sys, time, itertools
from enum_pieces import search
SZ={'S4':1,'A4':2,'D4':3,'S3':4,'C4':6,"V4'":6,'C3':8}
order=['S4','A4','D4','S3','C4',"V4'",'C3']
def multisets(maxn):
    out=[]
    def rec(i, cur, tot):
        if i==len(order):
            if cur: out.append(list(cur)); 
            return
        t=order[i]; k=0
        while tot+k*SZ[t]<=maxn:
            rec(i+1, cur+[t]*k, tot+k*SZ[t]); k+=1
    rec(0,[],0); return out
which=sys.argv[1].split(','); nmin,nmax=int(sys.argv[2]),int(sys.argv[3]); tl=float(sys.argv[4])
t0=time.time()
for W in sorted(multisets(nmax), key=lambda w:sum(SZ[t] for t in w)):
    n=sum(SZ[t] for t in W)
    if n<nmin: continue
    if time.time()-t0>165: print('arret (temps global)'); break
    import enum_pieces
    st,c,no,f=search(W, which, limit_time=tl)
    if f or st!='complet': print(n, '+'.join(W), no, c, st, [(w,len(s)) for w,s in f][:3], flush=True)
print('fin', round(time.time()-t0))
