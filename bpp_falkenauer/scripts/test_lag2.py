import math, random

def parse():
    with open("/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt") as f:
        lines=[l.rstrip("\n") for l in f]
    P=int(lines[0].strip()); idx=1
    name=lines[idx].strip(); idx+=1
    parts=lines[idx].split(); idx+=1
    C=int(parts[0]); n=int(parts[1]); best=int(parts[2])
    sizes=[int(lines[idx+i]) for i in range(n)]
    return name,C,n,best,sizes
name,C,n,best,sizes=parse(); C=150

def knap(profits):
    n=len(profits); dp=[0.0]*(C+1); keep=[[False]*(C+1) for _ in range(n)]
    for i in range(n):
        w=sizes[i]; p=profits[i]
        for c in range(C,w-1,-1):
            cand=dp[c-w]+p
            if cand>dp[c]+1e-12:
                dp[c]=cand; keep[i][c]=True
    bestc=max(range(C+1),key=lambda c:dp[c]); val=dp[bestc]
    items=[]; c=bestc
    for i in range(n-1,-1,-1):
        if keep[i][c]:
            items.append(i); c-=sizes[i]
    return val, set(items)

K=49
# try several strategies
import time
def run(niter=5000, steprule="dim", target=47.2):
    lam=[1.0]*n
    bestL=-1e9; bestlam=None
    for it in range(niter):
        v,S=knap(lam)
        L=sum(lam)+K*min(0.0,1.0-v)
        if L>bestL: bestL=L; bestlam=lam[:]
        if v>1.0:
            g=[1.0-(K if i in S else 0.0) for i in range(n)]
        else:
            g=[1.0]*n
        gnorm=math.sqrt(sum(gi*gi for gi in g))
        if gnorm<1e-12: break
        if steprule=="dim":
            step=50.0/(it+1)
        elif steprule=="polyak":
            step=max(0.0,(target-L))/max(gnorm*gnorm,1e-12)
        elif steprule=="polyak2":
            step=max(0.0,(bestL+1.0-L))/max(gnorm*gnorm,1e-12)
        for i in range(n): lam[i]+=step*g[i]
    return bestL
for rule in ["dim","polyak","polyak2"]:
    t0=time.time()
    bl=run(5000,rule)
    print(rule,"bestL",round(bl,5),"time",round(time.time()-t0,1))
