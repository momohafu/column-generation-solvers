import math, time

def parse():
    with open("/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt") as f:
        lines=[l.rstrip("\n") for l in f]
    P=int(lines[0].strip()); idx=1
    name=lines[idx].strip(); idx+=1
    parts=lines[idx].split(); idx+=1
    C=int(parts[0]); n=int(parts[1]); best=int(parts[2])
    sizes=[int(lines[idx+i]) for i in range(n)]
    return name,C,n,best,sizes
name,C,n,best,sizes=parse()
C=150

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

def ffd():
    order=sorted(range(n),key=lambda i:-sizes[i])
    bins=[]; asg={}
    for i in order:
        s=sizes[i]; pl=False
        for b in range(len(bins)):
            if bins[b]>=s: bins[b]-=s; asg[i]=b; pl=True; break
        if not pl: asg[i]=len(bins); bins.append(C-s)
    return len(bins),asg

ub,_=ffd()
print("FFD UB", ub)
for K in [48, 49]:
    lam=[1.0]*n
    bestL=-1e9
    Lhist=[]
    for it in range(2000):
        v,S = knap(lam)
        L = sum(lam) + K*min(0.0, 1.0-v)
        if L>bestL: bestL=L
        # subgradient
        g=[0.0]*n
        if v>1.0:
            for i in range(n):
                g[i] = 1.0 - (K if i in S else 0.0)
        else:
            g=[1.0]*n
        gnorm=math.sqrt(sum(gi*gi for gi in g))
        if gnorm<1e-9: break
        step = 0.5*(ub-L)/gnorm if False else 1.0/(it+1)*100
        # polyak step with target ub
        step = (ub - L)/(gnorm*gnorm) if gnorm>1e-12 else 0.0
        for i in range(n):
            lam[i]+=step*g[i]
        if it%200==0:
            print("K",K,"it",it,"L",round(L,4),"best",round(bestL,4),"knap",round(v,4))
    print("K",K,"final bestL",round(bestL,4),"ub",ub,"gap%",round((ub-bestL)/ub*100,3))
