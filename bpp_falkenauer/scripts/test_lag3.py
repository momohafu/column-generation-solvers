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
name,C,n,best,sizes=parse(); C=150; n=120

def knap(profits):
    dp=[0.0]*(C+1); keep=[[False]*(C+1) for _ in range(n)]
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
target=47.27
# volume-like: lambda_avg evaluated; damping on lambda
lam=[1.0]*n
lam_avg=[1.0]*n
bestL=-1e9
t0=time.time()
for it in range(20000):
    v,S=knap(lam_avg)
    L=sum(lam_avg)+K*min(0.0,1.0-v)
    if L>bestL: bestL=L
    # subgradient at current lambda (not avg)
    vc,Sc=knap(lam)
    if vc>1.0:
        g=[1.0-(K if i in Sc else 0.0) for i in range(n)]
    else:
        g=[1.0]*n
    gnorm2=sum(gi*gi for gi in g)
    if gnorm2<1e-12: break
    step=max(0.0,(target-L))/gnorm2
    step*=0.5  # damping
    for i in range(n):
        lam[i]+=step*g[i]
        lam_avg[i]=0.9*lam_avg[i]+0.1*lam[i]
    if it%2000==0:
        print("it",it,"L",round(L,5),"best",round(bestL,5),"v",round(v,4))
print("bestL",round(bestL,5),"time",round(time.time()-t0,1))
