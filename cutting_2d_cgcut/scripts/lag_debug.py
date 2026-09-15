import sys, time
sys.setrecursionlimit(1000000)
from functools import lru_cache
def parse(path):
    lines=[l.strip() for l in open(path).read().strip().splitlines() if l.strip()]
    m=int(lines[0]); W,H=map(int,lines[1].split())
    items=[]
    for l in lines[2:2+m]:
        a,b,q,v=map(int,l.split()); items.append((a,b,q,v))
    return m,W,H,items
m,W,H,items=parse("/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt")
q=[it[2] for it in items]; v=[it[3] for it in items]
# unbounded guillotine DP with argmax counts, profits c
def best_pattern(c):
    @lru_cache(maxsize=None)
    def F(w,h):
        best_val=0.0; best_counts=(0,)*m
        for i,(a,b,qq,vv) in enumerate(items):
            if a<=w and b<=h:
                val=c[i]
                if val>best_val:
                    cnt=[0]*m; cnt[i]=1
                    best_val=val; best_counts=tuple(cnt)
        for x in range(1,w):
            v1,c1=F(x,h); v2,c2=F(w-x,h)
            val=v1+v2
            if val>best_val:
                best_val=val; best_counts=tuple(c1[j]+c2[j] for j in range(m))
        for y in range(1,h):
            v1,c1=F(w,y); v2,c2=F(w,h-y)
            val=v1+v2
            if val>best_val:
                best_val=val; best_counts=tuple(c1[j]+c2[j] for j in range(m))
        return best_val,best_counts
    return F(W,H)
# exact LB via capped pool max
@lru_cache(maxsize=None)
def capped_patterns(w,h):
    res={(0,)*m}
    for i,(a,b,qq,vv) in enumerate(items):
        if a<=w and b<=h:
            t=[0]*m; t[i]=1; res.add(tuple(t))
    for x in range(1,w):
        P1=capped_patterns(x,h); P2=capped_patterns(w-x,h)
        for p1 in P1:
            for p2 in P2:
                t=tuple(a+b for a,b in zip(p1,p2))
                if all(t[j]<=q[j] for j in range(m)): res.add(t)
    for y in range(1,h):
        P1=capped_patterns(w,y); P2=capped_patterns(w,h-y)
        for p1 in P1:
            for p2 in P2:
                t=tuple(a+b for a,b in zip(p1,p2))
                if all(t[j]<=q[j] for j in range(m)): res.add(t)
    return res
P=capped_patterns(W,H); LB=max(sum(p[i]*v[i] for i in range(m)) for p in P)
print("capped pool",len(P),"exact LB",LB)
pi=[0.0]*m
best_UB=1e18; best_pi=None; best_UB_iter=0
T=500
UB_hist=[]
LB_best=0
for t in range(T):
    c=[v[i]-pi[i] for i in range(m)]
    val,counts=best_pattern(c)
    UB=sum(pi[i]*q[i] for i in range(m)) + val
    if UB<best_UB:
        best_UB=UB; best_pi=list(pi); best_UB_iter=t
    # repair: clip counts
    clip=[min(counts[i],q[i]) for i in range(m)]
    lb=sum(v[i]*clip[i] for i in range(m))
    if lb>LB_best: LB_best=lb
    g=[q[i]-counts[i] for i in range(m)]
    gnorm2=sum(gg*gg for gg in g)
    if gnorm2<1e-12: break
    # Polyak-ish step; use LB (exact 244) as target to show convergence? use current LB_best
    step=(best_UB - LB_best + 1.0)/gnorm2
    step=min(step, 2.0)
    step*=0.5
    pi=[max(0.0, pi[i]+step*g[i]) for i in range(m)]
    if t%50==0 or t==T-1:
        print(f"t={t} UB={UB:.6f} bestUB={best_UB:.6f} LB_best={LB_best} step={step:.4f} counts={counts} pi={[round(x,3) for x in pi]}")
print("final best_UB",best_UB,"best_LB",LB_best,"gap",(best_UB-LB)/LB*100)
print("best_pi",[round(x,4) for x in best_pi])
