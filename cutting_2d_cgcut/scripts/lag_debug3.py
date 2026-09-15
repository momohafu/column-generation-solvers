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
def best_pattern(c):
    @lru_cache(maxsize=None)
    def F(w,h):
        best_val=0.0; best_counts=(0,)*m
        for i,(a,b,qq,vv) in enumerate(items):
            if a<=w and b<=h:
                if c[i]>best_val:
                    cnt=[0]*m; cnt[i]=1; best_val=c[i]; best_counts=tuple(cnt)
        for x in range(1,w):
            v1,c1=F(x,h); v2,c2=F(w-x,h); val=v1+v2
            if val>best_val: best_val=val; best_counts=tuple(c1[j]+c2[j] for j in range(m))
        for y in range(1,h):
            v1,c1=F(w,y); v2,c2=F(w,h-y); val=v1+v2
            if val>best_val: best_val=val; best_counts=tuple(c1[j]+c2[j] for j in range(m))
        return best_val,best_counts
    return F(W,H)
LB=244.0
for cap in [0.5, 2.0, 10.0]:
    pi=[0.0]*m; best_UB=1e18
    for t in range(5000):
        c=[v[i]-pi[i] for i in range(m)]
        val,counts=best_pattern(c)
        UB=sum(pi[i]*q[i] for i in range(m))+val
        if UB<best_UB: best_UB=UB
        g=[q[i]-counts[i] for i in range(m)]
        gnorm2=sum(gg*gg for gg in g)
        if gnorm2<1e-12: break
        step=(UB-LB)/gnorm2
        step=max(0.0,min(step,cap))
        pi=[max(0.0,pi[i]-step*g[i]) for i in range(m)]
    print("cap",cap,"best_UB",best_UB,"gap%",(best_UB-LB)/LB*100)
