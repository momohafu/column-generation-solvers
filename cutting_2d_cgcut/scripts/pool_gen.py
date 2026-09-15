import sys, time, pickle, datetime
sys.setrecursionlimit(1000000)
from functools import lru_cache

def parse(path):
    lines=[l.strip() for l in open(path).read().strip().splitlines() if l.strip()]
    m=int(lines[0]); W,H=map(int,lines[1].split())
    items=[]
    for l in lines[2:2+m]:
        a,b,q,v=map(int,l.split())
        items.append((a,b,q,v))
    return m,W,H,items

def build(path):
    m,W,H,items=parse(path)
    q=tuple(it[2] for it in items)
    v=tuple(it[3] for it in items)
    orients=[]
    for i,(a,b,qq,vv) in enumerate(items):
        orients.append(tuple(sorted(set([(a,b),(b,a)]))))
    @lru_cache(maxsize=None)
    def patterns(w,h):
        res=set()
        res.add((0,)*m)
        for i in range(m):
            for (dw,dh) in orients[i]:
                if dw<=w and dh<=h:
                    t=[0]*m; t[i]=1
                    res.add(tuple(t))
        for x in range(1,w):
            P1=patterns(x,h); P2=patterns(w-x,h)
            for p1 in P1:
                for p2 in P2:
                    t=tuple(a+b for a,b in zip(p1,p2))
                    if all(t[j]<=q[j] for j in range(m)):
                        res.add(t)
        for y in range(1,h):
            P1=patterns(w,y); P2=patterns(w,h-y)
            for p1 in P1:
                for p2 in P2:
                    t=tuple(a+b for a,b in zip(p1,p2))
                    if all(t[j]<=q[j] for j in range(m)):
                        res.add(t)
        return res
    t=time.time()
    P=sorted(patterns(W,H))
    dt=time.time()-t
    # add zero pattern if not present (present)
    print("patterns",len(P),"time",dt)
    return m,W,H,items,P

if __name__=="__main__":
    path=sys.argv[1] if len(sys.argv)>1 else "/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt"
    out=sys.argv[2] if len(sys.argv)>2 else "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/scripts/pool.pkl"
    m,W,H,items,P=build(path)
    with open(out,"wb") as f:
        pickle.dump({"m":m,"W":W,"H":H,"items":items,"patterns":P}, f)
    print("saved",out)
