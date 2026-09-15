import sys, time
sys.setrecursionlimit(1000000)
from functools import lru_cache
from ortools.sat.python import cp_model
def parse(path):
    lines=[l.strip() for l in open(path).read().strip().splitlines() if l.strip()]
    m=int(lines[0]); W,H=map(int,lines[1].split())
    items=[]
    for l in lines[2:2+m]:
        a,b,q,v=map(int,l.split()); items.append((a,b,q,v))
    return m,W,H,items
m,W,H,items=parse("/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt")
q=[it[2] for it in items]; v=[it[3] for it in items]
@lru_cache(maxsize=None)
def patterns(w,h):
    res={(0,)*m}
    for i,(a,b,qq,vv) in enumerate(items):
        if a<=w and b<=h:
            t=[0]*m; t[i]=1; res.add(tuple(t))
    for x in range(1,w):
        P1=patterns(x,h); P2=patterns(w-x,h)
        for p1 in P1:
            for p2 in P2:
                t=tuple(a+b for a,b in zip(p1,p2))
                if all(t[j]<=q[j] for j in range(m)): res.add(t)
    for y in range(1,h):
        P1=patterns(w,y); P2=patterns(w,h-y)
        for p1 in P1:
            for p2 in P2:
                t=tuple(a+b for a,b in zip(p1,p2))
                if all(t[j]<=q[j] for j in range(m)): res.add(t)
    return res
P=patterns(W,H); Pset=set(P); print("pool",len(P))
def is_feasible(n):
    return tuple(n) in Pset
# master with cuts
model=cp_model.CpModel()
n=[model.NewIntVar(0,q[i],f"n{i}") for i in range(m)]
model.Maximize(sum(v[i]*n[i] for i in range(m)))
solver=cp_model.CpSolver()
solver.parameters.max_time_in_seconds=120
solver.parameters.num_search_workers=8
cuts=[]
t0=time.time()
for it in range(100):
    status=solver.Solve(model)
    ns=[solver.Value(n[i]) for i in range(m)]
    obj=sum(v[i]*ns[i] for i in range(m))
    feas=is_feasible(ns)
    print(f"it={it+1} obj={obj} ns={ns} feasible={feas} ncalls={len(cuts)} status={solver.StatusName(status)}")
    if feas:
        print("FEASIBLE found",obj); break
    if status!=cp_model.OPTIMAL:
        print("master not optimal"); break
    # add no-good: sum over i with ns[i]>0 of delta_i >=1, delta_i=1 iff n_i <= ns[i]-1
    ds=[]; constr=[]
    for i in range(m):
        if ns[i]>0:
            d=model.NewBoolVar(f"d_{it}_{i}")
            model.Add(n[i] <= (ns[i]-1) + q[i]*(1-d))
            model.Add(n[i] >= ns[i]*(1-d))
            ds.append(d)
    model.Add(sum(ds)>=1)
    cuts.append(tuple(ns))
print("time",round(time.time()-t0,2),"cuts",len(cuts))
