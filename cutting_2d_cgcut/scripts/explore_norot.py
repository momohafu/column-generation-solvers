import sys, time, pickle, datetime
sys.setrecursionlimit(1000000)
from functools import lru_cache
from ortools.math_opt.python import mathopt
def parse(path):
    lines=[l.strip() for l in open(path).read().strip().splitlines() if l.strip()]
    m=int(lines[0]); W,H=map(int,lines[1].split())
    items=[]
    for l in lines[2:2+m]:
        a,b,q,v=map(int,l.split())
        items.append((a,b,q,v))
    return m,W,H,items
m,W,H,items=parse("/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt")
q=tuple(it[2] for it in items); v=tuple(it[3] for it in items)
orients=[((a,b),) for (a,b,q,v) in items]
@lru_cache(maxsize=None)
def patterns(w,h):
    res={(0,)*m}
    for i in range(m):
        for (dw,dh) in orients[i]:
            if dw<=w and dh<=h:
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
P=sorted(patterns(W,H)); print("norot patterns",len(P))
model=mathopt.Model(name="ip")
xv=[model.add_variable(lb=0.0,ub=1.0,is_integer=True,name=f"x{k}") for k in range(len(P))]
for i in range(m):
    model.add_linear_constraint(sum(P[k][i]*xv[k] for k in range(len(P)))<=q[i])
model.add_linear_constraint(sum(xv)<=1)
model.maximize(sum(sum(P[k][j]*v[j] for j in range(m))*xv[k] for k in range(len(P))))
res=mathopt.solve(model, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120.0),enable_output=False))
print("norot obj",res.objective_value(),res.termination.reason)
sol=[k for k in range(len(P)) if res.variable_values()[xv[k]]>0.5]
print("used",[P[k] for k in sol])
