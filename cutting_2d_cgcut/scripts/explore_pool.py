import sys, time
sys.setrecursionlimit(100000)
def parse(path):
    lines=[l.strip() for l in open(path).read().strip().splitlines() if l.strip()]
    m=int(lines[0]); W,H=map(int,lines[1].split())
    items=[]
    for l in lines[2:2+m]:
        a,b,q,v=map(int,l.split())
        items.append((a,b,q,v))
    return m,W,H,items
m,W,H,items=parse("/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt")
q=[it[2] for it in items]
v=[it[3] for it in items]
orient=[]
for i,(a,b,qq,vv) in enumerate(items):
    orient.append(sorted(set([(a,b),(b,a)])))
print("orient",orient)

from functools import lru_cache
# enumerate all guillotine patterns within q limits
@lru_cache(maxsize=None)
def patterns(w,h):
    res=set()
    res.add((0,)*m)
    for i in range(m):
        for (dw,dh) in orient[i]:
            if dw<=w and dh<=h:
                t=[0]*m; t[i]=1
                res.add(tuple(t))
    for x in range(1,w):
        P1=patterns(x,h); P2=patterns(w-x,h)
        for p1 in P1:
            for p2 in P2:
                t=tuple(min(q[j], p1[j]+p2[j]) for j in range(m))
                res.add(t)
    for y in range(1,h):
        P1=patterns(w,y); P2=patterns(w,h-y)
        for p1 in P1:
            for p2 in P2:
                t=tuple(min(q[j], p1[j]+p2[j]) for j in range(m))
                res.add(t)
    return res
t=time.time()
P=patterns(W,H)
print("num patterns",len(P),"time",time.time()-t)
# filter dominated: keep only count vectors within q
# solve ILP over P
from ortools.math_opt.python import mathopt
model=mathopt.Model(name="pool_ip")
xv=[model.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{k}") for k in range(len(P))]
Pl=list(P)
for i in range(m):
    model.add_linear_constraint(sum(Pl[k][i]*xv[k] for k in range(len(P))) <= q[i])
model.add_linear_constraint(sum(xv) <= 1)
model.maximize(sum(sum(Pl[k][j]*v[j] for j in range(m))*xv[k] for k in range(len(P))))
res=mathopt.solve(model, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=120.0, enable_output=False))
print("IP status",res.termination.reason,"obj",res.objective_value())
sol=[k for k in range(len(P)) if res.variable_values()[xv[k]]>0.5]
print("used patterns",[Pl[k] for k in sol])
