import sys, time, datetime
sys.setrecursionlimit(1000000)
from functools import lru_cache
from ortools.math_opt.python import mathopt
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
P=sorted(patterns(W,H)); val={p:sum(p[i]*v[i] for i in range(m)) for p in P}
def solve_rmp(cols):
    model=mathopt.Model(name="rmp")
    x=[model.add_variable(lb=0.0,ub=float("inf"),is_integer=False,name=f"x{k}") for k in range(len(cols))]
    con_items=[model.add_linear_constraint(sum(cols[k][i]*x[k] for k in range(len(cols)))<=q[i]) for i in range(m)]
    con_sheet=model.add_linear_constraint(sum(x[k] for k in range(len(cols)))<=1)
    model.maximize(sum(val[cols[k]]*x[k] for k in range(len(cols))))
    res=mathopt.solve(model, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60),enable_output=False))
    return model,x,con_items,con_sheet,res
initial=[(0,)*m]+[tuple(1 if i==j else 0 for j in range(m)) for i in range(m)]
cols=list(initial)
for it in range(1000):
    model,x,con_items,con_sheet,res=solve_rmp(cols)
    lp=res.objective_value(); dv=res.dual_values()
    pi=[dv[con_items[i]] for i in range(m)]; mu=dv[con_sheet]
    best_rc=-1e18; best_p=None
    for p in P:
        rc=val[p]-sum(pi[i]*p[i] for i in range(m))-mu
        if rc>best_rc: best_rc=rc; best_p=p
    print(f"it={it+1} lp={lp:.6f} rc={best_rc:.3e} ncols={len(cols)} added={best_p if best_rc>1e-7 else None}")
    if best_rc<=1e-7: break
    cols.append(best_p)
print("generated cols:")
for c in cols:
    print(c, "val", val[c])
# IP over generated
model3=mathopt.Model(name="gip")
xg=[model3.add_variable(lb=0.0,ub=1.0,is_integer=True,name=f"x{k}") for k in range(len(cols))]
for i in range(m):
    model3.add_linear_constraint(sum(cols[k][i]*xg[k] for k in range(len(cols)))<=q[i])
model3.add_linear_constraint(sum(xg)<=1)
model3.maximize(sum(val[cols[k]]*xg[k] for k in range(len(cols))))
res3=mathopt.solve(model3, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60),enable_output=False))
print("gen IP",res3.termination.reason, res3.objective_value())
print("sol",[cols[k] for k in range(len(cols)) if res3.variable_values()[xg[k]]>0.5])
