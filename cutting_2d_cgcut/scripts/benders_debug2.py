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
P_all=sorted(patterns(W,H)); val={p:sum(p[i]*v[i] for i in range(m)) for p in P_all}
cand=set([(0,)*m])
for i in range(m):
    t=[0]*m; t[i]=1; cand.add(tuple(t))
top=sorted(P_all, key=lambda p:-val[p])[:40]
cand.update(top); cand=sorted(cand); C=cand; n=len(C)
def solve_sub(ybar):
    model=mathopt.Model(name="sub")
    x=[model.add_variable(lb=0.0,ub=float("inf"),is_integer=False,name=f"x{k}") for k in range(n)]
    con_item=[model.add_linear_constraint(sum(C[k][i]*x[k] for k in range(n))<=q[i]) for i in range(m)]
    con_sheet=model.add_linear_constraint(sum(x[k] for k in range(n))<=1)
    con_ub=[model.add_linear_constraint(x[k]<=ybar[k]) for k in range(n)]
    model.maximize(sum(val[C[k]]*x[k] for k in range(n)))
    res=mathopt.solve(model, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60),enable_output=False))
    dv=res.dual_values(); pi=[dv[con_item[i]] for i in range(m)]; mu=dv[con_sheet]; gamma=[dv[con_ub[k]] for k in range(n)]
    return res.objective_value(), pi, mu, gamma
def solve_master(cuts, bigM=1e6):
    model=mathopt.Model(name="master")
    y=[model.add_variable(lb=0.0,ub=1.0,is_integer=True,name=f"y{k}") for k in range(n)]
    theta=model.add_variable(lb=0.0,ub=bigM,is_integer=False,name="theta")
    model.add_linear_constraint(sum(y[k] for k in range(n))<=1)
    model.add_linear_constraint(theta <= sum(val[C[k]]*y[k] for k in range(n)) + bigM*(1 - sum(y[k] for k in range(n))))
    for (rhs_pi, rhs_mu, rhs_gamma) in cuts:
        model.add_linear_constraint(theta <= rhs_pi + rhs_mu + sum(rhs_gamma[k]*y[k] for k in range(n)))
    model.maximize(theta)
    res=mathopt.solve(model, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60),enable_output=False))
    yv=[round(res.variable_values()[y[k]]) for k in range(n)]
    return yv, res.objective_value()
cuts=[]; best_LB=0; best_y=None
for it in range(100):
    ystar, ub=solve_master(cuts)
    LB, pi, mu, gamma=solve_sub(ystar)
    if LB>best_LB: best_LB=LB; best_y=ystar
    print(f"it={it+1} UB={ub:.6f} LB={LB:.6f} gap={ub-LB:.3e} ncuts={len(cuts)} sel={[k for k in range(n) if ystar[k]==1]}")
    if ub-LB<=1e-6 or it>=60: break
    cuts.append((sum(q[i]*pi[i] for i in range(m)), mu, gamma))
print("best_LB",best_LB,"pattern",[C[k] for k in range(n) if best_y and best_y[k]==1])
