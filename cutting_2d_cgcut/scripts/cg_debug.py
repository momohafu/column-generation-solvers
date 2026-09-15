import sys, time, datetime
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
path="/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt"
m,W,H,items=parse(path)
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
t=time.time(); P=sorted(patterns(W,H)); print("pool",len(P),"enum_time",round(time.time()-t,2))
val={p: sum(p[i]*v[i] for i in range(m)) for p in P}

def solve_rmp(cols):
    model=mathopt.Model(name="rmp")
    x=[model.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in range(len(cols))]
    con_items=[]
    for i in range(m):
        con_items.append(model.add_linear_constraint(sum(cols[k][i]*x[k] for k in range(len(cols))) <= q[i]))
    con_sheet=model.add_linear_constraint(sum(x[k] for k in range(len(cols))) <= 1)
    model.maximize(sum(val[cols[k]]*x[k] for k in range(len(cols))))
    res=mathopt.solve(model, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
    return model,x,con_items,con_sheet,res

def run_cg(max_iter=1000, tol=1e-7):
    initial=[(0,)*m] + [tuple(1 if i==j else 0 for j in range(m)) for i in range(m)]
    cols=list(initial)
    hist=[]
    t0=time.time()
    for it in range(max_iter):
        model,x,con_items,con_sheet,res=solve_rmp(cols)
        lp=res.objective_value()
        dv=res.dual_values()
        pi=[dv[con_items[i]] for i in range(m)]
        mu=dv[con_sheet]
        best_rc=-1e18; best_p=None
        for p in P:
            rc=val[p] - sum(pi[i]*p[i] for i in range(m)) - mu
            if rc>best_rc:
                best_rc=rc; best_p=p
        hist.append((it+1, lp, best_rc, len(cols)))
        if it<3 or it%10==0 or best_rc<=tol:
            print(f"it={it+1} lp={lp:.6f} best_rc={best_rc:.3e} ncols={len(cols)}")
        if best_rc<=tol:
            print("CONVERGED")
            break
        if best_p not in cols:
            cols.append(best_p)
        else:
            print("no new col, rc",best_rc); break
    dt=time.time()-t0
    return cols,hist,dt

cols,hist,dt=run_cg()
print("cg time",round(dt,2),"ncols",len(cols))
# full pool LP
model=mathopt.Model(name="full_lp")
x=[model.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in range(len(P))]
for i in range(m):
    model.add_linear_constraint(sum(P[k][i]*x[k] for k in range(len(P))) <= q[i])
model.add_linear_constraint(sum(x) <= 1)
model.maximize(sum(val[P[k]]*x[k] for k in range(len(P))))
res=mathopt.solve(model, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
print("full LP obj",res.objective_value(),res.termination.reason)
# IP over full pool
model2=mathopt.Model(name="full_ip")
xi=[model2.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{k}") for k in range(len(P))]
for i in range(m):
    model2.add_linear_constraint(sum(P[k][i]*xi[k] for k in range(len(P))) <= q[i])
model2.add_linear_constraint(sum(xi) <= 1)
model2.maximize(sum(val[P[k]]*xi[k] for k in range(len(P))))
res2=mathopt.solve(model2, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
print("full IP obj",res2.objective_value(),res2.termination.reason)
sol=[k for k in range(len(P)) if res2.variable_values()[xi[k]]>0.5]
print("full IP pattern",[P[k] for k in sol])
# IP over generated cols
model3=mathopt.Model(name="gen_ip")
xg=[model3.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{k}") for k in range(len(cols))]
for i in range(m):
    model3.add_linear_constraint(sum(cols[k][i]*xg[k] for k in range(len(cols))) <= q[i])
model3.add_linear_constraint(sum(xg) <= 1)
model3.maximize(sum(val[cols[k]]*xg[k] for k in range(len(cols))))
res3=mathopt.solve(model3, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
print("gen IP obj",res3.objective_value(),res3.termination.reason)
sol3=[k for k in range(len(cols)) if res3.variable_values()[xg[k]]>0.5]
print("gen IP pattern",[cols[k] for k in sol3])
