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
def solve_rmp(cols):
    model=mathopt.Model(name="rmp")
    x=[model.add_variable(lb=0.0,ub=float("inf"),is_integer=False,name=f"x{k}") for k in range(len(cols))]
    con_item=[model.add_linear_constraint(sum(cols[k][i]*x[k] for k in range(len(cols)))<=q[i]) for i in range(m)]
    con_sheet=model.add_linear_constraint(sum(x[k] for k in range(len(cols)))<=1)
    model.maximize(sum(sum(cols[k][i]*v[i] for i in range(m))*x[k] for k in range(len(cols))))
    res=mathopt.solve(model, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60),enable_output=False))
    dv=res.dual_values(); pi=[dv[con_item[i]] for i in range(m)]; mu=dv[con_sheet]
    return res.objective_value(), pi, mu
cols=[(0,)*m]+[tuple(1 if i==j else 0 for j in range(m)) for i in range(m)]
for it in range(100):
    lp,pi,mu=solve_rmp(cols)
    c=[v[i]-pi[i] for i in range(m)]
    pval,counts=best_pattern(c)
    rc=pval - mu
    print(f"it={it+1} lp={lp:.6f} mu={mu:.3f} pval={pval:.3f} rc={rc:.3e} ncols={len(cols)} counts={counts}")
    if rc<=1e-7:
        print("converged lp",lp); break
    if counts not in cols: cols.append(counts)
    else:
        print("dup"); break
