import math, time, datetime
from ortools.math_opt.python import mathopt

def parse():
    with open("/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt") as f:
        lines=[l.rstrip("\n") for l in f]
    P=int(lines[0].strip()); idx=1
    name=lines[idx].strip(); idx+=1
    parts=lines[idx].split(); idx+=1
    C=int(parts[0]); n=int(parts[1]); best=int(parts[2])
    sizes=[int(lines[idx+i]) for i in range(n)]
    return name,C,n,best,sizes

name,C,n,best,sizes = parse()
columns=[tuple([i]) for i in range(n)]
colset=set(columns)
def build_rmp(columns):
    model=mathopt.Model(name="rmp")
    lam=[model.add_variable(lb=0.0,is_integer=False) for _ in columns]
    cons=[]
    for i in range(n):
        cons.append(model.add_linear_constraint(mathopt.fast_sum(lam[j] for j,col in enumerate(columns) if i in col)==1.0))
    model.minimize(mathopt.fast_sum(lam))
    return model,lam,cons
def knapsack_exact(profits,weights,capacity):
    n=len(weights); dp=[0.0]*(capacity+1); keep=[[False]*(capacity+1) for _ in range(n)]
    for i in range(n):
        w=weights[i]; p=profits[i]
        for c in range(capacity,w-1,-1):
            cand=dp[c-w]+p
            if cand>dp[c]+1e-12:
                dp[c]=cand; keep[i][c]=True
    bestc=max(range(capacity+1),key=lambda c:dp[c]); val=dp[bestc]
    items=[]; c=bestc
    for i in range(n-1,-1,-1):
        if keep[i][c]:
            items.append(i); c-=weights[i]
    return val,items
TOL=1e-7
for it in range(2000):
    model,lam,cons=build_rmp(columns)
    res=mathopt.solve(model,mathopt.SolverType.GLOP,params=mathopt.SolveParameters(enable_output=False))
    obj=res.objective_value(); pi=[res.dual_values()[cons[i]] for i in range(n)]
    profit,items=knapsack_exact(pi,sizes,C)
    rc=1.0-profit
    if rc<-TOL:
        newcol=tuple(sorted(items))
        if newcol in colset: break
        columns.append(newcol); colset.add(newcol)
    else: break
print("ncol",len(columns))
# IP
model=mathopt.Model(name="ip")
lam=[model.add_variable(lb=0.0,ub=float(n),is_integer=True) for _ in columns]
for i in range(n):
    model.add_linear_constraint(mathopt.fast_sum(lam[j] for j,col in enumerate(columns) if i in col)==1.0)
model.minimize(mathopt.fast_sum(lam))
res=mathopt.solve(model,mathopt.SolverType.HIGHS,params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120),enable_output=False))
vals=res.variable_values()
used=[j for j in range(len(columns)) if vals[lam[j]]>0.5]
print("used patterns",len(used))
# build bins: replicate each pattern by its integer usage
bins=[]
for j in used:
    cnt=round(vals[lam[j]])
    col=columns[j]
    for _ in range(cnt):
        bins.append(list(col))
print("total bins",len(bins))
# verify each item exactly once
from collections import Counter
c=Counter()
for b in bins:
    for i in b: c[i]+=1
print("items not exactly once:", [i for i in range(n) if c[i]!=1])
print("loads min/max", min(sum(sizes[i] for i in b) for b in bins), max(sum(sizes[i] for i in b) for b in bins))
# assignment
asg={}
for b,items in enumerate(bins):
    for i in items: asg[i]=b
with open("/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer/scripts/packing48.txt","w") as f:
    f.write("\n".join(str(asg[i]) for i in range(n)))
print("saved packing48.txt")
