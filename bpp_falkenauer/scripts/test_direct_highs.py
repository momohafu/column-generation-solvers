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
name,C,n,best,sizes=parse()

def ffd():
    order=sorted(range(n),key=lambda i:-sizes[i])
    bins=[];asg={}
    for i in order:
        s=sizes[i];pl=False
        for b in range(len(bins)):
            if bins[b]>=s: bins[b]-=s;asg[i]=b;pl=True;break
        if not pl: asg[i]=len(bins);bins.append(C-s)
    return len(bins),asg
ub,_=ffd()
K=ub
model=mathopt.Model(name="direct_mip")
x=[[model.add_variable(lb=0.0,ub=1.0,is_integer=True) for b in range(K)] for i in range(n)]
y=[model.add_variable(lb=0.0,ub=1.0,is_integer=True) for b in range(K)]
for i in range(n):
    model.add_linear_constraint(mathopt.fast_sum(x[i][b] for b in range(K))==1.0)
for b in range(K):
    model.add_linear_constraint(mathopt.fast_sum(sizes[i]*x[i][b] for i in range(n))<=C*y[b])
for b in range(K-1):
    model.add_linear_constraint(y[b]>=y[b+1])
for i in range(n):
    for b in range(i+1,K):
        model.add_linear_constraint(x[i][b]==0.0)
model.add_linear_constraint(mathopt.fast_sum(y)>=48)
model.minimize(mathopt.fast_sum(y))
t0=time.time()
res=mathopt.solve(model,mathopt.SolverType.HIGHS,params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=90),enable_output=False))
print("status",res.termination.reason,"time",round(time.time()-t0,2))
if res.has_primal_feasible_solution:
    print("obj",res.objective_value(),"bound",res.best_objective_bound())
