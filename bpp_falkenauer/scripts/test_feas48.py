import math, time
from ortools.sat.python import cp_model

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
K=48
model=cp_model.CpModel()
x=[[model.NewBoolVar(f"x_{i}_{b}") for b in range(K)] for i in range(n)]
for i in range(n):
    model.Add(sum(x[i][b] for b in range(K))==1)
for b in range(K):
    model.Add(sum(sizes[i]*x[i][b] for i in range(n))<=C)
# symmetry: bin b only items i>=b
for i in range(n):
    for b in range(i+1,K):
        model.Add(x[i][b]==0)
# force bins non-decreasing by first item? already implied by x[i][b]=0 for b>i.
solver=cp_model.CpSolver()
solver.parameters.max_time_in_seconds=120.0
solver.parameters.num_search_workers=8
solver.parameters.log_search_progress=False
t0=time.time()
status=solver.Solve(model)
print("status", solver.StatusName(status), "time", round(time.time()-t0,2), "wall", solver.WallTime())
if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    print("feasible 48 found!")
    asg=[0]*n
    for i in range(n):
        for b in range(K):
            if solver.Value(x[i][b])>0.5:
                asg[i]=b
    loads=[0]*K
    for i in range(n): loads[asg[i]]+=sizes[i]
    print("loads min/max", min(loads), max(loads), "all<=C", all(l<=C for l in loads))
