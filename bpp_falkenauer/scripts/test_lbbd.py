import time
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
name,C,n,best,sizes=parse(); C=150; n=120

def ffd():
    order=sorted(range(n),key=lambda i:-sizes[i])
    bins=[];asg={}
    for i in order:
        s=sizes[i];pl=False
        for b in range(len(bins)):
            if bins[b]>=s: bins[b]-=s;asg[i]=b;pl=True;break
        if not pl: asg[i]=len(bins);bins.append(C-s)
    return len(bins),asg
K,_=ffd()
print("K (FFD)", K)

model=cp_model.CpModel()
x=[[model.NewBoolVar(f"x_{i}_{b}") for b in range(K)] for i in range(n)]
y=[model.NewBoolVar(f"y_{b}") for b in range(K)]
for i in range(n):
    model.Add(sum(x[i][b] for b in range(K))==1)
for b in range(K):
    for i in range(n):
        model.Add(x[i][b]<=y[b])
for b in range(K-1):
    model.Add(y[b]>=y[b+1])
for i in range(n):
    for b in range(i+1,K):
        model.Add(x[i][b]==0)
model.Minimize(sum(y))

solver=cp_model.CpSolver()
solver.parameters.max_time_in_seconds=30.0
solver.parameters.num_search_workers=8

t0=time.time()
for it in range(300):
    status=solver.Solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        print("master status", solver.StatusName(status)); break
    obj=round(solver.ObjectiveValue())
    # extract assignment
    asg=[-1]*n
    for i in range(n):
        for b in range(K):
            if solver.Value(x[i][b])>0.5:
                asg[i]=b
    loads=[0]*K
    for i in range(n):
        if asg[i]>=0: loads[asg[i]]+=sizes[i]
    viol=[b for b in range(K) if loads[b]>C]
    if not viol:
        print("LBBD feasible at iter", it, "obj", obj, "time", round(time.time()-t0,2))
        break
    # add no-good cuts for violated bins
    for b in viol:
        S=[i for i in range(n) if asg[i]==b]
        model.Add(sum(1-x[i][b] for i in S)>=1)
    if it%10==0:
        print("iter", it, "obj", obj, "viol", len(viol), "loads>", [ (b,loads[b]) for b in viol[:3] ], "time", round(time.time()-t0,2))
else:
    print("iteration cap reached")
print("final status", solver.StatusName(status), "total time", round(time.time()-t0,2))
