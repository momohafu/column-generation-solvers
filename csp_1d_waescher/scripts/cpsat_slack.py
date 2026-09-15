import time
from ortools.sat.python import cp_model
p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
length=[l for l,d in items]; dem=[d for l,d in items]
total=sum(l*d for l,d in items)
K=28
model=cp_model.CpModel()
a=[[model.NewIntVar(0,d,f'a_{i}_{k}') for k in range(K)] for i,(l,d) in enumerate(items)]
load=[model.NewIntVar(0,L,f'load_{k}') for k in range(K)]
slack=[model.NewIntVar(0,L,f'slack_{k}') for k in range(K)]
for k in range(K):
    model.Add(sum(items[i][0]*a[i][k] for i in range(m)) == load[k])
    model.Add(load[k]+slack[k] == L)
for i,(l,d) in enumerate(items):
    model.Add(sum(a[i][k] for k in range(K)) == d)
model.Add(sum(slack[k] for k in range(K)) == K*L-total)
# sort loads ascending (equivalently slack descending) to break symmetry
for k in range(K-1):
    model.Add(load[k] <= load[k+1])
# decision strategy on slack variables (small domain 0..65?) actually slack unbounded but sum 65, each<=65 by implication? add upper bound 65
for k in range(K):
    model.Add(slack[k] <= 65)
model.AddDecisionStrategy(slack, cp_model.CHOOSE_FIRST, cp_model.SELECT_MIN_VALUE)
solver=cp_model.CpSolver()
solver.parameters.max_time_in_seconds=120.0
solver.parameters.num_search_workers=4
solver.parameters.log_search_progress=False
t0=time.time(); status=solver.Solve(model); dt=time.time()-t0
print('K',K,'status',solver.StatusName(status),'time',dt,'branches',solver.NumBranches(),'conflicts',solver.NumConflicts())
if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    for k in range(K):
        pat=[solver.Value(a[i][k]) for i in range(m)]
        print(k,'load',solver.Value(load[k]),'slack',solver.Value(slack[k]))
