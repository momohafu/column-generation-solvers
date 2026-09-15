import time
from ortools.sat.python import cp_model
p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
length=[l for l,d in items]; dem=[d for l,d in items]
exp=[]
for i,(l,d) in enumerate(items): exp += [l]*d
exp.sort(reverse=True)
bins=[]; assign=[]
for it in exp:
    best=None
    for i,b in enumerate(bins):
        if b+it<=L and (best is None or b>bins[best]): best=i
    if best is None: bins.append(it); assign.append([it])
    else: bins[best]+=it; assign[best].append(it)
print('FFD bins',len(bins))
K=29
model=cp_model.CpModel()
y=[model.NewBoolVar(f'y_{k}') for k in range(K)]
a=[[model.NewIntVar(0,d,f'a_{i}_{k}') for k in range(K)] for i,(l,d) in enumerate(items)]
for k in range(K):
    model.Add(sum(items[i][0]*a[i][k] for i in range(m)) <= L*y[k])
for i,(l,d) in enumerate(items):
    model.Add(sum(a[i][k] for k in range(K)) == d)
for k in range(K-1): model.Add(y[k]>=y[k+1])
model.Minimize(sum(y))
for k in range(K):
    if k < len(bins):
        model.AddHint(y[k], 1)
        counts=[0]*m
        for it in assign[k]:
            counts[length.index(it)] += 1
        for i in range(m): model.AddHint(a[i][k], counts[i])
    else:
        model.AddHint(y[k], 0)
        for i in range(m): model.AddHint(a[i][k], 0)
model.AddDecisionStrategy(y, cp_model.CHOOSE_FIRST, cp_model.SELECT_MAX_VALUE)
solver=cp_model.CpSolver()
solver.parameters.max_time_in_seconds=120.0
solver.parameters.num_search_workers=4
solver.parameters.log_search_progress=False
t0=time.time(); status=solver.Solve(model); dt=time.time()-t0
print('status',solver.StatusName(status),'time',dt,'obj',solver.ObjectiveValue(),'bound',solver.BestObjectiveBound(),'branches',solver.NumBranches())
if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    yv=[solver.Value(v) for v in y]
    print('used',sum(yv),'y',yv)
    for k in range(K):
        if yv[k]:
            pat=[solver.Value(a[i][k]) for i in range(m)]
            print(k,sum(items[i][0]*pat[i] for i in range(m)))
