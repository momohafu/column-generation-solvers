import time
from ortools.sat.python import cp_model
p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
total=sum(l*d for l,d in items); LB=(total+L-1)//L
print('LB',LB)
for K in [28,29]:
    model=cp_model.CpModel()
    y=[model.NewBoolVar(f'y_{k}') for k in range(K)]
    a=[[model.NewIntVar(0,d,f'a_{i}_{k}') for k in range(K)] for i,(l,d) in enumerate(items)]
    for k in range(K):
        model.Add(sum(items[i][0]*a[i][k] for i in range(m)) <= L*y[k])
    for i,(l,d) in enumerate(items):
        model.Add(sum(a[i][k] for k in range(K)) == d)
    for k in range(K-1):
        model.Add(y[k]>=y[k+1])
    model.Add(sum(y)==K)
    # decision strategy: branch on y to use rolls, then on largest item types
    model.AddDecisionStrategy(y, cp_model.CHOOSE_FIRST, cp_model.SELECT_MAX_VALUE)
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=60.0
    solver.parameters.num_search_workers=8
    solver.parameters.log_search_progress=False
    t0=time.time()
    status=solver.Solve(model)
    dt=time.time()-t0
    print('K',K,'status',solver.StatusName(status),'time',dt,'obj',solver.ObjectiveValue(),'bound',solver.BestObjectiveBound(),'branches',solver.NumBranches())
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        yv=[solver.Value(v) for v in y]
        print('y',yv)
        for k in range(K):
            if yv[k]:
                pat=[solver.Value(a[i][k]) for i in range(m)]
                ll=sum(items[i][0]*pat[i] for i in range(m))
                print(k,ll)
        break
