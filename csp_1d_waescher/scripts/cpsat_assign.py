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
for i,(l,d) in enumerate(items):
    exp += [(l,i)]*d
n=len(exp)
exp.sort(reverse=True, key=lambda t:t[0])
print('n items',n)
K=28
model=cp_model.CpModel()
x=[[model.NewBoolVar(f'x_{r}_{k}') for k in range(K)] for r in range(n)]
load=[model.NewIntVar(0,L,f'load_{k}') for k in range(K)]
for r,(w,i) in enumerate(exp):
    model.Add(sum(x[r][k] for k in range(K)) == 1)
for k in range(K):
    model.Add(sum(exp[r][0]*x[r][k] for r in range(n)) == load[k])
    model.Add(load[k] <= L)
# symmetry: load non-increasing by sorting? not strict. add bin used order by first item
# force bins used in order: bin k used iff contains item 0? no.
# stronger: lexicographic by the largest item bin assignment? skip for now.
# decision strategy: assign largest items first, first bin
vars_strat=[]
for r in range(n):
    vars_strat += x[r]
model.AddDecisionStrategy(vars_strat, cp_model.CHOOSE_FIRST, cp_model.SELECT_MIN_VALUE)
solver=cp_model.CpSolver()
solver.parameters.max_time_in_seconds=120.0
solver.parameters.num_search_workers=8
solver.parameters.log_search_progress=False
t0=time.time(); status=solver.Solve(model); dt=time.time()-t0
print('K',K,'status',solver.StatusName(status),'time',dt,'branches',solver.NumBranches(),'conflicts',solver.NumConflicts())
if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    for k in range(K):
        pat=[0]*m
        for r in range(n):
            if solver.Value(x[r][k]): pat[exp[r][1]]+=1
        ll=sum(items[i][0]*pat[i] for i in range(m))
        print('bin',k,'load',ll)
