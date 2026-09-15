import time, json
from ortools.sat.python import cp_model

DATA='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(DATA).read().splitlines()
N,T=map(int,lines[0].split())
tasks=[None]
for line in lines[1:1+N]:
    s,f=map(int,line.split()); tasks.append((s,f))
arcs={}
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); arcs[(i,j)]=c
arc_list=list(arcs.keys())
K=27
m=cp_model.CpModel()
x={}
for c in range(K):
    carcs=[]
    for i in range(1,N+1):
        x[(c,'0',i)]=m.NewBoolVar(f'c{c}_0_{i}')
        carcs.append((0,i,x[(c,'0',i)]))
        x[(c,'e',i)]=m.NewBoolVar(f'c{c}_{i}_e')
        carcs.append((i,N+1,x[(c,'e',i)]))
        x[(c,'s',i)]=m.NewBoolVar(f'c{c}_s_{i}')
        carcs.append((i,i,x[(c,'s',i)]))
    for (i,j) in arc_list:
        x[(c,'a',i,j)]=m.NewBoolVar(f'c{c}_{i}_{j}')
        carcs.append((i,j,x[(c,'a',i,j)]))
    x[(c,'empty')]=m.NewBoolVar(f'c{c}_empty')
    carcs.append((0,N+1,x[(c,'empty')]))
    x[(c,'close')]=m.NewBoolVar(f'c{c}_close')
    carcs.append((N+1,0,x[(c,'close')]))
    m.AddCircuit(carcs)
    m.Add(sum(tasks[i][1]*x[(c,'e',i)] for i in range(1,N+1)) - sum(tasks[i][0]*x[(c,'0',i)] for i in range(1,N+1)) <= T)
for i in range(1,N+1):
    m.Add(sum(x[(c,'s',i)] for c in range(K)) == K-1)
m.Minimize(sum(arcs[(i,j)]*x[(c,'a',i,j)] for c in range(K) for (i,j) in arc_list))
s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=120; s.parameters.num_workers=8; s.parameters.log_search_progress=True
t=time.time(); st=s.Solve(m); wall=time.time()-t
print('status',s.StatusName(st),'wall',wall,'obj',s.ObjectiveValue() if st in (cp_model.FEASIBLE,cp_model.OPTIMAL) else None,'bound',s.BestObjectiveBound())
if st in (cp_model.FEASIBLE,cp_model.OPTIMAL):
    routes=[]
    for c in range(K):
        first=None
        for i in range(1,N+1):
            if s.Value(x[(c,'0',i)])==1: first=i; break
        if first is None: continue
        seq=[]; cur=first
        while cur is not None:
            seq.append(cur)
            nxt=None
            if s.Value(x[(c,'e',cur)])==1: nxt=None
            else:
                for (i,j) in arc_list:
                    if i==cur and s.Value(x[(c,'a',i,j)])==1: nxt=j; break
            cur=nxt
        routes.append({'crew':c,'tasks':seq,'cost':sum(arcs[(seq[k],seq[k+1])] for k in range(len(seq)-1)),'span':tasks[seq[-1]][1]-tasks[seq[0]][0]})
    print('routes',len(routes),'total cost',sum(r['cost'] for r in routes),'span max',max(r['span'] for r in routes))
    for r in routes: print(r)
