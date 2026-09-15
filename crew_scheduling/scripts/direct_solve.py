import json, time, sys
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
print('N',N,'T',T,'arcs',len(arc_list), flush=True)

# lower bound K
total_dur=sum(f-s for s,f in tasks[1:])
LB=max(1, (total_dur + T -1)//T)
print('total_dur',total_dur,'T',T,'LB_K_duration',(total_dur + T -1)//T, flush=True)
# overlap LB
ev=[]
for i,(s,f) in enumerate(tasks[1:],1):
    ev.append((s,1)); ev.append((f,-1))
ev.sort(key=lambda x:(x[0],-x[1]))
cur=0; overlap=0
for t,d in ev:
    cur+=d; overlap=max(overlap,cur)
print('overlap_LB',overlap, flush=True)
K=LB

def build(K, objective=True):
    m=cp_model.CpModel()
    # per crew nodes: 0 start, 1..N tasks, N+1 end
    x={}
    def nv(c,name):
        v=m.NewBoolVar(f'c{c}_{name}')
        x[(c,name)]=v
        return v
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
        # span <= T (depot travel times zero)
        start_expr=sum(tasks[i][0]*x[(c,'0',i)] for i in range(1,N+1))
        end_expr=sum(tasks[i][1]*x[(c,'e',i)] for i in range(1,N+1))
        m.Add(end_expr - start_expr <= T)
    # each task exactly once globally
    for i in range(1,N+1):
        m.Add(sum(x[(c,'s',i)] for c in range(K)) == K-1)
    if objective:
        m.Minimize(sum(arcs[(i,j)]*x[(c,'a',i,j)] for c in range(K) for (i,j) in arc_list))
    return m, x

def extract(solver, K, x):
    routes=[]
    for c in range(K):
        # find first task
        first=None; last=None; seq=[]
        for i in range(1,N+1):
            if solver.Value(x[(c,'0',i)])==1:
                first=i; break
        if first is None:
            continue
        cur=first
        while cur is not None:
            seq.append(cur)
            # next
            nxt=None
            if solver.Value(x[(c,'e',cur)])==1:
                nxt=None
            else:
                for (i,j) in arc_list:
                    if i==cur and solver.Value(x[(c,'a',i,j)])==1:
                        nxt=j; break
            cur=nxt
        span=tasks[seq[-1]][1]-tasks[seq[0]][0] if seq else 0
        cost=sum(arcs[(seq[k],seq[k+1])] for k in range(len(seq)-1))
        routes.append({'crew':c,'tasks':seq,'span':span,'cost':cost,'start':tasks[seq[0]][0] if seq else 0,'finish':tasks[seq[-1]][1] if seq else 0})
    return routes

# feasibility search for min K (no objective) quick
feasible_K=None
for k in range(LB, N+1):
    m,x=build(k, objective=False)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=30; s.parameters.num_workers=8; s.parameters.log_search_progress=False
    st=s.Solve(m)
    print('feas K',k,'->',s.StatusName(st),'wall',s.WallTime(), flush=True)
    if st in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        feasible_K=k; break
    if st==cp_model.UNKNOWN:
        print('UNKNOWN at K',k,'continuing search upward?', flush=True)
        break
print('feasible_K',feasible_K, flush=True)

# optimize for feasible_K
if feasible_K is not None:
    m,x=build(feasible_K, objective=True)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=120; s.parameters.num_workers=8; s.parameters.log_search_progress=True
    t0=time.time()
    st=s.Solve(m)
    wall=time.time()-t0
    print('OPT status',s.StatusName(st),'wall',wall,'obj',s.ObjectiveValue() if st in (cp_model.FEASIBLE,cp_model.OPTIMAL) else None,'best_bound',s.BestObjectiveBound() if hasattr(s,'BestObjectiveBound') else None, flush=True)
    if st in (cp_model.FEASIBLE,cp_model.OPTIMAL):
        routes=extract(s, feasible_K, x)
        print(json.dumps({'K':feasible_K,'objective':s.ObjectiveValue(),'status':s.StatusName(st),'routes':routes}, indent=2))
        with open('/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts/direct_result.json','w') as f:
            json.dump({'K':feasible_K,'objective':s.ObjectiveValue(),'status':s.StatusName(st),'wall':wall,'routes':routes}, f, indent=2)
