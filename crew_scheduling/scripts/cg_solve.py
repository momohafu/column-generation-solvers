import datetime, time, json
from ortools.math_opt.python import mathopt

DATA='/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt'
lines=open(DATA).read().splitlines()
N,T=map(int,lines[0].split())
tasks=[None]
for line in lines[1:1+N]:
    s,f=map(int,line.split()); tasks.append((s,f))
arcs={}
for line in lines[1+N:]:
    i,j,c=map(int,line.split()); arcs[(i,j)]=c
adj={i:[] for i in range(1,N+1)}
for (i,j) in arcs: adj[i].append(j)
paths=[]
for start in range(1,N+1):
    s0=tasks[start][0]
    stack=[(start,[start],0)]
    while stack:
        u,seq,c=stack.pop()
        paths.append({'tasks':list(seq),'cost':c,'span':tasks[u][1]-s0})
        for v in adj[u]:
            c2=c+arcs[(u,v)]
            if tasks[v][1]-s0<=T:
                stack.append((v,seq+[v],c2))
K=27
P=len(paths)
print('paths',P,'K',K)
bigM=1_000_000.0
cover=[[] for _ in range(N+1)]
for p,o in enumerate(paths):
    for i in o['tasks']: cover[i].append(p)
active=[]
t0=time.time()
for it in range(1,201):
    # build RMP with active cols + artificials
    m=mathopt.Model(name='rmp')
    y=[]
    idx={}
    for a,p in enumerate(active):
        y.append(m.add_variable(lb=0.0, ub=1.0, is_integer=False, name=f'y{p}'))
        idx[p]=a
    sp=[m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name=f'sp{i}') for i in range(1,N+1)]
    sm=[m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name=f'sm{i}') for i in range(1,N+1)]
    tp=m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='tp')
    tm=m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='tm')
    cov_cons={}
    for i in range(1,N+1):
        terms=sum(y[idx[p]] for p in cover[i] if p in idx)
        cov_cons[i]=m.add_linear_constraint(terms + sp[i-1] - sm[i-1] == 1)
    count_cons=m.add_linear_constraint(sum(y) + tp - tm == K)
    m.minimize(sum(paths[p]['cost']*y[idx[p]] for p in active) + bigM*(sum(sp)+sum(sm)+tp+tm))
    res=mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=30), enable_output=False))
    if res.termination.reason != mathopt.TerminationReason.OPTIMAL:
        print('LP fail',res.termination.reason); break
    obj=res.objective_value(); dv=res.dual_values()
    pi={i:dv[cov_cons[i]] for i in range(1,N+1)}
    mu=dv[count_cons]
    art=sum(res.variable_values()[v] for v in sp+sm+[tp,tm])
    if it<=3 or it%10==0:
        print(f'iter {it}: active {len(active)} lpobj {obj:.6f} art {art:.2e} mu {mu:.4f}', flush=True)
    # pricing over all pool columns
    best_rc=0.0; best_p=-1; neg=[]
    for p in range(P):
        rc=paths[p]['cost'] - sum(pi[i] for i in paths[p]['tasks']) - mu
        if rc < -1e-7:
            neg.append((rc,p))
        if rc < best_rc-1e-12:
            best_rc=rc; best_p=p
    if best_p<0:
        print('no negative reduced cost columns at iter',it,'best_rc',best_rc)
        print('LP objective original',sum(paths[p]['cost']*res.variable_values()[y[idx[p]]] for p in active),'art',art)
        break
    # add all negative columns
    added=0
    for rc,p in neg:
        if p not in idx:
            active.append(p); added+=1
    print(f'iter {it}: added {added} columns (min rc {best_rc:.6f})')
else:
    print('reached max iter')
wall=time.time()-t0
print('CG wall',wall,'iterations',it,'active',len(active))

# integer recovery on full pool (complete pool => exact)
m=mathopt.Model(name='ip')
y=[m.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'z{p}') for p in range(P)]
for i in range(1,N+1):
    m.add_linear_constraint(sum(y[p] for p in cover[i]) == 1)
m.add_linear_constraint(sum(y) == K)
m.minimize(sum(paths[p]['cost']*y[p] for p in range(P)))
t1=time.time()
res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
print('IP',res.termination.reason,'obj',res.objective_value(),'wall',time.time()-t1)
if res.termination.reason==mathopt.TerminationReason.OPTIMAL:
    sol=[p for p in range(P) if res.variable_values()[y[p]]>0.5]
    print('IP routes',len(sol),'cost',sum(paths[p]['cost'] for p in sol))
    with open('/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts/cg_result.json','w') as f:
        json.dump({'iterations':it,'lp_obj':None,'ip_obj':res.objective_value(),'active':active,'routes':[paths[p] for p in sol]}, f, indent=2)
