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
K=27; P=len(paths)
cover=[[] for _ in range(N+1)]
for p,o in enumerate(paths):
    for i in o['tasks']: cover[i].append(p)
# initial active: singleton columns (path with one task)
singletons=[p for p,o in enumerate(paths) if len(o['tasks'])==1]
active=set(singletons)
bigM=1_000_000.0
t0=time.time(); lp_objs=[]; art_vals=[]
for it in range(1,501):
    m=mathopt.Model(name='rmp')
    y=[]; idx={}
    for p in sorted(active):
        idx[p]=len(y); y.append(m.add_variable(lb=0.0,ub=1.0,is_integer=False,name=f'y{p}'))
    sp=[m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name=f'sp{i}') for i in range(1,N+1)]
    sm=[m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name=f'sm{i}') for i in range(1,N+1)]
    tp=m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name='tp')
    tm=m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name='tm')
    cov_cons={}
    for i in range(1,N+1):
        cov_cons[i]=m.add_linear_constraint(sum(y[idx[p]] for p in cover[i] if p in idx) + sp[i-1] - sm[i-1] == 1)
    count_cons=m.add_linear_constraint(sum(y) + tp - tm == K)
    m.minimize(sum(paths[p]['cost']*y[idx[p]] for p in idx) + bigM*(sum(sp)+sum(sm)+tp+tm))
    res=mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=30), enable_output=False))
    if res.termination.reason != mathopt.TerminationReason.OPTIMAL:
        print('LP fail',res.termination.reason); break
    obj=res.objective_value(); dv=res.dual_values()
    pi={i:dv[cov_cons[i]] for i in range(1,N+1)}
    mu=dv[count_cons]
    art=sum(res.variable_values()[v] for v in sp+sm+[tp,tm])
    orig=sum(paths[p]['cost']*res.variable_values()[y[idx[p]]] for p in idx)
    lp_objs.append(orig); art_vals.append(art)
    # pricing: find most negative among not active
    best_rc=1e9; best_p=-1
    for p in range(P):
        if p in idx: continue
        rc=paths[p]['cost'] - sum(pi[i] for i in paths[p]['tasks']) - mu
        if rc < best_rc - 1e-12:
            best_rc=rc; best_p=p
    if it<=5 or it%20==0:
        print(f'iter {it}: active {len(idx)} orig_obj {orig:.6f} art {art:.2e} mu {mu:.2f} best_rc {best_rc:.6f}', flush=True)
    if best_p is None or best_rc >= -1e-7:
        print('terminate at iter',it,'best_rc',best_rc,'art',art)
        break
    active.add(best_p)
wall=time.time()-t0
print('CG done wall',round(wall,3),'iters',it,'active',len(active),'final orig obj',lp_objs[-1] if lp_objs else None,'final art',art_vals[-1] if art_vals else None)

# integer recovery on full pool
m=mathopt.Model(name='ip')
z=[m.add_variable(lb=0.0,ub=1.0,is_integer=True,name=f'z{p}') for p in range(P)]
for i in range(1,N+1):
    m.add_linear_constraint(sum(z[p] for p in cover[i]) == 1)
m.add_linear_constraint(sum(z) == K)
m.minimize(sum(paths[p]['cost']*z[p] for p in range(P)))
t1=time.time()
res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
print('IP',res.termination.reason,'obj',res.objective_value(),'wall',round(time.time()-t1,3))
if res.termination.reason==mathopt.TerminationReason.OPTIMAL:
    sol=[p for p in range(P) if res.variable_values()[z[p]]>0.5]
    print('IP routes',len(sol),'cost',sum(paths[p]['cost'] for p in sol))
