import datetime, time
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
print('paths',P,'K',K)
# master constraints initially only sum y = K; coverage cuts added lazily
master_rows=[]  # list of (coeff dict, lb, ub) for linear rows over y
master_rows.append(({p:1.0 for p in range(P)}, K, K))
def solve_master(rows):
    m=mathopt.Model(name='master')
    y=[m.add_variable(lb=0.0,ub=1.0,is_integer=True,name=f'y{p}') for p in range(P)]
    for coeff,lb,ub in rows:
        terms=sum(coeff[p]*y[p] for p in range(P) if coeff.get(p,0)!=0)
        if lb==ub:
            m.add_linear_constraint(terms == lb)
        else:
            if lb is not None:
                m.add_linear_constraint(terms >= lb)
            if ub is not None:
                m.add_linear_constraint(terms <= ub)
    m.minimize(sum(paths[p]['cost']*y[p] for p in range(P)))
    res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=30), enable_output=False))
    if res.termination.reason not in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
        return None,res
    yv=[round(res.variable_values()[y[p]]) for p in range(P)]
    return yv,res
t0=time.time()
for it in range(1,200):
    yv,res=solve_master(master_rows)
    if yv is None:
        print('master fail',res.termination.reason); break
    obj=sum(paths[p]['cost']*yv[p] for p in range(P))
    # subproblem: check coverage of selected paths (LP feasibility)
    viol=[]
    for i in range(1,N+1):
        cnt=sum(1 for p in cover[i] if yv[p]==1)
        if cnt!=1:
            viol.append((i,cnt))
    if not viol:
        print('feasible at iter',it,'obj',obj)
        break
    # add violated coverage rows (Benders feasibility cuts)
    added=0
    for i,cnt in viol:
        coeff={p:1.0 for p in cover[i]}
        if cnt==0:
            row=(coeff,1,None)
        else:
            row=(coeff,None,1)
        # only add if not already present (simple duplicate check)
        key=(tuple(sorted(coeff)), row[1], row[2])
        existing=[ (tuple(sorted(rr[0].keys())), rr[1], rr[2]) for rr in master_rows ]
        if key not in existing:
            master_rows.append(row); added+=1
    print(f'iter {it}: obj {obj:.1f} selected {sum(yv)} violations {len(viol)} added {added} rows', flush=True)
    if added==0:
        # all violations already have rows but master still violates? shouldn't happen
        print('no new cuts, stop'); break
wall=time.time()-t0
print('Benders-lazy wall',round(wall,3),'iters',it,'master_rows',len(master_rows),'final_obj',obj if 'obj' in locals() else None)
