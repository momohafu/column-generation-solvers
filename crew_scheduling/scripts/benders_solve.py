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
print('paths',P,'K',K)

def solve_master(cuts):
    m=mathopt.Model(name='bm')
    y=[m.add_variable(lb=0.0,ub=1.0,is_integer=True,name=f'y{p}') for p in range(P)]
    theta=m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name='theta')
    m.add_linear_constraint(sum(y) == K)
    for const,coef in cuts:
        m.add_linear_constraint(theta - const + sum(coef[p]*y[p] for p in range(P)) >= 0)
    m.minimize(sum(paths[p]['cost']*y[p] for p in range(P)) + theta)
    res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=30), enable_output=False))
    if res.termination.reason not in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
        return None,None,res
    yv=[round(res.variable_values()[y[p]]) for p in range(P)]
    return yv,res.variable_values()[theta],res

def solve_subproblem(yhat):
    m=mathopt.Model(name='bs')
    x=[]; xcons=[]
    for p in range(P):
        xv=m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name=f'x{p}')
        x.append(xv); xcons.append(m.add_linear_constraint(xv <= yhat[p]))
    sp=[m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name=f'sp{i}') for i in range(1,N+1)]
    sm=[m.add_variable(lb=0.0,ub=float('inf'),is_integer=False,name=f'sm{i}') for i in range(1,N+1)]
    eq=[]
    for i in range(1,N+1):
        eq.append(m.add_linear_constraint(sum(x[p] for p in cover[i]) + sp[i-1] - sm[i-1] == 1))
    m.minimize(sum(sp)+sum(sm))
    res=mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=30), enable_output=False))
    if res.termination.reason != mathopt.TerminationReason.OPTIMAL:
        return None,None,None,res
    dv=res.dual_values()
    pi={i:dv[eq[i-1]] for i in range(1,N+1)}
    alpha={p:-dv[xcons[p]] for p in range(P)}
    return res.objective_value(),pi,alpha,res

cuts=[]; best_ub=None; best_routes=None
t0=time.time()
for it in range(1,81):
    yv,th,mres=solve_master(cuts)
    if yv is None:
        print('master fail',mres.termination.reason); break
    V,pi,alpha,sres=solve_subproblem(yv)
    selected=[p for p in range(P) if yv[p]==1]
    obj=sum(paths[p]['cost'] for p in selected)
    print(f'iter {it}: master_obj {obj+th:.3f} theta {th:.3f} selected {len(selected)} sub_penalty {V:.6f}', flush=True)
    # repair selected set to integer partition if every task covered by at least one selected path
    covered_all=all(any(p in selected for p in cover[i]) for i in range(1,N+1))
    if covered_all:
        m=mathopt.Model(name='repair')
        z=[m.add_variable(lb=0.0,ub=1.0,is_integer=True,name=f'z{p}') for p in selected]
        smap={p:k for k,p in enumerate(selected)}
        for i in range(1,N+1):
            m.add_linear_constraint(sum(z[smap[p]] for p in cover[i] if p in smap) == 1)
        m.add_linear_constraint(sum(z) == K)
        m.minimize(sum(paths[p]['cost']*z[smap[p]] for p in selected))
        rr=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=20), enable_output=False))
        if rr.termination.reason==mathopt.TerminationReason.OPTIMAL:
            ub=rr.objective_value(); routes=[p for p in selected if rr.variable_values()[z[smap[p]]]>0.5]
            if best_ub is None or ub<best_ub:
                best_ub=ub; best_routes=routes
                print('  repair ub',ub,'routes',len(routes))
    if V is not None and V < 1e-6:
        print('  subproblem LP feasible (theta convergence)')
        if best_ub is not None:
            break
    const=sum(pi.values()); coef={p:alpha[p] for p in range(P)}
    cuts.append((const,coef))
    if it<=5 or it%10==0:
        print('  cut',len(cuts),'const',round(const,3),'sum_alpha',round(sum(alpha.values()),3))
wall=time.time()-t0
print('Benders wall',round(wall,3),'cuts',len(cuts),'ub',best_ub)
if best_ub is not None:
    print('best cost',sum(paths[p]['cost'] for p in best_routes),'routes',len(best_routes))
    for p in best_routes: print(paths[p])
