import time, datetime, math
from ortools.math_opt.python import mathopt

DATA='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(DATA).read().split(); m,n=map(int,toks[:2]); costs=list(map(int,toks[2:2+n])); idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1; rows.append([int(t)-1 for t in toks[idx:idx+k]]); idx+=k
colrows=[[] for _ in range(n)]
for i,row in enumerate(rows):
    for j in row: colrows[j].append(i)

def ratio_greedy(selected):
    sel=set(selected); covered=[False]*m
    for j in sel:
        for i in colrows[j]: covered[i]=True
    while not all(covered):
        best=None
        for j in range(n):
            if j in sel: continue
            new=sum(1 for i in colrows[j] if not covered[i])
            if new==0: continue
            key=costs[j]/new
            if best is None or key<best[0]: best=(key,j,new)
        if best is None: return None
        j=best[1]; sel.add(j)
        for i in colrows[j]: covered[i]=True
    return list(sel)

def rc_greedy(pi, selected):
    sel=set(selected); covered=[False]*m
    for j in sel:
        for i in colrows[j]: covered[i]=True
    while not all(covered):
        best=None
        for j in range(n):
            if j in sel: continue
            new=0; s=0.0
            for i in colrows[j]:
                if not covered[i]: new+=1
                s+=pi[i]
            if new==0: continue
            rc=costs[j]-s
            key=(rc, costs[j], -new)
            if best is None or key<best[0]: best=(key,j)
        if best is None: return None
        j=best[1]; sel.add(j)
        for i in colrows[j]: covered[i]=True
    return list(sel)

UB=sum(costs[j] for j in ratio_greedy([]))
print('initial ratio greedy UB', UB)
pi=[0.0]*m; lam=2.0; best_L=float('-inf'); best_pi=None; best_x=None; no_impr=0; max_iter=600
t0=time.perf_counter()
for it in range(max_iter):
    rc=[0.0]*n
    for j in range(n):
        s=0.0
        for i in colrows[j]: s+=pi[i]
        rc[j]=costs[j]-s
    x=[1 if rc[j]<0 else 0 for j in range(n)]
    L=sum(pi)+sum(min(0.0, rc[j]) for j in range(n))
    if L>best_L+1e-9:
        best_L=L; best_pi=pi[:]; best_x=x[:]; no_impr=0
    else:
        no_impr+=1
    g=[1.0]*m
    for j in range(n):
        if x[j]:
            for i in colrows[j]: g[i]-=1.0
    norm2=sum(gi*gi for gi in g)
    if it % 50 == 0 or it == max_iter-1:
        sel=ratio_greedy([j for j in range(n) if x[j]])
        if sel is not None:
            c=sum(costs[j] for j in sel)
            if c<UB: UB=c
    step=0.0 if norm2<1e-12 else lam*(UB-L)/norm2
    for i in range(m): pi[i]=max(0.0, pi[i]+step*g[i])
    if no_impr>=60:
        lam=lam/2.0; no_impr=0
wall=time.perf_counter()-t0
sel_ratio=ratio_greedy([j for j in range(n) if best_x[j]])
cost_ratio=sum(costs[j] for j in sel_ratio) if sel_ratio else None
if cost_ratio is not None and cost_ratio<UB: UB=cost_ratio
sel_rc=rc_greedy(best_pi, [j for j in range(n) if best_x[j]])
cost_rc=sum(costs[j] for j in sel_rc) if sel_rc else None
print('Lagrangian best lower bound', best_L)
print('ratio greedy repair', cost_ratio, 'cols', len(sel_ratio) if sel_ratio else None)
print('reduced-cost greedy repair', cost_rc, 'cols', len(sel_rc) if sel_rc else None)
print('best feasible UB during subgradient', UB)
print('subgradient wall', round(wall,3), 'iterations', max_iter, 'final lambda', lam)

# MIP repair on full pool
t1=time.perf_counter()
mip=mathopt.Model(name='scp41_lag_repair')
x=[mip.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'x{j}') for j in range(n)]
mip.minimize_linear_objective(sum(costs[j]*x[j] for j in range(n)))
for i,row in enumerate(rows):
    mip.add_linear_constraint(sum(x[j] for j in row) >= 1.0, name=f'cov{i}')
mres=mathopt.solve(mip, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
mv=mres.variable_values(x)
selm=[j for j in range(n) if mv[j]>0.5]
print('MIP repair', mres.termination.reason, mres.objective_value(), 'best_bound', mres.best_objective_bound(), 'cols', len(selm), 'wall', round(time.perf_counter()-t1,3))
mobj=mres.objective_value()
print('dual gap vs MIP repair', round((mobj-best_L)/mobj, 6))
print('selected_mip', sorted(selm))
