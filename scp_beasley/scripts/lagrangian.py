import time, datetime, math
from ortools.math_opt.python import mathopt

p='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(p).read().split()
m,n=map(int,toks[:2])
costs=list(map(int,toks[2:2+n]))
idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1
    rows.append([int(t)-1 for t in toks[idx:idx+k]]); idx+=k
colrows=[[] for _ in range(n)]
for i,row in enumerate(rows):
    for j in row: colrows[j].append(i)

def greedy_repair(selected):
    sel=set(selected)
    covered=[False]*m
    for j in sel:
        for i in colrows[j]: covered[i]=True
    # add cheapest per newly covered
    while not all(covered):
        best=None
        for j in range(n):
            if j in sel: continue
            new=[i for i in colrows[j] if not covered[i]]
            if not new: continue
            ratio=costs[j]/len(new)
            if best is None or ratio<best[0]:
                best=(ratio,j,len(new))
        if best is None:
            return None
        j=best[1]; sel.add(j)
        for i in colrows[j]: covered[i]=True
    return list(sel)

# initial greedy UB
sel0=greedy_repair([])
UB=sum(costs[j] for j in sel0)
print('greedy initial UB', UB, 'cols', len(sel0))

pi=[0.0]*m
lam=2.0
best_L=float('-inf')
best_pi=None
best_x=None
no_impr=0
max_iter=600
t0=time.perf_counter()
hist=[]
for it in range(max_iter):
    rc=[0.0]*n
    # compute rc_j = c_j - sum_i pi_i a_ij
    for j in range(n):
        s=0.0
        for i in colrows[j]:
            s += pi[i]
        rc[j]=costs[j]-s
    x=[1 if rc[j]<0 else 0 for j in range(n)]
    L=sum(pi)+sum(min(0.0, rc[j]) for j in range(n))
    if L>best_L+1e-9:
        best_L=L; best_pi=pi[:]; best_x=x[:]; no_impr=0
    else:
        no_impr+=1
    # subgradient
    g=[1.0]*m
    for j in range(n):
        if x[j]:
            for i in colrows[j]: g[i]-=1.0
    norm2=sum(gi*gi for gi in g)
    # periodically repair current x
    if it % 50 == 0 or it == max_iter-1:
        sel=greedy_repair([j for j in range(n) if x[j]])
        if sel is not None:
            cost=sum(costs[j] for j in sel)
            if cost<UB: UB=cost
    if norm2<1e-12:
        step=0.0
    else:
        step=lam*(UB-L)/norm2
    for i in range(m):
        pi[i]=max(0.0, pi[i]+step*g[i])
    if no_impr>=60:
        lam=lam/2.0; no_impr=0
    hist.append((it,L,UB,best_L,lam,norm2))
    if time.perf_counter()-t0>110:
        print('time limit at', it)
        break
wall=time.perf_counter()-t0
# final greedy from best_x
sel_best=greedy_repair([j for j in range(n) if best_x[j]])
greedy_best=sum(costs[j] for j in sel_best)
if greedy_best<UB: UB=greedy_best
print('Lagrangian best_L', best_L, 'UB', UB, 'gap', (UB-best_L)/UB if UB else None)
print('greedy_best cost', greedy_best, 'cols', len(sel_best))
print('wall', round(wall,3), 'iters', len(hist))
print('first hist', hist[:3])
print('last hist', hist[-3:])
# final MIP repair full pool
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
print('selected_mip', sorted(selm))
