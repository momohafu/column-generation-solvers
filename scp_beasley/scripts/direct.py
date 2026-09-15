import time, datetime
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

t0=time.perf_counter()
model=mathopt.Model(name='scp41_direct')
x=[model.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'x{j}') for j in range(n)]
model.minimize(sum(costs[j]*x[j] for j in range(n)))
for i,row in enumerate(rows):
    model.add_linear_constraint(sum(x[j] for j in row) >= 1.0, name=f'cov{i}')
params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False)
res=mathopt.solve(model, mathopt.SolverType.HIGHS, params=params)
dt=time.perf_counter()-t0
vals=res.variable_values(x)
sel=[j for j in range(n) if vals[j]>0.5]
covered=[False]*m
for j in sel:
    for i in colrows[j]: covered[i]=True
print('termination', res.termination.reason)
print('objective', res.objective_value(), 'best_bound', res.best_objective_bound())
print('solve_time', res.solve_time(), 'wall', round(dt,3))
print('num_selected', len(sel), 'obj_check', sum(costs[j] for j in sel), 'covered', sum(covered), '/', m)
print('selected_sorted', sorted(sel))
