import time, datetime
from ortools.math_opt.python import mathopt

p='/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt'
toks=open(p).read().split(); m,n=map(int,toks[:2]); costs=list(map(int,toks[2:2+n])); idx=2+n
rows=[]
for _ in range(m):
    k=int(toks[idx]); idx+=1; rows.append([int(t)-1 for t in toks[idx:idx+k]]); idx+=k
cut_rows=set()
t0=time.perf_counter()
max_iter=50
for it in range(max_iter):
    model=mathopt.Model(name='lbbd_master')
    y=[model.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'y{j}') for j in range(n)]
    model.minimize_linear_objective(sum(costs[j]*y[j] for j in range(n)))
    for r in cut_rows:
        model.add_linear_constraint(sum(y[j] for j in rows[r]) >= 1.0, name=f'cut{r}')
    res=mathopt.solve(model, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
    yv=res.variable_values(y)
    sel=[j for j in range(n) if yv[j]>0.5]
    uncovered=[]
    for r,row in enumerate(rows):
        if not any(yv[j]>0.5 for j in row):
            uncovered.append(r)
    print('it',it,'obj',res.objective_value(),'sel',len(sel),'uncovered',len(uncovered),'cuts',len(cut_rows))
    if not uncovered:
        print('covered all; break')
        break
    for r in uncovered:
        cut_rows.add(r)
    if time.perf_counter()-t0>110:
        print('time limit'); break
print('final obj', res.objective_value(), 'iterations', it+1, 'total cuts', len(cut_rows), 'wall', round(time.perf_counter()-t0,3), 'term', res.termination.reason)
print('selected', sorted(sel))
