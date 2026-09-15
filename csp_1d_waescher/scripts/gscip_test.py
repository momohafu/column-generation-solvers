import time, datetime
from ortools.math_opt.python import mathopt
p='/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt'
lines=open(p,encoding='utf-8').read().splitlines()
m=int(lines[0].strip()); L=int(lines[1].strip())
items=[]
for line in lines[2:2+m]:
    a=line.split(); items.append((int(a[0]),int(a[1])))
K=28
model=mathopt.Model(name='csp_direct')
y=[model.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'y_{k}') for k in range(K)]
x=[[model.add_variable(lb=0.0, ub=float(d), is_integer=True, name=f'x_{i}_{k}') for k in range(K)] for i,(l,d) in enumerate(items)]
for k in range(K):
    model.add_linear_constraint(sum(items[i][0]*x[i][k] for i in range(m)) <= L*y[k])
for i,(l,d) in enumerate(items):
    model.add_linear_constraint(sum(x[i][k] for k in range(K)) == d)
for k in range(K-1):
    model.add_linear_constraint(y[k] >= y[k+1])
model.minimize(sum(y[k] for k in range(K)))
t0=time.time()
res=mathopt.solve(model, mathopt.SolverType.GSCIP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120.0), enable_output=False))
dt=time.time()-t0
print('term',res.termination.reason,'time',dt,'primal',res.primal_bound,'dual',res.dual_bound,'best_bound',res.best_objective_bound)
print('solutions',len(res.solutions))
