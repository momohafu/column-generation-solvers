import json, datetime
from ortools.math_opt.python import mathopt
paths=[]
with open('/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts/paths.json') as f:
    for line in f:
        o=json.loads(line); paths.append(o)
N=50; K=27
m=mathopt.Model(name='spp_lp')
y=[m.add_variable(lb=0.0, ub=1.0, is_integer=False, name=f'y{p}') for p in range(len(paths))]
cover={i:[] for i in range(1,N+1)}
for p,o in enumerate(paths):
    for i in o['tasks']: cover[i].append(p)
for i in range(1,N+1):
    m.add_linear_constraint(sum(y[p] for p in cover[i]) == 1)
m.add_linear_constraint(sum(y) == K)
m.minimize(sum(paths[p]['cost']*y[p] for p in range(len(paths))))
res=mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=30), enable_output=False))
print('LP',res.termination.reason,'obj',res.objective_value())
