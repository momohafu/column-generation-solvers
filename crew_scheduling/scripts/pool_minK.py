import json, time, datetime
from ortools.math_opt.python import mathopt
paths=[]
with open('/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts/paths.json') as f:
    for line in f:
        o=json.loads(line); paths.append(o)
N=50
m=mathopt.Model(name='minK')
y=[m.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'y{p}') for p in range(len(paths))]
cover={i:[] for i in range(1,N+1)}
for p,o in enumerate(paths):
    for i in o['tasks']: cover[i].append(p)
for i in range(1,N+1):
    m.add_linear_constraint(sum(y[p] for p in cover[i]) == 1)
m.minimize(sum(y))
res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=30), enable_output=False))
print('minK term',res.termination.reason,'obj',res.objective_value() if res.termination.reason==mathopt.TerminationReason.OPTIMAL else None)
if res.termination.reason==mathopt.TerminationReason.OPTIMAL:
    sol=[p for p in range(len(paths)) if res.variable_values()[y[p]]>0.5]
    print('routes',len(sol),'tasks',sum(len(paths[p]['tasks']) for p in sol))
    # print path lengths
    from collections import Counter
    print(Counter(len(paths[p]['tasks']) for p in sol))
