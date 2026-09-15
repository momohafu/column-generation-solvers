import json, time, datetime
from ortools.math_opt.python import mathopt
paths=[]
with open('/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts/paths.json') as f:
    for line in f:
        o=json.loads(line); paths.append(o)
N=50; K=27
m=mathopt.Model(name='spp_opt')
y=[m.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f'y{p}') for p in range(len(paths))]
cover={i:[] for i in range(1,N+1)}
for p,o in enumerate(paths):
    for i in o['tasks']: cover[i].append(p)
for i in range(1,N+1):
    m.add_linear_constraint(sum(y[p] for p in cover[i]) == 1)
m.add_linear_constraint(sum(y) == K)
m.minimize(sum(paths[p]['cost']*y[p] for p in range(len(paths))))
t=time.time()
res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
wall=time.time()-t
print('term',res.termination.reason,'obj',res.objective_value() if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE) else None,'bound',res.best_objective_bound() if res.termination.reason!=mathopt.TerminationReason.INFEASIBLE else None,'wall',round(wall,3))
if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
    sol=[p for p in range(len(paths)) if res.variable_values()[y[p]]>0.5]
    print('routes',len(sol),'cost',sum(paths[p]['cost'] for p in sol))
    for p in sol: print(paths[p])
    with open('/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts/pool_opt_result.json','w') as f:
        json.dump({'K':K,'objective':sum(paths[p]['cost'] for p in sol),'status':str(res.termination.reason),'wall':wall,'routes':[paths[p] for p in sol]}, f, indent=2)
