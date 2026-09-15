import datetime
from ortools.math_opt.python import mathopt
m=mathopt.Model(name='lp')
x=[m.add_variable(lb=0,ub=1,is_integer=False,name=f'x{i}') for i in range(3)]
m.add_linear_constraint(x[0]+x[1] == 1)
m.add_linear_constraint(x[0]+x[1]+x[2] == 2)
m.minimize(2*x[0]+3*x[1]+x[2])
res=mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=5), enable_output=False))
print('obj',res.objective_value())
print('has dual_values', hasattr(res,'dual_values'))
if hasattr(res,'dual_values'):
    dv=res.dual_values()
    print(type(dv), dv)
    try: print('dict', dv)
    except Exception as e: print('err',e)
    print('dir sample', [a for a in dir(dv) if not a.startswith('_')][:20])
