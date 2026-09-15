from ortools.math_opt.python import mathopt
import datetime
print(hasattr(mathopt,'SolveParameters'))
print(mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120.0)).time_limit)
m=mathopt.Model(name='t')
x=m.add_variable(lb=0.0,ub=1.0,is_integer=True,name='x')
m.maximize(x)
res=mathopt.solve(m, mathopt.SolverType.HIGHS, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=10.0), enable_output=False))
print(res.termination.reason, res.objective_value(), res.variable_values(x))
print(type(res.variable_values(x)))
print([a for a in dir(res) if 'bound' in a.lower()])
