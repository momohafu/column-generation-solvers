from ortools.math_opt.python import mathopt
m=mathopt.Model(name='t')
x=m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='x')
m.add_linear_constraint(x>=1.0)
m.minimize_linear_objective(x)
try:
    res=mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=120.0, enable_output=False))
    print('float solve OK', res.objective_value())
except Exception as e:
    print('float solve ERROR:', type(e).__name__, e)
