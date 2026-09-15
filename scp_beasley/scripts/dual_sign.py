
from ortools.math_opt.python import mathopt
m=mathopt.Model(name='t')
x=m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='x')
c=m.add_linear_constraint(x <= 1.0, name='ub')
m.minimize_linear_objective(-2*x)
res=mathopt.solve(m, mathopt.SolverType.GLOP)
print('obj', res.objective_value(), 'x', res.variable_values([x])[0])
dv=res.dual_values([c])[0]
print('dual ub', dv)
m2=mathopt.Model(name='t2')
x2=m2.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='x')
c2=m2.add_linear_constraint(x2 >= 1.0, name='lb')
m2.minimize_linear_objective(2*x2)
res2=mathopt.solve(m2, mathopt.SolverType.GLOP)
print('obj2', res2.objective_value(), 'x2', res2.variable_values([x2])[0], 'dual lb', res2.dual_values([c2])[0])
