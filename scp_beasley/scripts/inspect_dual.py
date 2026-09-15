
from ortools.math_opt.python import mathopt
m=mathopt.Model(name='t')
x=m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='x')
c=m.add_linear_constraint(x >= 1.0, name='c')
m.minimize(2*x)
res=mathopt.solve(m, mathopt.SolverType.GLOP)
print('obj', res.objective_value())
dv=res.dual_values()
print('dual type', type(dv), dv)
print('dual for c', dv.get(c))
print('vars', res.variable_values())
print('reduced', res.reduced_costs())
