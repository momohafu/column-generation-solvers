
from ortools.math_opt.python import mathopt
m=mathopt.Model(name='t')
x=m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='x')
c1=m.add_linear_constraint(x >= 1.0, name='c1')
c2=m.add_linear_constraint(2*x >= 1.0, name='c2')
m.minimize(3*x)
res=mathopt.solve(m, mathopt.SolverType.GLOP)
print('dual dict', res.dual_values())
print('dual list', res.dual_values([c1,c2]))
print('vals list', res.variable_values([x]))
