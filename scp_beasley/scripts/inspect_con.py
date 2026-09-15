
from ortools.math_opt.python import mathopt
m=mathopt.Model(name='t')
x=m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='x')
c=m.add_linear_constraint(x >= 1.0, name='c')
print(type(c), [a for a in dir(c) if not a.startswith('_')])
