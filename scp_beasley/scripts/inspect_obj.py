
from ortools.math_opt.python import mathopt
import inspect
m=mathopt.Model(name='t')
print('set_linear_objective', inspect.signature(m.set_linear_objective))
print('minimize_linear_objective', inspect.signature(m.minimize_linear_objective))
print('minimize', inspect.signature(m.minimize))
