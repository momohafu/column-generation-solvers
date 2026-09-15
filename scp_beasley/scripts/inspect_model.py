
from ortools.math_opt.python import mathopt
m=mathopt.Model(name='t')
print([a for a in dir(m) if not a.startswith('_')])
