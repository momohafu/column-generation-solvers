
from ortools.math_opt.python import mathopt
print([a for a in dir(mathopt) if 'inf' in a.lower()])
print(mathopt.__dict__.get('INFINITY'))
