from ortools.math_opt.python import mathopt
print([a for a in dir(mathopt) if a.isupper() or 'inf' in a.lower()])
