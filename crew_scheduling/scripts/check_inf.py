from ortools.math_opt.python import mathopt
print([x for x in dir(mathopt) if 'inf' in x.lower()])
print([x for x in dir(mathopt) if 'INF' in x])
import ortools.math_opt.python.mathopt as mo
print('inf?', getattr(mo,'inf',None), 'INF?', getattr(mo,'INF',None))
