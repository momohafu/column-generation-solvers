from ortools.math_opt.python import mathopt
names=[a for a in dir(mathopt) if not a.startswith("_")]
print(len(names))
print(names[:120])
