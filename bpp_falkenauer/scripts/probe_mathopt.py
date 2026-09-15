from ortools.math_opt.python import mathopt
print([a for a in dir(mathopt) if not a.startswith("_")])
print("has inf:", hasattr(mathopt, "inf"))
try:
    print("inf:", mathopt.inf)
except Exception as e:
    print("inf err", e)
