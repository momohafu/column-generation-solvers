from ortools.math_opt.python import mathopt
import inspect
m = mathopt.Model(name="t")
print(inspect.signature(m.add_variable))
print(m.add_variable.__doc__[:800])
# try float inf
try:
    v = m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name="x")
    print("float inf ok, ub=", v.upper_bound)
except Exception as e:
    print("float inf err", repr(e))
try:
    v2 = m.add_variable(lb=0.0, name="y")
    print("no ub ok, ub=", v2.upper_bound)
except Exception as e:
    print("no ub err", repr(e))
