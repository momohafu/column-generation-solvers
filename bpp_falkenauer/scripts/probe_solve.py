from ortools.math_opt.python import mathopt
import datetime
m = mathopt.Model(name="t")
x = m.add_variable(lb=0.0, is_integer=False, name="x")
m.minimize(x)
m.add_linear_constraint(x >= 3.0)
res = mathopt.solve(m, mathopt.SolverType.GLOP, params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
print("attrs:", [a for a in dir(res) if not a.startswith("_")])
print("termination reason", res.termination.reason)
print("objective", res.objective_value())
print("variable_values type", type(res.variable_values()))
print("has_feasible", getattr(res, "has_feasible_solution", "NA"))
try:
    print("best bound", res.best_objective_bound())
except Exception as e:
    print("best bound err", repr(e))
try:
    print("dual keys", list(res.dual_values().keys()))
except Exception as e:
    print("dual err", repr(e))
