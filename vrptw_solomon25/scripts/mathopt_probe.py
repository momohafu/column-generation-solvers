from ortools.math_opt.python import mathopt
import inspect
m=mathopt.Model(name="probe")
x=m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name="x")
y=m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name="y")
c1=m.add_linear_constraint(x+y>=1.0, name="c1")
m.minimize(x+2*y)
r=mathopt.solve(m, mathopt.SolverType.GLOP)
print("obj", r.objective_value())
print("varvals", r.variable_values())
print("dual", r.dual_values())
print("term", r.termination.reason)
print("has dual", hasattr(r, "dual_values"))
