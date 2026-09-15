
from ortools.math_opt.python import mathopt
import inspect
print('solve result methods:')
print([m for m in dir(mathopt.SolveResult) if not m.startswith('_')])
print('dual_values sig', inspect.signature(mathopt.SolveResult.dual_values))
print('reduced_costs sig', inspect.signature(mathopt.SolveResult.reduced_costs))
print('variable_values sig', inspect.signature(mathopt.SolveResult.variable_values))
