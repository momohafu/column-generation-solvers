from ortools.sat.python import cp_model
m=cp_model.CpModel()
print('circuit methods:', [x for x in dir(m) if 'Circuit' in x])
import ortools
print('ortools', ortools.__version__)
