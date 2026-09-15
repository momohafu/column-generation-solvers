from ortools.sat.python import cp_model
import inspect
m=cp_model.CpModel()
print(inspect.getdoc(m.AddCircuit))
print('---MULTI---')
print(inspect.getdoc(m.AddMultipleCircuit))
