from ortools.math_opt.python import mathopt
m=mathopt.Model(name='t')
try:
    x=m.add_variable(lb=0.0, ub=float('inf'), is_integer=False, name='x')
    print('ok', x)
except Exception as e:
    print('err', type(e), e)
