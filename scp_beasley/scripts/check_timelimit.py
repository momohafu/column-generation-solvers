import datetime
from ortools.math_opt.python import mathopt
import inspect
print('sig', inspect.signature(mathopt.SolveParameters))
try:
    p = mathopt.SolveParameters(time_limit=120.0, enable_output=False)
    print('float time_limit OK', p.time_limit)
except Exception as e:
    print('float time_limit ERROR:', type(e).__name__, e)
try:
    p = mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False)
    print('timedelta time_limit OK', p.time_limit)
except Exception as e:
    print('timedelta ERROR:', type(e).__name__, e)
