import time, datetime, math
from ortools.math_opt.python import mathopt
import sys
sys.path.insert(0,"/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg import column_generation, integer_recovery
info=column_generation(verbose=False)
print("pool", len(info["paths"]), "cg selected", len(info["selected"]))
t0=time.time()
ip=integer_recovery(info, time_limit=120, full_pool=True, verbose=True)
print("full IP result", ip)
print("wall", time.time()-t0)
