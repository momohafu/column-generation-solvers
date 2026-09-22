# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer/scripts")
import bpp_core as bc
from ortools.math_opt.python import mathopt
import time

patterns = [bc.pattern_of([i]) for i in range(bc.N)]
sel = list(range(bc.N))
for it in range(300):
    try:
        mm = mathopt.Model()
        x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in sel]
        covers = []
        for i in range(bc.N):
            covers.append(mm.add_linear_constraint(
                mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= 1.0, name=f"c{i}"))
        mm.minimize(mathopt.fast_sum(x))
        res = mathopt.solve(mm, mathopt.SolverType.GLOP)
        dv = res.dual_values()
        pi = [max(0.0, dv[covers[i]]) for i in range(bc.N)]
        xvals = [res.variable_values()[x[k]] for k in range(len(sel))]
        added = 0
        pi_eff = list(pi)
        for g in range(8):
            v, isel = bc.knap_rebuild(pi_eff)
            rc = 1.0 - v
            if rc >= -1e-7:
                break
            pat = bc.pattern_of(isel)
            if pat not in patterns:
                patterns.append(pat)
                sel.append(len(patterns)-1)
                added += 1
            for i in isel:
                pi_eff[i] -= 0.05
        print(f"iter {it+1}: LP={round(res.objective_value(),5)} 加列 {added} 模式 {len(sel)}", flush=True)
        if added == 0:
            break
        if len(sel) > 400:
            keep = [k for k in range(len(xvals)) if xvals[k] > 1e-6] + list(range(len(xvals), len(sel)))
            if len(keep) < len(sel):
                sel = keep
                patterns = [patterns[k] for k in keep]
    except Exception as e:
        print(f"iter {it+1} CRASH: {type(e).__name__}: {str(e)[:150]}", flush=True)
        break
print("end", len(sel))
