# -*- coding: utf-8 -*-
"""u120_00：CG 充分收敛（1500 轮）+ 划分池 MIP 找 48。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, time, datetime
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer/scripts")
import bpp_core as bc
from ortools.math_opt.python import mathopt

lp, patterns, sel, iters = bc.cg_min_rolls(tol=1e-7, max_iter=1500, verbose=False)
print(f"CG 收敛: {iters} 轮 | LP = {round(lp, 6)} | 模式 {len(sel)}")

# 划分（每件恰一次）池 MIP
mm = mathopt.Model()
x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=True, name=f"x{k}") for k in sel]
for i in range(bc.N):
    mm.add_linear_constraint(mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) == 1.0, name=f"c{i}")
mm.minimize(mathopt.fast_sum(x))
t0 = time.time()
res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                    params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=180), enable_output=False))
print(f"划分池 MIP: {res.termination.reason.name} | 箱数 {res.objective_value() if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE) else None} | {round(time.time()-t0,1)}s")
if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
    cnt = [0]*bc.N
    nz = 0
    for k in range(len(sel)):
        vv = res.variable_values()[x[k]]
        if vv > 0.5:
            nz += 1
            for i in range(bc.N):
                cnt[i] += patterns[sel[k]][i]
    print("每件恰一次:", all(c == 1 for c in cnt), "| 使用箱数:", nz, "| 总负载:", sum(bc.SIZES[i]*cnt[i] for i in range(bc.N)))
    print("最优性: 长度下界 48 =", res.objective_value(), "->", abs(res.objective_value() - 48) < 1e-6)
