# -*- coding: utf-8 -*-
"""GSCIP K=28：item->bin MIP + 每箱负载 ∈ [9935,10000] + min 空隙。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import time, datetime
from ortools.math_opt.python import mathopt

p = "/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt"
lines = open(p, encoding="utf-8").read().splitlines()
m = int(lines[0].strip()); L = int(lines[1].strip())
items = []
for line in lines[2:2+m]:
    a = line.split(); items.append((int(a[0]), int(a[1])))
w = []
for l, d in items:
    w += [l]*d
n = len(w)
TOTAL = sum(w)
SLACK_TOTAL = 28*L - TOTAL
K = 28
mm = mathopt.Model()
x = [[mm.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{i}_{b}") for b in range(K)] for i in range(n)]
slack = [mm.add_variable(lb=0.0, ub=SLACK_TOTAL, is_integer=True, name=f"s{b}") for b in range(K)]
for i in range(n):
    mm.add_linear_constraint(mathopt.fast_sum(x[i]) == 1.0, name=f"one{i}")
for b in range(K):
    mm.add_linear_constraint(mathopt.fast_sum([w[i]*x[i][b] for i in range(n)]) + slack[b] == L, name=f"cap{b}")
mm.minimize(mathopt.fast_sum(slack))
t0 = time.time()
res = mathopt.solve(mm, mathopt.SolverType.GSCIP,
                    params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=550), enable_output=False))
print(f"GSCIP K=28: {res.termination.reason.name} | 时间 {round(time.time()-t0,1)}s")
if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
    print(f"目标(总空隙) = {res.objective_value()} | 28 可行: {res.objective_value() <= SLACK_TOTAL}")
elif res.termination.reason == mathopt.TerminationReason.INFEASIBLE:
    print("28 箱不可行（证明）-> 最优 29")
else:
    print("未决")
