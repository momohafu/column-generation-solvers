# -*- coding: utf-8 -*-
"""K=28 攻击 v2：每箱负载 ∈ [9935,10000]（总空隙 65 的必然推论）——强传播。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import time
from ortools.sat.python import cp_model

p = "/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt"
lines = open(p, encoding="utf-8").read().splitlines()
m = int(lines[0].strip()); L = int(lines[1].strip())
items = []
for line in lines[2:2+m]:
    a = line.split(); items.append((int(a[0]), int(a[1])))
w = []
for l, d in items:
    w += [l]*d
w.sort(reverse=True)
n = len(w)
TOTAL = sum(w)
SLACK_TOTAL = 28*L - TOTAL
print(f"件数 {n} | 总空隙 {SLACK_TOTAL} -> 每箱负载 ∈ [{L-SLACK_TOTAL}, {L}]")

K = 28
model = cp_model.CpModel()
x = [[model.NewBoolVar(f"x{i}_{b}") for b in range(K)] for i in range(n)]
load = [model.NewIntVar(L-SLACK_TOTAL, L, f"ld{b}") for b in range(K)]
for i in range(n):
    model.Add(sum(x[i]) == 1)
for b in range(K):
    model.Add(sum(w[i]*x[i][b] for i in range(n)) == load[b])
for b in range(K-1):
    model.Add(load[b] >= load[b+1])
model.Minimize(K*L - sum(load))
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = 560
solver.parameters.num_search_workers = 8
t0 = time.time()
st = solver.Solve(model)
el = time.time()-t0
print(f"状态: {st} | 时间 {round(el,1)}s | 最优空隙: {solver.ObjectiveValue() if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None}")
if st == cp_model.OPTIMAL:
    print("结论: 28 箱可行!" if solver.ObjectiveValue() <= SLACK_TOTAL else "结论: 28 箱不可行 -> 最优 29")
elif st == cp_model.INFEASIBLE:
    print("结论: 28 箱不可行（证明）-> 最优 29")
else:
    print("未决")
