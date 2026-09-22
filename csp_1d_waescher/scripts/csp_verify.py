# -*- coding: utf-8 -*-
"""验证 CG 池 MIP 的 28 辊解 + LBBD 惰性割 K=28 测试。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, math, time, datetime
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/csp_1d_waescher/scripts")
from csp_core import (ITEMS, LENGTH, DEM, L, TOTAL, M, N1, cg_min_rolls, ip_recovery,
                      knap_rebuild)
from ortools.math_opt.python import mathopt

# 1) 28 解验证
lp, patterns, sel, iters = cg_min_rolls()
# 重建整数解：直接解池 MIP 并取 x 值
mm = mathopt.Model()
x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=True, name=f"x{k}") for k in sel]
for i in range(M):
    mm.add_linear_constraint(mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= DEM[i], name=f"c{i}")
mm.minimize(mathopt.fast_sum(x))
res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                    params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
Ksol = res.objective_value()
print(f"整数解: {Ksol} 辊 | {res.termination.reason.name}")
used = [(k, patterns[k]) for k in sel if res.variable_values()[x[k]] > 0.5]
print(f"使用模式数: {len(used)}")
total_check = 0
cnt = [0]*M
ok = True
for k, pat in used:
    plen = sum(pat[i]*LENGTH[i] for i in range(M))
    assert plen <= L, ("超长", plen)
    for i in range(M):
        cnt[i] += pat[i]
for i in range(M):
    if cnt[i] != DEM[i]:
        ok = False
        print(f"  需求不符 类型{i}: {cnt[i]} vs {DEM[i]}")
print("覆盖校验: 每类型计数=需求 ->", ok)
print("总件数:", sum(cnt), "= 114 ->", sum(cnt) == N1)
print("总长度:", sum(cnt[i]*LENGTH[i] for i in range(M)), "=", TOTAL)
print(f"最优性: LB={math.ceil(TOTAL/L)} = UB={Ksol} -> 最优 {Ksol} 已证明:",
      math.ceil(TOTAL/L) == Ksol)

# 2) LBBD 惰性容量割 K=28（带紧负载界）
items = sorted([l for l, d in ITEMS for _ in range(d)], reverse=True)
n = len(items)
K = 28
SLACK = K*L - TOTAL
cuts = []
t0 = time.time()
feasible = False
for it in range(30):
    mm2 = mathopt.Model()
    x2 = [[mm2.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{i}_{b}") for b in range(K)] for i in range(n)]
    for i in range(n):
        mm2.add_linear_constraint(mathopt.fast_sum(x2[i]) == 1.0, name=f"one{i}")
    for b in range(K):
        mm2.add_linear_constraint(mathopt.fast_sum([items[i]*x2[i][b] for i in range(n)]) <= L, name=f"cap{b}")
        mm2.add_linear_constraint(mathopt.fast_sum([items[i]*x2[i][b] for i in range(n)]) >= L - SLACK, name=f"lb{b}")
    for (S, b) in cuts:
        mm2.add_linear_constraint(mathopt.fast_sum([x2[i][b] for i in S]) <= len(S)-1, name=f"cut")
    mm2.minimize(0)
    res2 = mathopt.solve(mm2, mathopt.SolverType.HIGHS,
                         params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=90), enable_output=False))
    if res2.termination.reason == mathopt.TerminationReason.INFEASIBLE:
        print(f"LBBD iter {it+1}: 主问题不可行（K=28 不可行）？")
        break
    if res2.termination.reason not in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
        print(f"LBBD iter {it+1}: {res2.termination.reason.name} 未决")
        break
    xv = res2.variable_values()
    added = 0
    for b in range(K):
        S = [i for i in range(n) if xv[x2[i][b]] > 0.5]
        if sum(items[i] for i in S) > L:
            S2 = S[:]
            for i in S:
                if len(S2) <= 1:
                    break
                S3 = [j for j in S2 if j != i]
                if sum(items[j] for j in S3) > L:
                    S2 = S3
            cuts.append((tuple(S2), b))
            added += 1
    print(f"LBBD iter {it+1}: 超载箱 {added}（累计割 {len(cuts)}）")
    if added == 0:
        feasible = True
        print(f"LBBD: 28 辊可行分配找到！(iter {it+1})")
        break
    if time.time()-t0 > 150:
        break
print("LBBD 总时间:", round(time.time()-t0, 1), "s | 可行:", feasible)
