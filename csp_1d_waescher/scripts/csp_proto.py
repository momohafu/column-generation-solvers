# -*- coding: utf-8 -*-
"""Waescher_TEST0005 基线：总长度下界 + 完整 CG(LP) + CP-SAT 直接模型。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import math, time, datetime
from ortools.math_opt.python import mathopt
from ortools.sat.python import cp_model

p = "/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt"
lines = open(p, encoding="utf-8").read().splitlines()
m = int(lines[0].strip()); L = int(lines[1].strip())
items = []
for line in lines[2:2+m]:
    a = line.split(); items.append((int(a[0]), int(a[1])))
length = [l for l, d in items]; dem = [d for l, d in items]
total_len = sum(l*d for l, d in items)
LB = math.ceil(total_len / L)
print(f"m={m} L={L} | 总长度 {total_len} | 长度下界 LB=ceil({total_len}/{L})={LB}")

# ---------- 列生成（Gilmore-Gomory）：RMP min Σx，定价 = 背包 DP ----------
exp = []   # 展开的 0/1 件（每单位一件）
for i, (l, d) in enumerate(items):
    for _ in range(d):
        exp.append((l, i))
n1 = len(exp)
print("展开件数 n =", n1)

def knapsack(pi):
    """定价子问题：max Σ pi_i a_i s.t. Σ l_i a_i <= L（0/1 背包 DP，含重建）。"""
    dp = [-1e300]*(L+1)
    dp[0] = 0.0
    keep = [bytearray(L+1) for _ in range(n1)]
    for idx, (w, typ) in enumerate(exp):
        val = pi[typ]
        k = keep[idx]
        for cap in range(L, w-1, -1):
            if dp[cap-w] > -1e299:
                nv = dp[cap-w] + val
                if nv > dp[cap]:
                    dp[cap] = nv
                    k[cap] = 1
    best_cap = max(range(L+1), key=lambda c: dp[c])
    best_val = dp[best_cap]
    pat = [0]*m
    cap = best_cap
    for idx in range(n1-1, -1, -1):
        if keep[idx][cap]:
            pat[exp[idx][1]] += 1
            cap -= exp[idx][0]
    return best_val, pat

# 初始模式：单件最大复制
patterns = []
for i, (l, d) in enumerate(items):
    a = [0]*m
    a[i] = min(d, L//l)
    patterns.append(a)

sel = list(range(len(patterns)))
t0 = time.time()
for it in range(200):
    # RMP: min Σx s.t. Σ a_ip x_p >= d_i
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in sel]
    covers = []
    for i in range(m):
        covers.append(mm.add_linear_constraint(
            mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= dem[i], name=f"c{i}"))
    mm.minimize(mathopt.fast_sum(x))
    res = mathopt.solve(mm, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
    lp = res.objective_value()
    dv = res.dual_values()
    pi = [max(0.0, dv[covers[i]]) for i in range(m)]
    best_val, pat = knapsack(pi)
    rc = 1.0 - best_val
    if rc >= -1e-7:
        print(f"CG 收敛: {it+1} 轮 | LP 下界(辊数) = {round(lp, 6)} | 模式 {len(sel)} | 时间 {round(time.time()-t0,1)}s")
        break
    patterns.append(pat)
    sel.append(len(patterns)-1)
    if it % 20 == 0 or it < 5:
        print(f"iter {it+1}: LP={round(lp,6)} | 定价值 {round(best_val,4)} rc={round(rc,4)}")

# ---------- 整数恢复：生成模式池 MIP ----------
mip = mathopt.Model()
x = [mip.add_variable(lb=0.0, ub=float("inf"), is_integer=True, name=f"x{k}") for k in sel]
for i in range(m):
    mip.add_linear_constraint(mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= dem[i], name=f"c{i}")
mip.minimize(mathopt.fast_sum(x))
t1 = time.time()
res = mathopt.solve(mip, mathopt.SolverType.HIGHS,
                    params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
print(f"整数恢复: {res.termination.reason.name} | 辊数 {res.objective_value()} | 时间 {round(time.time()-t1,1)}s")
print(f"LP 下界 {round(lp,4)} vs 整数 {res.objective_value()} | gap {round((res.objective_value()-lp)/lp*100,4)}%")

# ---------- 直接 CP-SAT（item->bin 分配）----------
def cpsat_bins(K, time_limit=120.0):
    model = cp_model.CpModel()
    a = [[model.NewIntVar(0, dem[i], f"a{i}_{k}") for k in range(K)] for i in range(m)]
    for i in range(m):
        model.Add(sum(a[i][k] for k in range(K)) == dem[i])
    for k in range(K):
        model.Add(sum(a[i][k]*length[i] for i in range(m)) <= L)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    st = solver.Solve(model)
    return st, solver.WallTime()

for K in (28, 29):
    t1 = time.time()
    st, wt = cpsat_bins(K, 180)
    print(f"CP-SAT K={K}: {st} | 墙钟 {round(time.time()-t1,1)}s")
    if st == cp_model.OPTIMAL:
        print(f"  -> K={K} 可行最优，最少辊数 = {K}")
        break
    if st == cp_model.INFEASIBLE:
        print(f"  -> K={K} 不可行，最少辊数 >= {K+1}")
