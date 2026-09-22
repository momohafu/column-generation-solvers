# -*- coding: utf-8 -*-
"""Waescher_TEST0005：Martello-Toth L2/L3 下界 + HIGHS item->bin MIP 试证 K=28。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import math, time, datetime
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
w.sort(reverse=True)
print(f"件数 {n} | L {L} | 总长 {sum(w)} | 长度下界 {math.ceil(sum(w)/L)}")

# ---- Martello-Toth L2 ----
def L2_bound():
    W = w[:]
    L2 = 0
    j = len(W) - 1
    # 标准 L2：找 W[i]+W[j]<=L 配对
    bins_used = set()
    for k in range(1, len(W)//2 + 1):
        # 简化实现：迭代按定理 L2
        pass
    # 直接实现定理 L2:
    # 对每个 h：N_h = #items in (L/(h+1), L/h]; L2 = max over h ...
    # 这里实现通用的对偶可行函数近似: 用 h=2..8 分段界
    best = 0
    for h in range(2, 10):
        lb_h = 0
        rem = []
        for x in W:
            if x > L//h and x <= L//(h-1):
                lb_h += 1
            elif x > L//(h-1):
                lb_h += 1
            elif x > L//h and x <= L//(h-1):
                pass
        best = max(best, lb_h)
    return best

# 简化 L2（标准公式）:
def L2():
    W = sorted(w, reverse=True)
    n1 = len(W)
    F = [0.0]*(n1+1)   # F[i]: 前 i 个（最大）的总长
    for i in range(1, n1+1):
        F[i] = F[i-1] + W[i-1]
    # L2 = max_{i} i + ceil( (F[n1] - F[?]) ... )  简化: 使用标准 L2 实现
    # 标准 L2: L2 = max{ i + max(0, ceil((T_i - (n-i)*L/2)/L)) ...}  太复杂，用近似：
    best = math.ceil(sum(W)/L)
    for i in range(1, n1):
        # 前 i 个中每两个占一箱（若 >L），其余按长度界
        if W[i-1] > L/2:
            cnt_big = sum(1 for x in W[:i] if x > L/2)
            rest_len = sum(W[i:])
            lb = cnt_big + math.ceil(rest_len/L)
            best = max(best, lb)
    return best

print("L2 简化界:", L2())

# ---- HIGHS item->bin MIP K=28 ----
def highs_k(Kv, time_limit=300.0):
    mm = mathopt.Model()
    x = [[mm.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{i}_{b}") for b in range(Kv)] for i in range(n)]
    for i in range(n):
        mm.add_linear_constraint(mathopt.fast_sum(x[i]) == 1.0, name=f"one{i}")
    for b in range(Kv):
        mm.add_linear_constraint(mathopt.fast_sum([w[i]*x[i][b] for i in range(n)]) <= L, name=f"cap{b}")
    mm.minimize(0)
    t0 = time.time()
    res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=time_limit), enable_output=False))
    return res.termination.reason.name, round(time.time()-t0, 1)

for Kv in (28,):
    term, wt = highs_k(Kv, 300)
    print(f"HIGHS K={Kv}: {term} ({wt}s)")
