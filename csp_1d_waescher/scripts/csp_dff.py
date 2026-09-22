# -*- coding: utf-8 -*-
"""DFF 下界搜索：找对偶可行函数 f 使 Σ f(l_i)*d_i >= 29（证最少 29 辊）。
f 对偶可行: 任意模式 Σ f(l_i)*a_i <= 1 <=> 背包(利润 f) 最优值 <= 1。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import numpy as np, math

p = "/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt"
lines = open(p, encoding="utf-8").read().splitlines()
m = int(lines[0].strip()); L = int(lines[1].strip())
items = []
for line in lines[2:2+m]:
    a = line.split(); items.append((int(a[0]), int(a[1])))
length = np.array([l for l, d in items]); dem = [d for l, d in items]
exp = []
for i, (l, d) in enumerate(items):
    for _ in range(d):
        exp.append((l, i))
n1 = len(exp)
W = np.array([e[0] for e in exp], dtype=np.int64)

def knap_ok(vals_by_type, ub=1.0 + 1e-9):
    """检查 max 模式价值 <= ub（0/1 背包）。"""
    vv = np.array([vals_by_type[e[1]] for e in exp])
    dp = np.zeros(L+1, dtype=np.float64)
    for pos in range(n1):
        w = W[pos]
        cand = dp[:L+1-w] + vv[pos]
        better = cand > dp[w:]
        dp[w:][better] = cand[better]
        if dp.max() > ub:
            return False
    return dp.max() <= ub

best = ("", 0.0)
# 候选 DFF 族
def f_floor(k):
    return [math.floor(k*l/L)/1.0 for l in length]   # f = floor(kx/L)  值域整数，按 1 检查? 应检查 <=1: f=floor(kx/L) 最大可能 k-1 不行
cands = []
for k in range(2, 12):
    # f(x) = floor(k*x/L) / k    (值域 [0,1])
    cands.append((f"floor(kx/L)/k k={k}", [math.floor(k*l/L)/k for l in length]))
    # f(x) = min(1, k*x/L)
    cands.append((f"min(1,kx/L) k={k}", [min(1.0, k*l/L) for l in length]))
    # f(x) = x/L * k/(k+?); f(x)=ceil(kx/L)/k 通常不可行
    cands.append((f"ceil(kx/L)/k k={k}", [math.ceil(k*l/L)/k for l in length]))
# 组合：两个 DFF 的 max/min
import itertools
base = [c for c in cands][:12]
for (n1_, f1), (n2_, f2) in itertools.combinations(base, 2):
    cands.append((f"max({n1_},{n2_})", [max(a, b) for a, b in zip(f1, f2)]))
    cands.append((f"min({n1_},{n2_})", [min(a, b) for a, b in zip(f1, f2)]))

for name, f in cands:
    if knap_ok(f):
        bnd = sum(f[i]*dem[i] for i in range(m))
        if bnd > best[1]:
            best = (name, bnd)
print(f"最佳 DFF 下界: {best[0]} -> {round(best[1], 6)}")
print(f"结论: 下界 {round(best[1],4)} -> 最少辊数 >= {math.ceil(best[1] - 1e-9)}")
