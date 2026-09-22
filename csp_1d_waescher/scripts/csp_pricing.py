# -*- coding: utf-8 -*-
"""numpy 向量化 0/1 背包定价（含重建）+ 基准。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import numpy as np, time

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
print(f"件数 {n1} | L {L}")

def knapsack(pi, exclude_item=None):
    """numpy 0/1 背包：max Σ pi * a s.t. Σ l*a <= L。exclude_item: 排除该类型的全部件。
    返回 (最优值, 各类型计数)。快照式重建保证正确性。"""
    vals = np.array([pi[i] for i in range(m)], dtype=np.float64)
    mask = np.ones(n1, dtype=bool)
    if exclude_item is not None:
        mask &= np.array([e[1] != exclude_item for e in exp])
    ww = W[mask]
    vv = np.array([vals[e[1]] for e in exp])[mask]
    order = np.arange(n1)[mask]
    dp = np.zeros(L+1, dtype=np.float64)
    snaps = []
    for pos, idx in enumerate(order):
        snaps.append(dp.copy())
        w = ww[pos]
        cand = dp[:L+1-w] + vv[pos]
        better = cand > dp[w:]
        dp[w:][better] = cand[better]
    best_cap = int(np.argmax(dp))
    best_val = dp[best_cap]
    counts = [0]*m
    cap = best_cap
    cur = dp                               # 回溯时随步切换"当前"dp 状态
    for pos in range(len(order)-1, -1, -1):
        w = ww[pos]
        idx = order[pos]
        if cap >= w and cur[cap] == snaps[pos][cap-w] + vv[pos]:
            counts[exp[idx][1]] += 1
            cur = snaps[pos]
            cap -= w
    return best_val, counts

# 基准：对一组对偶试算
pi = [1.0 + 0.01*i for i in range(m)]
t0 = time.time()
for _ in range(5):
    v, c = knapsack(pi)
dt = time.time()-t0
print(f"5 次定价耗时 {round(dt,3)}s（单次 {round(dt/5*1000,1)}ms）| 价值 {round(v,4)} | 件数 {sum(c)} | 长度 {sum(c[i]*length[i] for i in range(m))}")
# forbid-pair 定价正确性：不含 (i,j) 的最大价值 = max(排除 i, 排除 j)
i, j = 0, 1
vi, ci = knapsack(pi, exclude_item=i)
vj, cj = knapsack(pi, exclude_item=j)
vfull, cfull = knapsack(pi)
print(f"forbid({i},{j}): max(排除i={round(vi,3)}, 排除j={round(vj,3)}) = {round(max(vi,vj),3)} vs 全局 {round(vfull,3)}")
