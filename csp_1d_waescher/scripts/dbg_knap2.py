# -*- coding: utf-8 -*-
"""真实数据上的重建调试。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import numpy as np

p = "/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt"
lines = open(p, encoding="utf-8").read().splitlines()
m = int(lines[0].strip()); L = int(lines[1].strip())
items = []
for line in lines[2:2+m]:
    a = line.split(); items.append((int(a[0]), int(a[1])))
length = [l for l, d in items]; dem = [d for l, d in items]
exp = []
for i, (l, d) in enumerate(items):
    for _ in range(d):
        exp.append((l, i))
n1 = len(exp)
W = np.array([e[0] for e in exp], dtype=np.int64)
pi = [1.0 + 0.01*i for i in range(m)]
vals = np.array([pi[i] for i in range(m)], dtype=np.float64)
vv = np.array([vals[e[1]] for e in exp])
mask = np.ones(n1, dtype=bool)
ww = W[mask]
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
print("best_cap", best_cap, "val", dp[best_cap], "nitems", n1)
counts = [0]*m
cap = best_cap
steps = 0
for pos in range(len(order)-1, -1, -1):
    w = ww[pos]
    if cap >= w:
        lhs = snaps[pos][cap-w] + vv[pos]
        if lhs == dp[cap]:
            counts[exp[order[pos]][1]] += 1
            cap -= w
            steps += 1
        elif steps < 3:
            print(f"  失败 pos={pos} cap={cap} w={w} lhs={lhs:.12f} dp={dp[cap]:.12f}")
    if steps == 0 and pos < len(order)-3:
        pass
print("steps", steps, "cap", cap, "总件数", sum(counts), "总长", sum(counts[i]*length[i] for i in range(m)))
