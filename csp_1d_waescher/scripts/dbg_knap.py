# -*- coding: utf-8 -*-
"""最小复现：快照重建为何只取 1 件。"""
import numpy as np
W = np.array([40, 74, 80], dtype=np.int64)   # 3 件
vv = np.array([1.01, 1.02, 1.03])
L = 200
dp = np.zeros(L+1)
snaps = []
for pos in range(3):
    snaps.append(dp.copy())
    w = W[pos]
    cand = dp[:L+1-w] + vv[pos]
    better = cand > dp[w:]
    dp[w:][better] = cand[better]
    print(f"item{pos}(w={w}) 后 dp[200]={dp[200]:.4f}")
best_cap = int(np.argmax(dp))
print("best_cap", best_cap, "val", dp[best_cap])
counts = [0]*3
cap = best_cap
for pos in range(2, -1, -1):
    w = W[pos]
    if cap >= w:
        lhs = snaps[pos][cap-w] + vv[pos]
        print(f"pos={pos} cap={cap} lhs={lhs:.10f} dp[cap]={dp[cap]:.10f} eq={lhs == dp[cap]}")
        if lhs == dp[cap]:
            counts[pos] += 1
            cap -= w
print("counts", counts, "cap", cap)
