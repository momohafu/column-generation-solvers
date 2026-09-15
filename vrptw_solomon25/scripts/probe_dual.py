# -*- coding: utf-8 -*-
"""探针：打印第 1 轮 RMP 的对偶值、单客户列 rc、CP-SAT 定价结果（供建模文档引用）。"""
import sys, math
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg_cpsat import (build_data, enumerate_pool, solve_rmp, cpsat_pricing,
                      exact_rc)

n, x, y, dem, ready, due, svc, cap, depot_due, dist, d_scaled = build_data()
K = 25
paths, costs, masks, loads = enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due)
path_to_idx = {p: i for i, p in enumerate(paths)}
selected = [path_to_idx[(i,)] for i in range(1, n + 1)]

obj, pi, mu = solve_rmp(n, K, paths, costs, masks, selected)
print("== 第 1 轮 RMP（初始 25 条单客户列）==")
print("RMP 目标值:", round(obj, 6))
print("mu(车辆数约束对偶):", round(mu, 6))
print("pi(覆盖约束对偶, 客户1..25):")
print([round(v, 4) for v in pi[1:]])
print("校验: 每条单客户列 rc = c_p - pi_i - mu 应 = 0 (绑定列):")
for i in range(1, n + 1):
    p = path_to_idx[(i,)]
    rc = costs[p] - pi[i] - mu
    if i <= 5 or abs(rc) > 1e-9:
        print(f"  客户{i:2d}: c={costs[p]:9.4f}  pi={pi[i]:9.4f}  rc={rc:+.2e}")
print()
path, scaled, st, wt = cpsat_pricing(n, dem, ready, due, svc, cap, depot_due, d_scaled, pi, mu)
rc = exact_rc(path, dist, pi, mu)
d = 0.0
prev = 0
for j in path:
    d += dist(prev, j); prev = j
d += dist(prev, 0)
s = sum(pi[i] for i in path)
print("== CP-SAT 定价子问题（第 1 轮）==")
print("找到路径:", path)
print("路径距离 c_p =", round(d, 4), "| Σ pi =", round(s, 4), "| mu =", round(mu, 4))
print("reduced cost = c - Σpi - mu =", round(rc, 6))
print("缩放目标值(scaled obj):", scaled, "(= rc * 1000 取整)", round(rc*1000, 1))
print("CP-SAT 状态:", st, "| 耗时", round(wt, 3), "s")
