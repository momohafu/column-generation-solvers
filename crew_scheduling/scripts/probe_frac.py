# -*- coding: utf-8 -*-
"""探测：K=[27,27] 下禁用单个弧，找节点 LP 分数（触发 Ryan-Foster 弧分支）的情形。"""
import sys
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts")
from crew_bnp import cg_node, arc_cost
found = 0
for (i, j) in sorted(arc_cost):
    res = cg_node(27, 27, (), {(i, j)}, time_limit=20.0)
    if res["lp"] is None:
        print(f"禁用 ({i},{j}): 节点不可行")
        continue
    if not res["integral"]:
        # 找最分数弧
        flows = {}
        from crew_bnp import col_arcs, pseq
        for p, val in res["xvals"].items():
            if val <= 1e-6:
                continue
            for arc in col_arcs(pseq[p]):
                flows[arc] = flows.get(arc, 0.0) + val
        frac = [(a, f) for a, f in flows.items() if a[0] >= 1 and a[1] >= 1 and a[1] <= 50 and 1e-6 < f < 1-1e-6]
        best = max(frac, key=lambda t: min(t[1], 1-t[1])) if frac else None
        print(f"禁用 ({i},{j}): LP={res['lp']:.2f} 分数 | Σx={sum(res['xvals'].values()):.3f} | 最分数弧 {best}")
        found += 1
        if found >= 4:
            break
print("探测完成, 分数情形数(前4):", found)
