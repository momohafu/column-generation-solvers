# -*- coding: utf-8 -*-
"""探测：固定 K∈[3,3]，分别禁用不同弧，找能产生分数弧流量的情形（触发 Ryan-Foster 弧分支）。"""
import sys, math, time
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg_cpsat import build_data, enumerate_pool
import importlib.util
spec = importlib.util.spec_from_file_location("bp", "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts/branch_and_price.py")
# 避免重复跑 main：手动加载所需函数 —— 直接复制 import 后取 cg_node：
bp = importlib.util.module_from_spec(spec)
import sys as _s
# branch_and_price.py 有 __main__ 吗？没有 guard！它会执行实验……改用直接导入会跑三个实验。
# 因此这里独立实现探测（复用 02 的函数思路），不 import bp。
from cg_cpsat import build_data, enumerate_pool, solve_rmp, pool_scan
from ortools.math_opt.python import mathopt

n, xc, yc, dem, ready, due, svc, cap, depot_due, dist, d_scaled = build_data()
paths, costs, masks, loads = enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due)
path_to_idx = {p: i for i, p in enumerate(paths)}
M_DUMMY = 1e6
EPS = 1e-7
K = 25

def col_arcs(p):
    a = []
    prev = 0
    for j in p:
        a.append((prev, j)); prev = j
    a.append((prev, 0))
    return a

def cg3(forb, K_ub=3, K_lb=3):
    sel = []
    sel_set = set()
    for i in range(1, n+1):
        idx = path_to_idx[(i,)]
        arcs = col_arcs(paths[idx])
        if all(a not in forb for a in arcs):
            sel_set.add(idx); sel.append(idx)
    for it in range(200):
        m = mathopt.Model()
        vs = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in sel]
        yv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(1, n+1)]
        covers = []
        for i in range(1, n+1):
            covers.append(m.add_linear_constraint(
                mathopt.fast_sum([vs[k] for k, p in enumerate(sel) if (masks[p] >> i) & 1]) + yv[i-1] >= 1.0, name=f"c{i}"))
        ubc = m.add_linear_constraint(mathopt.fast_sum(vs) <= K_ub, name="vub")
        lbc = m.add_linear_constraint(mathopt.fast_sum(vs) >= K_lb, name="vlb")
        m.minimize(mathopt.fast_sum([costs[p]*vs[k] for k, p in enumerate(sel)]) + M_DUMMY*mathopt.fast_sum(yv))
        res = mathopt.solve(m, mathopt.SolverType.GLOP)
        dv = res.dual_values()
        pi = [0.0]*(n+1)
        for i in range(1, n+1):
            pi[i] = max(0.0, dv[covers[i-1]])
        mu_sum = dv[ubc] + dv[lbc]
        xv = {p: res.variable_values()[vs[k]] for k, p in enumerate(sel)}
        y_used = sum(res.variable_values()[yv[i-1]] for i in range(1, n+1))
        add = []
        for idx in range(len(paths)):
            if idx in sel_set:
                continue
            arcs = col_arcs(paths[idx])
            if any(a in forb for a in arcs):
                continue
            s = 0.0
            mm = masks[idx]
            while mm:
                lb = mm & -mm
                s += pi[lb.bit_length()-1]
                mm -= lb
            rc = costs[idx] - s - mu_sum
            if rc < -EPS:
                add.append((rc, idx))
        add.sort()
        for _, idx in add[:2000]:
            if idx not in sel_set:
                sel_set.add(idx); sel.append(idx)
        if not add:
            break
    return res.objective_value(), y_used, xv, sel

for arc in [(13,17),(14,12),(15,16),(16,14),(17,18),(4,2),(2,1),(9,6),(6,4),(11,9),(24,25),(23,22),(22,21),(5,3)]:
    obj, y_used, xv, sel = cg3({arc})
    if y_used > 1e-6:
        print(f"禁用 {arc}: LP 不可行（虚拟列使用）")
        continue
    fV = sum(xv.values())
    flows = {}
    for p, val in xv.items():
        if val <= 1e-6: continue
        for a in col_arcs(paths[p]):
            flows[a] = flows.get(a, 0.0) + val
    frac = [(a, f) for a, f in flows.items() if a[0] >= 1 and a[1] >= 1 and 1e-6 < f < 1 - 1e-6]
    print(f"禁用 {arc}: LP={obj:.4f} Σx={fV:.3f} 整数={all(abs(v-round(v))<1e-6 for v in xv.values())} 分数弧数={len(frac)} 最分数={max(frac, key=lambda t: min(t[1],1-t[1])) if frac else None}")
