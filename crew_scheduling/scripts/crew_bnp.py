# -*- coding: utf-8 -*-
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import math, time, datetime
from ortools.math_opt.python import mathopt
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt"
_lines = open(DATA).read().splitlines()
N, T = map(int, _lines[0].split())
tasks = [None] + [tuple(map(int, l.split())) for l in _lines[1:1+N]]
arc_cost = {}
for l in _lines[1+N:]:
    i, j, c = map(int, l.split()); arc_cost[(i, j)] = c
K_MIN = 27
EPS = 1e-7
M_DUMMY = 10**6

def enumerate_pool():
    adj = {}
    for (i, j), c in arc_cost.items():
        adj.setdefault(i, []).append((j, c))
    paths = []
    for start in range(1, N+1):
        s0 = tasks[start][0]
        stack = [(start, [start], 0)]
        while stack:
            u, seq, c = stack.pop()
            paths.append((tuple(seq), c, tasks[u][1]-s0))
            for v, cv in adj.get(u, []):
                span = tasks[v][1] - s0
                if span <= T:
                    stack.append((v, seq+[v], c+cv))
    return paths

paths = enumerate_pool()
P = len(paths)
pseq = [p[0] for p in paths]
pcost = [p[1] for p in paths]
pmask = []
for s2 in pseq:
    m2 = 0
    for i in s2:
        m2 |= (1 << i)
    pmask.append(m2)

def col_arcs(seq):
    a = []
    prev = 0
    for j in seq:
        a.append((prev, j))
        prev = j
    a.append((prev, N+1))
    return a

def solve_rmp_node(selected, K_ub, K_lb, forced, forbidden):
    m = mathopt.Model()
    x = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in selected]
    yv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(1, N+1)]
    covers = []
    for i in range(1, N+1):
        covers.append(m.add_linear_constraint(
            mathopt.fast_sum([x[k] for k, p in enumerate(selected) if (pmask[p] >> i) & 1]) + yv[i-1] >= 1.0, name=f"c{i}"))
    ub_c = m.add_linear_constraint(mathopt.fast_sum(x) <= K_ub, name="vub")
    lb_c = m.add_linear_constraint(mathopt.fast_sum(x) >= K_lb, name="vlb")
    m.minimize(mathopt.fast_sum([pcost[p]*x[k] for k, p in enumerate(selected)]) + M_DUMMY*mathopt.fast_sum(yv))
    res = mathopt.solve(m, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
    dv = res.dual_values()
    pi = [0.0]*(N+1)
    for i in range(1, N+1):
        pi[i] = max(0.0, dv[covers[i-1]])
    mu_ub = dv[ub_c]
    mu_lb = dv[lb_c]
    xvals = {p: res.variable_values()[x[k]] for k, p in enumerate(selected)}
    y_used = sum(res.variable_values()[yv[i-1]] for i in range(1, N+1))
    return res.objective_value(), pi, mu_ub, mu_lb, xvals, y_used

def eligible(p, forced, forbidden):
    for arc in col_arcs(pseq[p]):
        if arc in forbidden:
            return False
    for arc in forced:
        if arc not in col_arcs(pseq[p]):
            return False
    return True

def cg_node(K_ub, K_lb, forced, forbidden, time_limit=60.0):
    forced = tuple(forced)
    forbidden = set(forbidden)
    t0 = time.time()
    sel = [p for p in range(P) if len(pseq[p]) == 1 and eligible(p, forced, forbidden)]
    sel_set = set(sel)
    if forced:
        cnt = 0
        for p in range(P):
            if p in sel_set:
                continue
            if eligible(p, forced, forbidden):
                sel_set.add(p)
                sel.append(p)
                cnt += 1
                if cnt >= 50:
                    break
    if not sel:
        return dict(lp=None, infeasible=True, integral=False, xvals={}, iters=0)
    lp = None
    xvals = {}
    y_used = 0.0
    for it in range(200):
        obj, pi, mu_ub, mu_lb, xvals, y_used = solve_rmp_node(sel, K_ub, K_lb, forced, forbidden)
        lp = obj
        mu_sum = mu_ub + mu_lb
        add = []
        for p in range(P):
            if p in sel_set or not eligible(p, forced, forbidden):
                continue
            s2 = 0.0
            mm = pmask[p]
            while mm:
                lb = mm & -mm
                s2 += pi[lb.bit_length()-1]
                mm -= lb
            rc = pcost[p] - s2 - mu_sum
            if rc < -EPS:
                add.append((rc, p))
        add.sort()
        add = [p for _, p in add[:1000]]
        if not add:
            break
        for p in add:
            if p not in sel_set:
                sel_set.add(p)
                sel.append(p)
    infeasible = (y_used > 1e-6)
    integral = (not infeasible) and all(abs(v-round(v)) < 1e-6 for v in xvals.values())
    return dict(lp=None if infeasible else lp, infeasible=infeasible, integral=integral,
                xvals=xvals, iters=it+1, time=time.time()-t0)

def run_bnp(extra_forbidden=(), root_K_ub=27, root_K_lb=27, time_limit=110.0):
    wall0 = time.time()
    incumbent = None
    inc_routes = None
    root = dict(K_ub=root_K_ub, K_lb=root_K_lb, forced=(), forbidden=set(extra_forbidden))
    stack = [root]
    nodes = 0
    branches = 0
    while stack and time.time()-wall0 < time_limit:
        node = stack.pop()
        nodes += 1
        res = cg_node(node["K_ub"], node["K_lb"], node["forced"], node["forbidden"])
        tag = f"node{nodes} K[{node['K_lb']},{node['K_ub']}]"
        if node["forced"]:
            tag += f" 强制{list(node['forced'])}"
        if node["forbidden"]:
            tag += f" 禁用{len(node['forbidden'])}弧"
        if res["infeasible"]:
            print(f"{tag}: 不可行 -> 剪枝 | CG {res['iters']} 轮 {res['time']:.1f}s")
            continue
        v = res["lp"]
        if incumbent is not None and v >= incumbent - 1e-7:
            print(f"{tag}: LP={v:.4f} >= incumbent {incumbent:.4f} -> 定界剪枝")
            continue
        if res["integral"]:
            routes = [pseq[p] for p, val in res["xvals"].items() if val > 0.5]
            if incumbent is None or v < incumbent:
                incumbent = v
                inc_routes = routes
            print(f"{tag}: LP={v:.4f} 整数解 ({len(routes)} crew) -> 节点解决, 更新 incumbent")
            continue
        fV = sum(res["xvals"].values())
        if abs(fV - round(fV)) > 1e-6:
            branches += 1
            stack.append(dict(node, K_lb=int(math.ceil(fV))))
            stack.append(dict(node, K_ub=int(math.floor(fV))))
            print(f"{tag}: LP={v:.4f} 分数(Σx={fV:.3f}) -> 车辆数分支")
            continue
        flows = {}
        for p, val in res["xvals"].items():
            if val <= 1e-6:
                continue
            for arc in col_arcs(pseq[p]):
                flows[arc] = flows.get(arc, 0.0) + val
        cands = [a for a, f in flows.items() if a[0] >= 1 and a[1] >= 1 and a[1] <= N and 1e-6 < f < 1-1e-6]
        if not cands:
            print(f"{tag}: 无分数弧 -> 停止")
            break
        arc = max(cands, key=lambda a: min(flows[a], 1-flows[a]))
        branches += 1
        stack.append(dict(node, forced=node["forced"] + (arc,)))
        stack.append(dict(node, forbidden=node["forbidden"] | {arc}))
        print(f"{tag}: LP={v:.4f} 分数 -> 弧分支 {arc} (流量 {flows[arc]:.3f})")
    print("== B&P 汇总 ==")
    print("节点:", nodes, "| 分支:", branches, "| 墙钟:", round(time.time()-wall0, 2), "s")
    print("incumbent:", None if incumbent is None else round(incumbent, 4))
    if inc_routes:
        for r in sorted(inc_routes, key=lambda s: -len(s)):
            print("  crew:", r, "成本", sum(arc_cost.get((r[k], r[k+1]), 0) for k in range(len(r)-1)))
    return dict(incumbent=incumbent, routes=inc_routes, nodes=nodes, branches=branches)

if __name__ == "__main__":
    print("== 实验一：原始 csp50 ==")
    run_bnp()
    print()
    print("== 实验二：分支验证（禁用最优弧 (1,10)）==")
    run_bnp(extra_forbidden={(1, 10)})
