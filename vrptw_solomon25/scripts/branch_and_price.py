# -*- coding: utf-8 -*-
"""VRPTW c101 branch-and-price 原型：
节点 = 受限主问题(覆盖 LP，含车辆数上下界 + 强制/禁用弧) 列生成求解（CP-SAT 定价 + 池扫描精确校验）
分支 = 车辆数分支（Σx 分数）或 Ryan-Foster 弧分支（流量最接近 0.5 的客户弧：禁用 vs 强制）
剪枝 = LP 下界 ≥ incumbent；节点 LP 整数（x∈{0,1}）或不可行
"""
import sys, math, time, datetime
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg_cpsat import build_data, enumerate_pool
from ortools.sat.python import cp_model
from ortools.math_opt.python import mathopt

n, xc, yc, dem, ready, due, svc, cap, depot_due, dist, d_scaled = build_data()
K = 25
LB_K = 3                      # ceil(460/200)
EPS = 1e-7
M_DUMMY = 1e6
SCALE = 1000
DUAL_SCALE = 1000
BIG = 10**9

t0 = time.time()
paths, costs, masks, loads = enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due)
path_to_idx = {p: i for i, p in enumerate(paths)}
print("完整池:", len(paths), "列, 枚举耗时", round(time.time()-t0, 1), "s")

def col_arcs(p):
    a = []
    prev = 0
    for j in p:
        a.append((prev, j))
        prev = j
    a.append((prev, 0))
    return a

def eligible(idx, forced, forbidden):
    for arc in col_arcs(paths[idx]):
        if arc in forbidden:
            return False
    for arc in forced:
        if arc not in col_arcs(paths[idx]):
            return False
    return True

def solve_rmp_node(selected, K_ub, K_lb):
    """RMP：覆盖 LP + 车辆数上下界 + 虚拟列。返回 (obj, pi, mu_ub, mu_lb, xvals, y_used)。"""
    m = mathopt.Model(name="rmp")
    vs = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in selected]
    yv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(1, n+1)]
    covers = []
    for i in range(1, n+1):
        covers.append(m.add_linear_constraint(
            mathopt.fast_sum([vs[k] for k, p in enumerate(selected) if (masks[p] >> i) & 1]) + yv[i-1] >= 1.0, name=f"c{i}"))
    ub_con = m.add_linear_constraint(mathopt.fast_sum(vs) <= K_ub, name="vub")
    lb_con = m.add_linear_constraint(mathopt.fast_sum(vs) >= K_lb, name="vlb")
    m.minimize(mathopt.fast_sum([costs[p]*vs[k] for k, p in enumerate(selected)]) + M_DUMMY*mathopt.fast_sum(yv))
    res = mathopt.solve(m, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL, res.termination.reason
    dv = res.dual_values()
    pi = [0.0]*(n+1)
    for i in range(1, n+1):
        pi[i] = max(0.0, dv[covers[i-1]])
    mu_ub = dv[ub_con]
    mu_lb = dv[lb_con]
    xvals = {p: res.variable_values()[vs[k]] for k, p in enumerate(selected)}
    y_used = sum(res.variable_values()[yv[i-1]] for i in range(1, n+1))
    return res.objective_value(), pi, mu_ub, mu_lb, xvals, y_used

def cpsat_pricing(pi, mu_sum, forced, forbidden, time_limit=10.0):
    model = cp_model.CpModel()
    xv = [[model.NewBoolVar(f"x{i}_{j}") for j in range(n+1)] for i in range(n+1)]
    model.Add(xv[0][0] == 0)
    model.AddCircuit([(i, j, xv[i][j]) for i in range(n+1) for j in range(n+1)])
    visited = [model.NewBoolVar(f"v{i}") for i in range(n+1)]
    model.Add(visited[0] == 1)
    for i in range(1, n+1):
        model.Add(visited[i] + xv[i][i] == 1)
    model.Add(sum(visited) >= 2)
    for (i, j) in forced:
        model.Add(xv[i][j] == 1)
    for (i, j) in forbidden:
        model.Add(xv[i][j] == 0)
    T = depot_due * SCALE
    t = [model.NewIntVar(0, T, f"t{i}") for i in range(n+1)]
    model.Add(t[0] == 0)
    for i in range(1, n+1):
        model.Add(t[i] >= ready[i]*SCALE - BIG*xv[i][i])
        model.Add(t[i] <= due[i]*SCALE + BIG*xv[i][i])
        model.Add(t[i] <= BIG*(1 - xv[i][i]))
    for i in range(n+1):
        for j in range(1, n+1):
            if i == j:
                continue
            model.Add(t[j] >= t[i] + svc[i]*SCALE + d_scaled[i][j] - BIG*(1 - xv[i][j]))
    for i in range(1, n+1):
        model.Add(t[i] + svc[i]*SCALE + d_scaled[i][0] <= depot_due*SCALE + BIG*(1 - xv[i][0]))
    q = [model.NewIntVar(0, cap, f"q{i}") for i in range(n+1)]
    model.Add(q[0] == 0)
    for i in range(n+1):
        for j in range(1, n+1):
            if i == j:
                continue
            model.Add(q[j] >= q[i] + dem[j] - BIG*(1 - xv[i][j]))
    pi_s = [int(round(pi[i]*DUAL_SCALE)) for i in range(n+1)]
    mu_s = int(round(mu_sum*DUAL_SCALE))
    model.Minimize(sum(d_scaled[i][j]*xv[i][j] for i in range(n+1) for j in range(n+1) if i != j)
                   - sum(pi_s[i]*visited[i] for i in range(1, n+1)) - mu_s)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    st = solver.Solve(model)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None
    path = []
    cur = 0
    while True:
        nxt = None
        for j in range(n+1):
            if solver.Value(xv[cur][j]):
                nxt = j
                break
        if nxt is None or nxt == 0:
            break
        path.append(nxt)
        cur = nxt
        if len(path) > n:
            break
    return tuple(path)

def exact_rc(path, pi, mu_sum):
    c = 0.0
    prev = 0
    for j in path:
        c += dist(prev, j)
        prev = j
    c += dist(prev, 0)
    return c - sum(pi[i] for i in path) - mu_sum

def cg_node(K_ub, K_lb, forced, forbidden, time_limit=60.0, verbose=False):
    """节点列生成。返回 dict。"""
    forced = tuple(forced)
    forbidden = set(forbidden)
    t0 = time.time()
    selected = []
    sel_set = set()
    for i in range(1, n+1):
        idx = path_to_idx[(i,)]
        if eligible(idx, forced, forbidden):
            sel_set.add(idx)
            selected.append(idx)
    if forced:
        # 强制弧节点：单客户列不含客户-客户弧，须用含强制弧的池列做种子，保证 Σx>=K_lb 可行
        cnt = 0
        for idx in range(len(paths)):
            if idx in sel_set:
                continue
            if eligible(idx, forced, forbidden):
                sel_set.add(idx)
                selected.append(idx)
                cnt += 1
                if cnt >= 50:
                    break
    if not selected:
        return dict(lp_obj=None, infeasible=True, integral=False, iters=0,
                    ncols=0, xvals={}, time=0.0)
    it = 0
    lp_obj = None
    xvals = {}
    y_used = 0.0
    while it < 200 and time.time()-t0 < time_limit:
        it += 1
        obj, pi, mu_ub, mu_lb, xvals, y_used = solve_rmp_node(selected, K_ub, K_lb)
        lp_obj = obj
        mu_sum = mu_ub + mu_lb
        cpath = cpsat_pricing(pi, mu_sum, forced, forbidden)
        cpsat_rc = exact_rc(cpath, pi, mu_sum) if cpath is not None else None
        add = []
        for idx in range(len(paths)):
            if idx in sel_set:
                continue
            if not eligible(idx, forced, forbidden):
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
        add = [i for _, i in add[:2000]]
        for idx in add:
            if idx not in sel_set:
                sel_set.add(idx)
                selected.append(idx)
        if not add and (cpsat_rc is None or cpsat_rc >= -EPS):
            break
    infeasible = (y_used > 1e-6)
    integral = (not infeasible) and all(abs(v - round(v)) < 1e-6 for v in xvals.values())
    return dict(lp_obj=None if infeasible else lp_obj, infeasible=infeasible, integral=integral,
                iters=it, ncols=len(selected), xvals=xvals, time=time.time()-t0)

def branch_and_price(extra_forbidden=(), root_K_ub=None, root_K_lb=None, time_limit=110.0):
    wall0 = time.time()
    incumbent = None
    inc_routes = None
    root = dict(K_ub=root_K_ub if root_K_ub is not None else K,
                K_lb=root_K_lb if root_K_lb is not None else LB_K,
                forced=(), forbidden=set(extra_forbidden))
    stack = [root]
    nodes = 0
    branches = 0
    log = []
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
            print(f"{tag}: 节点不可行 -> 剪枝 | CG {res['iters']} 轮 {res['time']:.1f}s")
            log.append((tag, "infeasible", None, res["iters"]))
            continue
        v = res["lp_obj"]
        if incumbent is not None and v >= incumbent - 1e-7:
            print(f"{tag}: LP={v:.6f} >= incumbent {incumbent:.6f} -> 定界剪枝 | CG {res['iters']} 轮 {res['time']:.1f}s")
            log.append((tag, "bound", v, res["iters"]))
            continue
        if res["integral"]:
            routes = [paths[p] for p, val in res["xvals"].items() if val > 0.5]
            if incumbent is None or v < incumbent:
                incumbent = v
                inc_routes = routes
            print(f"{tag}: LP={v:.6f} 整数解 ({len(routes)} 车) -> 节点解决, 更新 incumbent | CG {res['iters']} 轮 {res['time']:.1f}s")
            log.append((tag, "integral", v, res["iters"]))
            continue
        # ---- 分支 ----
        fV = sum(res["xvals"].values())
        if abs(fV - round(fV)) > 1e-6:
            branches += 1
            c1 = dict(node, K_ub=int(math.floor(fV)))
            c2 = dict(node, K_lb=int(math.ceil(fV)))
            stack.append(c2)
            stack.append(c1)
            print(f"{tag}: LP={v:.6f} 分数(Σx={fV:.3f}) -> 车辆数分支 [{node['K_lb']},{int(math.floor(fV))}] vs [{int(math.ceil(fV))},{node['K_ub']}]")
            log.append((tag, "branch_veh", v, res["iters"]))
            continue
        flows = {}
        for p, val in res["xvals"].items():
            if val <= 1e-6:
                continue
            for arc in col_arcs(paths[p]):
                flows[arc] = flows.get(arc, 0.0) + val
        cands = [a for a, f in flows.items() if a[0] >= 1 and a[1] >= 1 and 1e-6 < f < 1 - 1e-6]
        if not cands:
            print(f"{tag}: LP={v:.6f} 无分数弧（数值异常）-> 停止")
            log.append((tag, "stop", v, res["iters"]))
            break
        arc = max(cands, key=lambda a: min(flows[a], 1 - flows[a]))
        branches += 1
        c_forb = dict(node, forbidden=node["forbidden"] | {arc})
        c_force = dict(node, forced=node["forced"] + (arc,))
        stack.append(c_force)
        stack.append(c_forb)
        print(f"{tag}: LP={v:.6f} 分数 -> 弧分支 {arc} (流量 {flows[arc]:.3f}): 禁用 vs 强制")
        log.append((tag, f"branch_arc {arc}", v, res["iters"]))
    wall = time.time() - wall0
    print()
    print("== B&P 汇总 ==")
    print("节点数:", nodes, "| 分支数:", branches, "| 墙钟:", round(wall, 2), "s")
    print("incumbent:", None if incumbent is None else round(incumbent, 10),
          "（两位小数", None if incumbent is None else round(incumbent, 2), "）")
    if inc_routes:
        for r in inc_routes:
            seq = [0] + list(r) + [0]
            c = sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1))
            print("  路线", r, "距离", round(c, 6))
    return dict(incumbent=incumbent, routes=inc_routes, nodes=nodes, branches=branches, wall=wall, log=log)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "exp1"):
        print("========== 实验一：原始 c101（完整 B&P）==========")
        r1 = branch_and_price()
        print("BKS 对比: 3 车 / 191.81 -> match:",
              r1["incumbent"] is not None and abs(round(r1["incumbent"], 2) - 191.81) < 1e-9)
    if which in ("all", "exp2"):
        print()
        print("========== 实验二：分支机制验证（临时禁用弧 13->17）==========")
        r2 = branch_and_price(extra_forbidden={(13, 17)})
    if which in ("all", "exp3"):
        print()
        print("========== 实验三：Ryan-Foster 弧分支验证（固定 3 车 + 禁用 15->16）==========")
        r3 = branch_and_price(extra_forbidden={(15, 16)}, root_K_ub=3, root_K_lb=3)
