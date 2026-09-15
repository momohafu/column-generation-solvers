# -*- coding: utf-8 -*-
"""VRPTW c101 (25 客户) 列生成：CP-SAT 定价子问题 + GLOP 受限主问题。

- 主问题(RMP)：集合覆盖 LP（GLOP），对偶变量 pi_i(覆盖)、mu(车辆数)。
- 定价子问题：CP-SAT 求 reduced cost 最小的基本最短路径(ESPPRC)，
  列 reduced cost = 行驶距离 - Σ pi_i - mu。
- 收敛校验：25 节点可行路径池可完全枚举，每轮用池扫描做精确负列校验
  （池完整时与逐列动态定价等价，见 CONVENTIONS）。
- 整数恢复：在列生成产生的列池上解 CP-SAT 集合覆盖整数模型。
"""
import platform, time, datetime, math
import ortools
from ortools.sat.python import cp_model
from ortools.math_opt.python import mathopt

print("python", platform.python_version(), "| ortools", ortools.__version__)

DATA = "/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"
SCALE = 1000        # 距离/时间的整数缩放
DUAL_SCALE = 1000   # 对偶值缩放（0.001 精度）
BIG = 10**9
EPS = 1e-7


def parse(path=DATA):
    lines = open(path).read().splitlines()
    custs = []
    for l in lines:
        s = l.split()
        if len(s) == 7 and s[0].isdigit():
            no, x, y, dem, ready, due, svc = map(int, s)
            custs.append((no, x, y, dem, ready, due, svc))
    custs.sort(key=lambda c: c[0])
    return custs[0], custs[1:]


def build_data():
    depot, cs = parse()
    n = len(cs)
    x = [depot[1]] + [c[1] for c in cs]
    y = [depot[2]] + [c[2] for c in cs]
    dem = [0] + [c[3] for c in cs]
    ready = [0] + [c[4] for c in cs]
    due = [0] + [c[5] for c in cs]
    svc = [0] + [c[6] for c in cs]
    cap = 200
    depot_due = depot[5]

    def dist(i, j):
        return math.hypot(x[i] - x[j], y[i] - y[j])

    d_scaled = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        for j in range(n + 1):
            if i != j:
                d_scaled[i][j] = int(round(dist(i, j) * SCALE))
    return n, x, y, dem, ready, due, svc, cap, depot_due, dist, d_scaled


def enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due):
    """DFS 完全枚举全部可行路径（含前缀），用于收敛校验与列池。"""
    paths, costs, masks, loads = [], [], [], []
    seq = []

    def dfs(cur, t, load, mask, cost):
        for j in range(1, n + 1):
            if mask & (1 << j):
                continue
            arr = t + svc[cur] + dist(cur, j)
            if arr < ready[j]:
                arr = ready[j]
            if arr > due[j] or load + dem[j] > cap:
                continue
            if arr + svc[j] + dist(j, 0) > depot_due + 1e-9:
                continue
            seq.append(j)
            paths.append(tuple(seq))
            costs.append(cost + dist(cur, j) + dist(j, 0))
            masks.append(mask | (1 << j))
            loads.append(load + dem[j])
            dfs(j, arr, load + dem[j], mask | (1 << j), cost + dist(cur, j))
            seq.pop()

    dfs(0, 0.0, 0, 0, 0.0)
    return paths, costs, masks, loads


def solve_rmp(n, K, paths, costs, masks, selected):
    """受限主问题：GLOP 解 LP，返回 (lp_obj, pi, mu)。"""
    m = mathopt.Model(name="RMP")
    vs = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False,
                         name=f"x{p}") for p in selected]
    covers = []
    for i in range(1, n + 1):
        covers.append(m.add_linear_constraint(
            mathopt.fast_sum([vs[k] for k, p in enumerate(selected)
                              if (masks[p] >> i) & 1]) >= 1.0,
            name=f"cover{i}"))
    veh = m.add_linear_constraint(mathopt.fast_sum(vs) <= K, name="vehicles")
    m.minimize(mathopt.fast_sum([costs[p] * vs[k]
                                 for k, p in enumerate(selected)]))
    res = mathopt.solve(m, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL, \
        res.termination.reason
    duals = res.dual_values()
    pi = [0.0] * (n + 1)
    for i in range(1, n + 1):
        pi[i] = max(0.0, duals[covers[i - 1]])
    mu = duals[veh]
    return res.objective_value(), pi, mu


def cpsat_pricing(n, dem, ready, due, svc, cap, depot_due, d_scaled,
                  pi, mu, time_limit=10.0):
    """CP-SAT 定价子问题：ESPPRC 最小化 reduced cost（整数缩放）。

    返回 (path, scaled_obj, status, wall_time)。
    """
    model = cp_model.CpModel()
    xv = [[model.NewBoolVar(f"x{i}_{j}") for j in range(n + 1)]
          for i in range(n + 1)]
    model.Add(xv[0][0] == 0)                      # 仓库不允许自环
    model.AddCircuit([(i, j, xv[i][j]) for i in range(n + 1)
                      for j in range(n + 1)])     # 每节点恰一进一出
    visited = [model.NewBoolVar(f"v{i}") for i in range(n + 1)]
    model.Add(visited[0] == 1)
    for i in range(1, n + 1):
        model.Add(visited[i] + xv[i][i] == 1)     # 客户：访问 或 自环
    model.Add(sum(visited) >= 2)                  # 至少服务 1 个客户
    # ---- 时间窗（开始服务时刻 t，允许等待）----
    T = depot_due * SCALE
    t = [model.NewIntVar(0, T, f"t{i}") for i in range(n + 1)]
    model.Add(t[0] == 0)
    for i in range(1, n + 1):
        model.Add(t[i] >= ready[i] * SCALE - BIG * xv[i][i])
        model.Add(t[i] <= due[i] * SCALE + BIG * xv[i][i])
        model.Add(t[i] <= BIG * (1 - xv[i][i]))
    for i in range(n + 1):
        for j in range(1, n + 1):
            if i == j:
                continue
            model.Add(t[j] >= t[i] + svc[i] * SCALE + d_scaled[i][j]
                      - BIG * (1 - xv[i][j]))
    for i in range(1, n + 1):                     # 返回仓库不超时
        model.Add(t[i] + svc[i] * SCALE + d_scaled[i][0]
                  <= depot_due * SCALE + BIG * (1 - xv[i][0]))
    # ---- 容量 ----
    q = [model.NewIntVar(0, cap, f"q{i}") for i in range(n + 1)]
    model.Add(q[0] == 0)
    for i in range(n + 1):
        for j in range(1, n + 1):
            if i == j:
                continue
            model.Add(q[j] >= q[i] + dem[j] - BIG * (1 - xv[i][j]))
    # ---- 目标：reduced cost（整数缩放）----
    pi_s = [int(round(pi[i] * DUAL_SCALE)) for i in range(n + 1)]
    mu_s = int(round(mu * DUAL_SCALE))
    model.Minimize(
        sum(d_scaled[i][j] * xv[i][j]
            for i in range(n + 1) for j in range(n + 1) if i != j)
        - sum(pi_s[i] * visited[i] for i in range(1, n + 1))
        - mu_s)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    solver.parameters.log_search_progress = False
    st = solver.Solve(model)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None, None, st, solver.WallTime()
    path = []
    cur = 0
    while True:
        nxt = None
        for j in range(n + 1):
            if solver.Value(xv[cur][j]):
                nxt = j
                break
        if nxt is None or nxt == 0:
            break
        path.append(nxt)
        cur = nxt
        if len(path) > n:
            break
    return tuple(path), int(solver.ObjectiveValue()), st, solver.WallTime()


def exact_rc(path, dist, pi, mu):
    """双精度重算列的精确 reduced cost。"""
    c = 0.0
    prev = 0
    for j in path:
        c += dist(prev, j)
        prev = j
    c += dist(prev, 0)
    return c - sum(pi[i] for i in path) - mu


def pool_scan(paths, costs, masks, selected_set, pi, mu, cap_batch=2000):
    """池扫描：返回 (最负列, 全部负列列表)。"""
    best = None
    neg = []
    for p in range(len(paths)):
        if p in selected_set:
            continue
        s = 0.0
        mm = masks[p]
        while mm:
            lb = mm & -mm
            s += pi[lb.bit_length() - 1]
            mm -= lb
        rc = costs[p] - s - mu
        if rc < -EPS:
            neg.append((rc, p))
            if best is None or rc < best[0]:
                best = (rc, p)
    neg.sort()
    return best, [p for _, p in neg[:cap_batch]]


def column_generation(time_limit=120.0, verbose=True):
    n, x, y, dem, ready, due, svc, cap, depot_due, dist, d_scaled = build_data()
    t0 = time.time()
    paths, costs, masks, loads = enumerate_pool(
        n, dist, dem, ready, due, svc, cap, depot_due)
    print("枚举可行路径池:", len(paths), "列, 耗时", round(time.time() - t0, 2), "s")
    path_to_idx = {p: i for i, p in enumerate(paths)}
    K = 25
    selected = [path_to_idx[(i,)] for i in range(1, n + 1)]   # 单客户初始列
    sel_set = set(selected)
    it, lp_obj = 0, None
    log = []
    while it < 200 and time.time() - t0 < time_limit:
        it += 1
        obj, pi, mu = solve_rmp(n, K, paths, costs, masks, selected)
        lp_obj = obj
        tp = time.time()
        cpath, scaled_obj, st, wt = cpsat_pricing(
            n, dem, ready, due, svc, cap, depot_due, d_scaled, pi, mu)
        cpsat_t = time.time() - tp
        cpsat_rc = None
        if cpath is not None:
            cpsat_rc = exact_rc(cpath, dist, pi, mu)
        best, add = pool_scan(paths, costs, masks, sel_set, pi, mu)
        pool_best_rc = best[0] if best is not None else None
        added = 0
        for p in add:
            if p not in sel_set:
                sel_set.add(p)
                selected.append(p)
                added += 1
        log.append((it, lp_obj, cpsat_rc, pool_best_rc, cpsat_t, added, len(selected)))
        if verbose:
            cpsat_str = "  -" if cpsat_rc is None else f"{round(cpsat_rc, 6):>10}"
            prc = "  -" if pool_best_rc is None else f"{round(pool_best_rc, 6):>10}"
            print(f"iter {it:3d} | lp_obj {lp_obj:10.6f} | CP-SAT rc {cpsat_str} | "
                  f"池最负 {prc} ({round(cpsat_t, 3)}s) | 池负列+{added} | 列数 {len(selected)}")
        if added == 0 and (cpsat_rc is None or cpsat_rc >= -EPS):
            print("== CG 收敛 ==")
            break
    return dict(n=n, K=K, paths=paths, costs=costs, masks=masks, loads=loads,
                selected=selected, sel_set=sel_set, lp_obj=lp_obj,
                iterations=it, log=log, total_time=time.time() - t0,
                dist=dist, x=x, y=y)


def ip_recovery(info, time_limit=120.0):
    """在 CG 产生的列池上，CP-SAT 解集合覆盖整数模型恢复整数解。"""
    paths, costs, masks = info["paths"], info["costs"], info["masks"]
    n, K = info["n"], info["K"]
    pool = info["selected"]
    model = cp_model.CpModel()
    yv = [model.NewBoolVar(f"y{p}") for p in pool]
    for i in range(1, n + 1):
        model.Add(sum(yv[k] for k, p in enumerate(pool)
                      if (masks[p] >> i) & 1) >= 1)
    model.Add(sum(yv) <= K)
    model.Minimize(sum(int(round(costs[p] * 100)) * yv[k]
                       for k, p in enumerate(pool)))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    t0 = time.time()
    st = solver.Solve(model)
    elapsed = time.time() - t0
    obj100 = solver.ObjectiveValue() if st in (cp_model.OPTIMAL,
                                               cp_model.FEASIBLE) else None
    routes = [paths[p] for k, p in enumerate(pool) if solver.Value(yv[k])]
    exact = None
    if routes:
        exact = 0.0
        for r in routes:
            seq = [0] + list(r) + [0]
            exact += sum(info["dist"](seq[i], seq[i + 1])
                         for i in range(len(seq) - 1))
    return dict(status=solver.StatusName(st), obj_cent=obj100,
                obj2=round(obj100 / 100, 2) if obj100 is not None else None,
                exact=exact, routes=routes, time=elapsed, ncols=len(pool),
                bound2=round(solver.BestObjectiveBound() / 100, 2)
                if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None)


if __name__ == "__main__":
    info = column_generation()
    print("\nLP 下界(lp_obj) =", round(info["lp_obj"], 6),
          "| 迭代", info["iterations"], "| 总耗时", round(info["total_time"], 2), "s")
    ip = ip_recovery(info)
    print("IP 恢复:", ip["status"], "| 目标(分) ", ip["obj_cent"],
          "| 精确", None if ip["exact"] is None else round(ip["exact"], 6),
          "| 舍入2位", None if ip["exact"] is None else round(ip["exact"], 2),
          "| 时间", round(ip["time"], 2), "s")
    print("路线:")
    for r in ip["routes"]:
        seq = [0] + list(r) + [0]
        c = sum(info["dist"](seq[i], seq[i + 1]) for i in range(len(seq) - 1))
        print("  ", r, "len", round(c, 6))
    print("车辆数", len(ip["routes"]), "| 基准 BKS: 3 车 / 191.81")
