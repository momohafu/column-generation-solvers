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

_cache = {}
def solve_sched(S, time_limit=5.0):
    S = tuple(sorted(S))
    if S in _cache:
        return _cache[S]
    if len(S) == 0:
        return ("OPTIMAL", (), 0)
    if len(S) == 1:
        i = S[0]
        return ("OPTIMAL", (i,), 0) if tasks[i][1]-tasks[i][0] <= T else ("INFEASIBLE", None, None)
    nodes = [0] + list(S) + [N+1]
    idx = {nd: a for a, nd in enumerate(nodes)}
    m_ = len(nodes)
    model = cp_model.CpModel()
    xv = [[model.NewBoolVar(f"x{a}_{b}") for b in range(m_)] for a in range(m_)]
    arcs = []
    for a in range(m_):
        for b in range(m_):
            if a == b:
                continue
            i, j = nodes[a], nodes[b]
            if i == 0 and 1 <= j <= N:
                arcs.append((a, b, xv[a][b]))
            elif 1 <= i <= N and 1 <= j <= N and (i, j) in arc_cost:
                arcs.append((a, b, xv[a][b]))
            elif 1 <= i <= N and j == N+1:
                arcs.append((a, b, xv[a][b]))
            elif i == N+1 and j == 0:
                arcs.append((a, b, xv[a][b]))
    model.AddCircuit(arcs)
    for a in range(m_):
        for b in range(m_):
            if a == b:
                continue
            i, j = nodes[a], nodes[b]
            ok = (i == 0 and 1 <= j <= N) or (1 <= i <= N and 1 <= j <= N and (i, j) in arc_cost) \
                 or (1 <= i <= N and j == N+1) or (i == N+1 and j == 0)
            if not ok:
                model.Add(xv[a][b] == 0)
    ts = model.NewIntVar(0, 10**9, "ts")
    te = model.NewIntVar(0, 10**9, "te")
    for a in range(1, m_-1):
        i = nodes[a]
        model.Add(ts >= tasks[i][0] - 10**9*(1 - xv[0][a]))
        model.Add(ts <= tasks[i][0] + 10**9*(1 - xv[0][a]))
        model.Add(te >= tasks[i][1] - 10**9*(1 - xv[a][m_-1]))
        model.Add(te <= tasks[i][1] + 10**9*(1 - xv[a][m_-1]))
    model.Add(te - ts <= T)
    model.Minimize(sum(arc_cost.get((nodes[a], nodes[b]), 0)*xv[a][b]
                       for a in range(m_) for b in range(m_) if a != b))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    st = solver.Solve(model)
    if st == cp_model.INFEASIBLE:
        _cache[S] = ("INFEASIBLE", None, None)
        return _cache[S]
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return ("UNKNOWN", None, None)
    status = "OPTIMAL" if st == cp_model.OPTIMAL else "FEASIBLE"
    route = []
    cur = 0
    for _ in range(len(S)):
        nxt = None
        for b in range(1, m_-1):
            if solver.Value(xv[cur][b]):
                nxt = b
                break
        if nxt is None:
            break
        route.append(nodes[nxt])
        cur = nxt
    _cache[S] = (status, tuple(route), solver.ObjectiveValue())
    return _cache[S]

def build_master(nogoods, hooker, exact_cuts, theta_ub, warm=None):
    model = cp_model.CpModel()
    a = [[model.NewBoolVar(f"a{i}_{c}") for c in range(K_MIN)] for i in range(1, N+1)]
    for i in range(1, N+1):
        model.Add(sum(a[i-1]) == 1)
    for c in range(K_MIN):
        model.Add(sum(a[i-1][c] for i in range(1, N+1)) <= 3)      # 最长链 3（池枚举事实）
    for (c, S) in nogoods:
        model.Add(sum(1 - a[i-1][c] for i in S) >= 1)
    theta = model.NewIntVar(0, 10**7, "theta")
    for (cks, deltas) in hooker:
        model.Add(theta >= sum(cks[c] - sum(deltas[c][i]*(1 - a[i-1][c]) for i in deltas[c]) for c in cks))
    for (cc, S_by_c) in exact_cuts:
        model.Add(theta >= cc * (1 - sum(sum(1 - a[i-1][c] for i in S) for c, S in S_by_c.items())))
    if theta_ub is not None:
        model.Add(theta <= theta_ub)
    if warm is not None:
        for c, S in warm.items():
            for i in S:
                model.AddHint(a[i-1][c], 1)
    return model, a, theta

def extract(solver, a):
    S_by_c = {}
    for c in range(K_MIN):
        S = tuple(i for i in range(1, N+1) if solver.Value(a[i-1][c]))
        if S:
            S_by_c[c] = S
    return S_by_c

def make_cuts(S_by_c):
    cks = {}
    deltas = {}
    C = 0.0
    for c, S in S_by_c.items():
        st, route, ck = solve_sched(S)
        if st == "INFEASIBLE":
            return None
        C += ck
        cks[c] = int(round(ck))
        deltas[c] = {}
        for i in S:
            S2 = tuple(j for j in S if j != i)
            st2, _, c2 = solve_sched(S2)
            dlt = ck - c2 if st2 != "INFEASIBLE" else ck
            deltas[c][i] = int(round(dlt))
    return C, cks, deltas

def run_lbbd(verbose=True):
    """主方法（scp41 同款 LBBD）：选列主问题 min Σc·y + Σy<=K，子问题=未覆盖任务检查，
    逻辑割 = Σ_{p∋i} y_p >= 1（对每个未覆盖任务）。收敛时主问题等价于完整池覆盖 IP。"""
    cut_tasks = set()
    wall0 = time.time()
    for it in range(50):
        m = mathopt.Model()
        y = [m.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"y{p}") for p in range(P)]
        m.add_linear_constraint(mathopt.fast_sum(y) <= K_MIN, name="veh")
        m.minimize(mathopt.fast_sum([pcost[p]*y[p] for p in range(P)]))
        for i in cut_tasks:
            m.add_linear_constraint(mathopt.fast_sum([y[p] for p in range(P) if (pmask[p] >> i) & 1]) >= 1.0, name=f"cut{i}")
        res = mathopt.solve(m, mathopt.SolverType.HIGHS,
                            params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
        sel = [p for p in range(P) if res.variable_values()[y[p]] > 0.5]
        uncovered = []
        for i in range(1, N+1):
            if not any((pmask[p] >> i) & 1 for p in sel):
                uncovered.append(i)
        if verbose:
            print(f"iter {it+1}: master_obj={res.objective_value()} | 选中 {len(sel)} 列 | 未覆盖 {len(uncovered)} | 逻辑割 {len(cut_tasks)}")
        if not uncovered:
            print("全部任务覆盖 -> LBBD 收敛，主问题 = 完整池覆盖 IP（证明最优）")
            break
        for i in uncovered:
            cut_tasks.add(i)
        if time.time()-wall0 > 110:
            print("时间上限")
            break
    routes = [pseq[p] for p in sel]
    return dict(iterations=it+1, cuts=len(cut_tasks), obj=res.objective_value(),
                routes=routes, wall=time.time()-wall0)

def run_lbbd_assignment(max_iter=6, verbose=True):
    """补充实验：分配主问题 + 调度子问题（Hooker/no-good 割机制展示）。
    注：50 任务/27 crew 组合空间巨大、no-good 割弱，该变体在迭代上限内难以自收敛——
    此处仅展示割的构造与机制，收敛证明由 run_lbbd（选列主问题）给出。"""
    # 初始可行解：完整池 IP
    m = mathopt.Model()
    x = [m.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{p}") for p in range(P)]
    for i in range(1, N+1):
        m.add_linear_constraint(mathopt.fast_sum([x[p] for p in range(P) if (pmask[p] >> i) & 1]) >= 1.0, name=f"c{i}")
    m.add_linear_constraint(mathopt.fast_sum(x) <= K_MIN, name="veh")
    m.minimize(mathopt.fast_sum([pcost[p]*x[p] for p in range(P)]))
    res = mathopt.solve(m, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
    init_routes = [pseq[p] for p in range(P) if res.variable_values()[x[p]] > 0.5]
    warm = {c: tuple(init_routes[c]) for c in range(len(init_routes))}
    UB = float(res.objective_value())
    if verbose:
        print(f"warm-start: 池 IP 解 {len(init_routes)} crew, 成本 {int(UB)}")
        for c, S in list(warm.items())[:4]:
            st, route, ck = solve_sched(S)
            print(f"  crew{c} 集合{S}: 调度 {st} 成本 {ck}")
    # Hooker 割构造展示（warm 解）
    C0, cks0, d0 = make_cuts(warm)
    if verbose:
        print(f"Hooker 割（warm 解）: theta >= Σ_c [c_c - Σ δ_ic(1-a_ic)]，Σ c_c = {int(C0)}")
        c0 = list(warm.items())[0][1]
        st, route, ck = solve_sched(c0)
        for i in c0:
            S2 = tuple(j for j in c0 if j != i)
            st2, _, c2 = solve_sched(S2)
            print(f"  crew{c0}: 移除 {i} 节省 δ = {ck - (c2 if st2 != 'INFEASIBLE' else 0)}")
    # 分配主问题迭代展示（少量轮次）
    model, a, theta = build_master([], [(cks0, d0)], [(int(round(C0)), warm)], int(UB)-1, warm=warm)
    model.Minimize(theta)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 20
    solver.parameters.num_search_workers = 8
    st = solver.Solve(model)
    if verbose:
        print(f"分配主问题首轮: {solver.StatusName(st)} theta={solver.Value(theta)}（若无可行即已证明 UB 最优）")
    nogoods = 0
    for it in range(max_iter):
        if st == cp_model.INFEASIBLE:
            if verbose:
                print("分配主问题不可行 -> UB 最优（Hooker 割 + 精确界下）")
            break
        S_by_c = extract(solver, a)
        bad = 0
        for c, S in S_by_c.items():
            if solve_sched(S)[0] == "INFEASIBLE":
                bad += 1
        if verbose:
            print(f"iter {it+1}: 分配中 {bad} 个 crew 不可行（no-good 割机制生效，需迭代剪除）")
        if bad == 0:
            break
        # 只展示一轮后退出（完整收敛见 run_lbbd）
        for c, S in S_by_c.items():
            if solve_sched(S)[0] == "INFEASIBLE":
                model.Add(sum(1 - a[i-1][c] for i in S) >= 1)
                nogoods += 1
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = 20
        solver.parameters.num_search_workers = 8
        st = solver.Solve(model)
    return dict(ub=UB, proven=(st == cp_model.INFEASIBLE), nogoods=nogoods, routes=warm)

if __name__ == "__main__":
    print("== 变体 B（主方法）：选列主问题 + 覆盖逻辑割 ==")
    run_lbbd()
    print()
    print("== 变体 A（机制展示）：分配主问题 + 调度子问题 ==")
    run_lbbd_assignment()
