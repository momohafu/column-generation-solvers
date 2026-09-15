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

def solve_rmp(selected):
    m = mathopt.Model()
    x = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in selected]
    yv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(1, N+1)]
    covers = []
    for i in range(1, N+1):
        covers.append(m.add_linear_constraint(
            mathopt.fast_sum([x[k] for k, p in enumerate(selected) if (pmask[p] >> i) & 1]) + yv[i-1] >= 1.0, name=f"c{i}"))
    veh = m.add_linear_constraint(mathopt.fast_sum(x) <= K_MIN, name="veh")
    m.minimize(mathopt.fast_sum([pcost[p]*x[k] for k, p in enumerate(selected)]) + M_DUMMY*mathopt.fast_sum(yv))
    res = mathopt.solve(m, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
    dv = res.dual_values()
    pi = [0.0]*(N+1)
    for i in range(1, N+1):
        pi[i] = max(0.0, dv[covers[i-1]])
    mu = dv[veh]
    xvals = {p: res.variable_values()[x[k]] for k, p in enumerate(selected)}
    y_used = sum(res.variable_values()[yv[i-1]] for i in range(1, N+1))
    return res.objective_value(), pi, mu, xvals, y_used

def cpsat_pricing(pi, mu_sum, time_limit=5.0):
    """单 crew 定价：可选任务路径 min(成本 - Σπ - μ)，跨度<=T。"""
    model = cp_model.CpModel()
    nn = N + 2
    xv = [[model.NewBoolVar(f"x{a}_{b}") for b in range(nn)] for a in range(nn)]
    arcs = []
    for a in range(nn):
        for b in range(nn):
            if a == b:
                continue
            i, j = (a if a == 0 else a), (b if b == nn-1 else b)
            if a == 0 and 1 <= b <= N:
                arcs.append((a, b, xv[a][b]))
            elif 1 <= a <= N and 1 <= b <= N and (a, b) in arc_cost:
                arcs.append((a, b, xv[a][b]))
            elif 1 <= a <= N and b == nn-1:
                arcs.append((a, b, xv[a][b]))
            elif a == nn-1 and b == 0:
                arcs.append((a, b, xv[a][b]))
    for a in range(1, N+1):
        arcs.append((a, a, xv[a][a]))          # 未访问任务自环
    model.AddCircuit(arcs)
    for a in range(nn):
        for b in range(nn):
            if a == b:
                continue
            ok = (a == 0 and 1 <= b <= N) or (1 <= a <= N and 1 <= b <= N and (a, b) in arc_cost) \
                 or (1 <= a <= N and b == nn-1) or (a == nn-1 and b == 0)
            if not ok:
                model.Add(xv[a][b] == 0)
    visited = [model.NewBoolVar(f"v{a}") for a in range(1, N+1)]
    for a in range(1, N+1):
        model.Add(visited[a-1] + xv[a][a] == 1)
    model.Add(sum(visited) >= 1)
    model.Add(xv[0][0] == 0)
    model.Add(xv[nn-1][nn-1] == 0)
    ts = model.NewIntVar(0, 10**9, "ts")
    te = model.NewIntVar(0, 10**9, "te")
    for a in range(1, N+1):
        model.Add(ts >= tasks[a][0] - 10**9*(1 - xv[0][a]))
        model.Add(ts <= tasks[a][0] + 10**9*(1 - xv[0][a]))
        model.Add(te >= tasks[a][1] - 10**9*(1 - xv[a][nn-1]))
        model.Add(te <= tasks[a][1] + 10**9*(1 - xv[a][nn-1]))
    model.Add(te - ts <= T)
    pi_s = [int(round(pi[a]*1000)) for a in range(1, N+1)]
    mu_s = int(round(mu_sum*1000))
    model.Minimize(sum(arc_cost.get((a, b), 0)*xv[a][b] for a in range(nn) for b in range(nn) if a != b)
                   - sum(pi_s[a-1]*visited[a-1] for a in range(1, N+1)) - mu_s)
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
        for b in range(1, N+1):
            if solver.Value(xv[cur][b]):
                nxt = b
                break
        if nxt is None:
            break
        path.append(nxt)
        cur = nxt
        if len(path) > N:
            break
    return tuple(path)

def exact_rc(path, pi, mu_sum):
    c = sum(arc_cost.get((path[k], path[k+1]), 0) for k in range(len(path)-1))
    return c - sum(pi[i] for i in path) - mu_sum

def run_cg(verbose=True):
    sel = [p for p in range(P) if len(pseq[p]) == 1]
    sel_set = set(sel)
    lp_obj = None
    for it in range(200):
        obj, pi, mu, xvals, y_used = solve_rmp(sel)
        lp_obj = obj
        cpath = cpsat_pricing(pi, mu)
        cpsat_rc = exact_rc(cpath, pi, mu) if cpath is not None else None
        add = []
        for p in range(P):
            if p in sel_set:
                continue
            s2 = 0.0
            mm = pmask[p]
            while mm:
                lb = mm & -mm
                s2 += pi[lb.bit_length()-1]
                mm -= lb
            rc = pcost[p] - s2 - mu
            if rc < -EPS:
                add.append((rc, p))
        add.sort()
        add = [p for _, p in add[:1000]]
        if verbose:
            print(f"iter {it+1}: LP={round(lp_obj,4)} | CP-SAT rc={'-' if cpsat_rc is None else round(cpsat_rc,4)} | 加列 {len(add)} | 列 {len(sel)}")
        if not add and (cpsat_rc is None or cpsat_rc >= -EPS):
            break
        for p in add:
            if p not in sel_set:
                sel_set.add(p)
                sel.append(p)
    # 整数恢复：列池 CP-SAT
    m = cp_model.CpModel()
    yv = [m.NewBoolVar(f"y{p}") for p in sel]
    for i in range(1, N+1):
        m.Add(sum(yv[k] for k, p in enumerate(sel) if (pmask[p] >> i) & 1) >= 1)
    m.Add(sum(yv) <= K_MIN)
    m.Minimize(sum(pcost[p]*yv[k] for k, p in enumerate(sel)))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = 60
    st = s.Solve(m)
    routes = [pseq[p] for k, p in enumerate(sel) if s.Value(yv[k])]
    ip_obj = s.ObjectiveValue() if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    if verbose:
        print(f"LP 下界 = {round(lp_obj,4)} | 整数恢复 {s.StatusName(st)} obj={ip_obj} | crew {len(routes)}")
    return dict(lp=lp_obj, ip=ip_obj, iterations=it+1, routes=routes, status=s.StatusName(st))

if __name__ == "__main__":
    run_cg()
