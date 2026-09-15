# -*- coding: utf-8 -*-
"""csp50 单 crew 调度子问题原型（CP-SAT）：给定任务集 S，求可行序列（弧兼容 + 跨度<=480）与最小成本。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import math, time
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt"
lines = open(DATA).read().splitlines()
N, T = map(int, lines[0].split())
tasks = [None] + [tuple(map(int, l.split())) for l in lines[1:1+N]]
arc_cost = {}
for l in lines[1+N:]:
    i, j, c = map(int, l.split()); arc_cost[(i, j)] = c
M = 10**9

_cache = {}
def solve_sched(S, time_limit=5.0):
    """返回 (status, route, cost)。status ∈ OPTIMAL/FEASIBLE/INFEASIBLE/UNKNOWN。"""
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
    # 允许弧：0->i（开始）、i->j（arc_list）、i->N+1（结束）、N+1->0（闭合）
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
    for a in range(m_):          # 未给出的弧禁用
        for b in range(m_):
            if a == b:
                continue
            i, j = nodes[a], nodes[b]
            ok = (i == 0 and 1 <= j <= N) or (1 <= i <= N and 1 <= j <= N and (i, j) in arc_cost) \
                 or (1 <= i <= N and j == N+1) or (i == N+1 and j == 0)
            if not ok:
                model.Add(xv[a][b] == 0)
    # 跨度：ts=首个任务开始, te=末个任务结束
    ts = model.NewIntVar(0, 10**9, "ts")
    te = model.NewIntVar(0, 10**9, "te")
    for a in range(1, m_-1):
        i = nodes[a]
        model.Add(ts >= tasks[i][0] - M*(1 - xv[0][a]))
        model.Add(ts <= tasks[i][0] + M*(1 - xv[0][a]))
        model.Add(te >= tasks[i][1] - M*(1 - xv[a][m_-1]))
        model.Add(te <= tasks[i][1] + M*(1 - xv[a][m_-1]))
    model.Add(te - ts <= T)
    model.Minimize(sum(arc_cost.get((nodes[a], nodes[b]), 0) * xv[a][b]
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
    cost = solver.ObjectiveValue()
    _cache[S] = (status, tuple(route), cost)
    return _cache[S]

# 验证：最优解中的若干 crew 集合（来自 pool_opt_result.json）
tests = [(1, 10), (2, 13), (3,), (5, 19), (6, 15), (7, 20), (1, 2, 3), (1, 10, 2)]
for S in tests:
    st, route, cost = solve_sched(S)
    print(f"S={S}: {st} route={route} cost={cost}")
# 不可行集示例：把 50 个任务全塞一个 crew
st, route, cost = solve_sched(tuple(range(1, 51)))
print("全部 50 任务一 crew:", st)
# 弧兼容性检查：任意 (i,j) 弧是否都满足 end_i <= start_j
bad = [(i, j) for (i, j) in arc_cost if tasks[i][1] > tasks[j][0]]
print("时间不兼容弧数:", len(bad), bad[:5])
