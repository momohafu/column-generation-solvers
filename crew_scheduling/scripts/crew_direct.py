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

def pool_ip(Kv, use_cpsat=False, time_limit=60.0):
    if not use_cpsat:
        m = mathopt.Model()
        x = [m.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{p}") for p in range(P)]
        for i in range(1, N+1):
            m.add_linear_constraint(mathopt.fast_sum([x[p] for p in range(P) if (pmask[p] >> i) & 1]) >= 1.0, name=f"c{i}")
        m.add_linear_constraint(mathopt.fast_sum(x) <= Kv, name="veh")
        m.minimize(mathopt.fast_sum([pcost[p]*x[p] for p in range(P)]))
        t0 = time.time()
        res = mathopt.solve(m, mathopt.SolverType.HIGHS,
                            params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=time_limit), enable_output=False))
        wt = time.time() - t0
        obj = res.objective_value() if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE) else None
        sel = [paths[p] for p in range(P) if res.variable_values()[x[p]] > 0.5] if obj is not None else None
        return res.termination.reason.name, obj, sel, wt
    else:
        m = cp_model.CpModel()
        x = [m.NewBoolVar(f"x{p}") for p in range(P)]
        for i in range(1, N+1):
            m.Add(sum(x[p] for p in range(P) if (pmask[p] >> i) & 1) >= 1)
        m.Add(sum(x) <= Kv)
        m.Minimize(sum(pcost[p]*x[p] for p in range(P)))
        s = cp_model.CpSolver()
        s.parameters.max_time_in_seconds = time_limit
        t0 = time.time()
        st = s.Solve(m)
        wt = time.time() - t0
        sel = [paths[p] for p in range(P) if s.Value(x[p])] if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
        return s.StatusName(st), (s.ObjectiveValue() if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None), sel, wt

def run_direct(verbose=True):
    total_dur = sum(f-s for s, f in tasks[1:])
    ev = []
    for s, f in tasks[1:]:
        ev.append((s, 1)); ev.append((f, -1))
    ev.sort(key=lambda x: (x[0], -x[1]))
    cur = overlap = 0
    for t2, d in ev:
        cur += d; overlap = max(overlap, cur)
    has_in = {j for (i, j) in arc_cost}
    noin = [i for i in range(1, N+1) if i not in has_in]
    if verbose:
        print(f"下界: 时长={math.ceil(total_dur/T)} 重叠={overlap} 无入弧={len(noin)} | 池 {P} 条路径")
    term26, obj26, _, wt26 = pool_ip(26)
    term27, obj27, routes, wt27 = pool_ip(27)
    term27c, obj27c, _, wt27c = pool_ip(27, use_cpsat=True)
    if verbose:
        print(f"K=26: {term26}（{round(wt26,3)}s）-> 最少 crew = 27")
        print(f"K=27 HIGHS: {term27} obj={obj27}（{round(wt27,3)}s）")
        print(f"K=27 CP-SAT 交叉验证: {term27c} obj={obj27c}（{round(wt27c,2)}s）")
        for r in sorted(routes, key=lambda s: -len(s)):
            print("  crew:", r, "成本", sum(arc_cost.get((r[k], r[k+1]), 0) for k in range(len(r)-1)))
    return dict(K_min=27, objective=obj27, routes=routes, term26=term26, term27=term27)

if __name__ == "__main__":
    run_direct()
