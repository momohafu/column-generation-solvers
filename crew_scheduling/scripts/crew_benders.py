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

def solve_sp(yvec):
    sp = mathopt.Model()
    xv = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
    sv = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"s{i}") for i in range(1, N+1)]
    covs = []
    for i in range(1, N+1):
        covs.append(sp.add_linear_constraint(
            mathopt.fast_sum([xv[p] for p in range(P) if (pmask[p] >> i) & 1]) + sv[i-1] >= 1.0, name=f"c{i}"))
    veh = sp.add_linear_constraint(mathopt.fast_sum(xv) <= K_MIN, name="veh")
    xub = [sp.add_linear_constraint(xv[p] <= float(yvec[p]), name=f"xub{p}") for p in range(P)]
    sp.minimize(mathopt.fast_sum([pcost[p]*xv[p] for p in range(P)]) + M_DUMMY*mathopt.fast_sum(sv))
    res = mathopt.solve(sp, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
    dv = res.dual_values()
    pi = [max(0.0, dv[covs[i-1]]) for i in range(1, N+1)]
    mu = dv[veh]
    alpha = sum(pi) + K_MIN*mu
    lam = [-dv[xub[p]] for p in range(P)]
    return alpha, lam, res.objective_value()

def solve_master(cuts):
    mp = mathopt.Model()
    yv = [mp.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"y{p}") for p in range(P)]
    theta = mp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name="theta")
    for kk, (alpha, lam) in enumerate(cuts):
        mp.add_linear_constraint(theta + mathopt.fast_sum([lam[p]*yv[p] for p in range(P) if lam[p] != 0.0]) >= alpha, name=f"bcut{kk}")
    mp.minimize(theta)
    res = mathopt.solve(mp, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=20), enable_output=False))
    vals = res.variable_values(yv)
    ystar = [1 if vals[p] > 0.5 else 0 for p in range(P)]
    th = res.variable_values([theta])[0]
    return ystar, th, res.objective_value()

def run_benders(verbose=True):
    current = [1]*P
    cuts = []
    wall0 = time.time()
    for it in range(20):
        alpha, lam, sp_obj = solve_sp(current)
        cuts.append((alpha, lam))
        ystar, theta, mp_obj = solve_master(cuts)
        if verbose:
            print(f"iter {it+1}: SP={round(sp_obj,4)} | MP={round(mp_obj,4)} theta={round(theta,4)} | y选中 {sum(ystar)} | 割 {len(cuts)}")
        if ystar == current:
            break
        current = ystar
        if time.time()-wall0 > 110:
            break
    # 整数修复（完整池 HIGHS）
    m = mathopt.Model()
    x = [m.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{p}") for p in range(P)]
    for i in range(1, N+1):
        m.add_linear_constraint(mathopt.fast_sum([x[p] for p in range(P) if (pmask[p] >> i) & 1]) >= 1.0, name=f"c{i}")
    m.add_linear_constraint(mathopt.fast_sum(x) <= K_MIN, name="veh")
    m.minimize(mathopt.fast_sum([pcost[p]*x[p] for p in range(P)]))
    res = mathopt.solve(m, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
    routes = [pseq[p] for p in range(P) if res.variable_values()[x[p]] > 0.5]
    if verbose:
        print(f"整数修复: {res.termination.reason.name} obj={res.objective_value()} crew={len(routes)}")
        print(f"Benders 下界 {round(mp_obj,4)} vs 修复 {res.objective_value()} | gap {(res.objective_value()-mp_obj)/mp_obj*100 if mp_obj else 0:.4f}%")
    return dict(iterations=it+1, cuts=len(cuts), lb=mp_obj, ip=res.objective_value(),
                routes=routes, status=res.termination.reason.name)

if __name__ == "__main__":
    run_benders()
