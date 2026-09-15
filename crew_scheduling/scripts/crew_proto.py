# -*- coding: utf-8 -*-
"""csp50 机组排班核心原型：数据解析 + 完整池枚举 + 池 IP(直接基准) + CG LP 下界。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import math, time, datetime
from ortools.math_opt.python import mathopt
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt"
lines = open(DATA).read().splitlines()
N, T = map(int, lines[0].split())
tasks = [None]
for line in lines[1:1+N]:
    s, f = map(int, line.split()); tasks.append((s, f))
arcs = {}
for line in lines[1+N:]:
    i, j, c = map(int, line.split()); arcs[(i, j)] = c
print(f"N={N} T={T} arcs={len(arcs)}")

# 完整池枚举：从每个任务起 DFS，跨度 <= T
adj = {}
for (i, j), c in arcs.items():
    adj.setdefault(i, []).append((j, c))
paths = []   # (tuple tasks, cost, span)
for start in range(1, N+1):
    s0 = tasks[start][0]
    stack = [(start, [start], 0)]
    while stack:
        u, seq, c = stack.pop()
        paths.append((tuple(seq), c, tasks[u][1]-s0))
        for v, cv in adj.get(u, []):
            c2 = c + cv
            span = tasks[v][1] - s0
            if span <= T:
                stack.append((v, seq+[v], c2))
P = len(paths)
print(f"完整池: {P} 条可行路径")
pseq, pcost, pspan = zip(*paths)
pmask = []
for s2 in pseq:
    m = 0
    for i in s2:
        m |= (1 << i)
    pmask.append(m)
print("最长链:", max(len(s) for s in pseq), "| 单任务:", sum(1 for s in pseq if len(s) == 1),
      "| 双任务:", sum(1 for s in pseq if len(s) == 2), "| 三任务:", sum(1 for s in pseq if len(s) == 3))

# 下界
total_dur = sum(f-s for s, f in tasks[1:])
ev = []
for s, f in tasks[1:]:
    ev.append((s, 1)); ev.append((f, -1))
ev.sort(key=lambda x: (x[0], -x[1]))
cur = overlap = 0
for t2, d in ev:
    cur += d; overlap = max(overlap, cur)
inc = set(); outg = set()
has_in = set(); has_out = set()
for (i, j), c in arcs.items():
    has_in.add(j); has_out.add(i)
noin = [i for i in range(1, N+1) if i not in has_in]
noout = [i for i in range(1, N+1) if i not in has_out]
print(f"下界: 时长={math.ceil(total_dur/T)}, 重叠={overlap}, 无入弧={len(noin)}, 无出弧={len(noout)}")

def pool_ip(Kv, solver_name, time_limit=60.0, use_cpsat=False):
    if not use_cpsat:
        m = mathopt.Model(name="ip")
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
        return res.termination.reason, obj, wt
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
        return s.StatusName(st), (s.ObjectiveValue() if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None), time.time()-t0

print()
print("== 直接基准：完整池集合覆盖 IP（HIGHS）==")
for Kv in (26, 27):
    term, obj, wt = pool_ip(Kv, "HIGHS")
    print(f"K={Kv}: {term} obj={obj} wall={round(wt,3)}s")
    if term == mathopt.TerminationReason.INFEASIBLE:
        print(f"  -> K={Kv} 不可行，最少 crew 数 = {Kv+1}")

print()
print("== 列生成 LP 下界（GLOP + 池定价）==")
sel = [p for p in range(P) if len(pseq[p]) == 1]
sel_set = set(sel)
lp_obj = None
M_DUMMY = 10**6
for it in range(200):
    m = mathopt.Model()
    x = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in sel]
    yv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(1, N+1)]
    covers = []
    for i in range(1, N+1):
        covers.append(m.add_linear_constraint(mathopt.fast_sum([x[k] for k, p in enumerate(sel) if (pmask[p] >> i) & 1]) + yv[i-1] >= 1.0, name=f"c{i}"))
    veh = m.add_linear_constraint(mathopt.fast_sum(x) <= 27, name="veh")
    m.minimize(mathopt.fast_sum([pcost[p]*x[k] for k, p in enumerate(sel)]) + M_DUMMY*mathopt.fast_sum(yv))
    res = mathopt.solve(m, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
    lp_obj = res.objective_value()
    dv = res.dual_values()
    pi = [0.0]*(N+1)
    for i in range(1, N+1):
        pi[i] = max(0.0, dv[covers[i-1]])
    mu = dv[veh]
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
        if rc < -1e-7:
            add.append((rc, p))
    add.sort()
    add = [p for _, p in add[:1000]]
    if not add:
        print(f"CG 收敛: {it} 轮 | LP 下界 = {round(lp_obj, 6)} | 列数 {len(sel)}")
        break
    for p in add:
        if p not in sel_set:
            sel_set.add(p)
            sel.append(p)
else:
    print("未收敛")
print(f"LP 下界 vs 整数最优 3139: gap = {round((3139-lp_obj)/3139*100, 4)}%")
