# -*- coding: utf-8 -*-
"""01_direct 原型：3-index CP-SAT 直接建模（车辆自环=未使用；每客户恰被一车访问）。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import math, time
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"
SCALE = 10000
M = 10**12

def parse(path=DATA):
    lines = open(path).read().splitlines()
    K = cap = None
    for i, l in enumerate(lines):
        if l.strip().startswith("NUMBER"):
            s = lines[i+1].split(); K, cap = int(s[0]), int(s[1])
    custs = []
    for l in lines:
        s = l.split()
        if len(s) == 7 and s[0].isdigit():
            no, x, y, dem, ready, due, svc = map(int, s)
            custs.append((no, x, y, dem, ready, due, svc))
    custs.sort(key=lambda c: c[0])
    return K, cap, custs[0], custs[1:]

KMAX, cap, depot, cs = parse()
n = len(cs)
xc = [depot[1]] + [c[1] for c in cs]
yc = [depot[2]] + [c[2] for c in cs]
dem = [0] + [c[3] for c in cs]
ready = [0] + [c[4] for c in cs]
due = [0] + [c[5] for c in cs]
svc = [0] + [c[6] for c in cs]
depot_due = depot[5]
total_dem = sum(dem[1:])
LB_K = math.ceil(total_dem / cap)
def dist(i, j):
    return math.hypot(xc[i]-xc[j], yc[i]-yc[j])
d_scaled = [[0]*(n+1) for _ in range(n+1)]
for i in range(n+1):
    for j in range(n+1):
        if i != j:
            d_scaled[i][j] = int(round(dist(i, j)*SCALE))
print(f"n={n} KMAX={KMAX} cap={cap} 总需求={total_dem} LB_K={LB_K}")

def solve(Kv, time_limit=120.0):
    model = cp_model.CpModel()
    x = [[[model.NewBoolVar(f"x{k}_{i}_{j}") for j in range(n+1)] for i in range(n+1)] for k in range(Kv)]
    t = [[model.NewIntVar(0, depot_due*SCALE, f"t{k}_{i}") for i in range(n+1)] for k in range(Kv)]
    q = [[model.NewIntVar(0, cap, f"q{k}_{i}") for i in range(n+1)] for k in range(Kv)]
    for k in range(Kv):
        model.AddCircuit([(i, j, x[k][i][j]) for i in range(n+1) for j in range(n+1)])
    for i in range(1, n+1):                                  # 每客户恰被一车访问（其余车自环）
        model.Add(sum(x[k][i][i] for k in range(Kv)) == Kv - 1)
    for k in range(Kv-1):                                    # 未使用车辆=后缀（仓库自环）
        model.Add(x[k][0][0] <= x[k+1][0][0])
    for k in range(Kv):
        model.Add(t[k][0] == 0)
        model.Add(q[k][0] == 0)
        for i in range(1, n+1):
            model.Add(t[k][i] >= ready[i]*SCALE - M*x[k][i][i])
            model.Add(t[k][i] <= due[i]*SCALE + M*x[k][i][i])
            model.Add(t[k][i] <= M*(1 - x[k][i][i]))
            model.Add(q[k][i] >= dem[i] - M*x[k][i][i])
            model.Add(q[k][i] <= M*(1 - x[k][i][i]))
        for i in range(n+1):
            for j in range(1, n+1):
                if i == j:
                    continue
                model.Add(t[k][j] >= t[k][i] + svc[i]*SCALE + d_scaled[i][j] - M*(1 - x[k][i][j]))
                model.Add(q[k][j] >= q[k][i] + dem[j] - M*(1 - x[k][i][j]))
        for i in range(1, n+1):
            model.Add(t[k][i] + svc[i]*SCALE + d_scaled[i][0] <= depot_due*SCALE + M*(1 - x[k][i][0]))
    model.Minimize(sum(d_scaled[i][j]*x[k][i][j] for k in range(Kv) for i in range(n+1) for j in range(n+1) if i != j))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    t0 = time.time()
    st = solver.Solve(model)
    el = time.time()-t0
    obj = solver.ObjectiveValue() if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    routes = []
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        for k in range(Kv):
            if solver.Value(x[k][0][0]):
                continue
            route = []
            cur = 0
            while True:
                nxt = None
                for j in range(n+1):
                    if solver.Value(x[k][cur][j]):
                        nxt = j
                        break
                if nxt is None or nxt == 0:
                    break
                route.append(nxt)
                cur = nxt
                if len(route) > n:
                    break
            if route:
                routes.append(tuple(route))
    return st, obj, routes, el, solver.BestObjectiveBound()

for Kv in (2, 3):
    st, obj, routes, el, bound = solve(Kv, 120)
    print(f"K={Kv}: {st} | obj_scaled={obj} | 时间 {round(el,1)}s | bound={bound} | 路线 {routes}")
    if st == cp_model.OPTIMAL:
        exact = 0.0
        for r in routes:
            seq = [0] + list(r) + [0]
            exact += sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1))
        print(f"  精确总距离 {exact:.10f} | 两位小数 {round(exact,2)} | BKS 191.81 match: {abs(round(exact,2)-191.81)<1e-9}")
