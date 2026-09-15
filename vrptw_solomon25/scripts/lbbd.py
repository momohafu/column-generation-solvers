# -*- coding: utf-8 -*-
"""VRPTW c101 LBBD（逻辑 Benders 分解）：
阶段0：车辆数下界 LB_K + 角度扫描 warm-start；阶段1：分配主问题 + TSP-TW 子问题 + 逻辑割迭代。
run_lbbd() 返回结果 dict；直接运行时执行完整流程。
"""
import sys, math, time
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"
SCALE = 1000
BIG = 10**9
MAX_ITER = 60
WALL = 110.0

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

K, cap, depot, cs = parse()
n = len(cs)
xc = [depot[1]] + [c[1] for c in cs]
yc = [depot[2]] + [c[2] for c in cs]
dem = [0] + [c[3] for c in cs]
ready = [0] + [c[4] for c in cs]
due = [0] + [c[5] for c in cs]
svc = [0] + [c[6] for c in cs]
depot_due = depot[5]
total_dem = sum(dem[1:])
def dist(i, j):
    return math.hypot(xc[i]-xc[j], yc[i]-yc[j])
d_scaled = [[0]*(n+1) for _ in range(n+1)]
for i in range(n+1):
    for j in range(n+1):
        if i != j:
            d_scaled[i][j] = int(round(dist(i, j)*SCALE))

_cache = {}
def solve_route(S, time_limit=5.0):
    """单车辆 TSP-TW 子问题。返回 (status, route, cost)。带缓存。"""
    S = tuple(sorted(S))
    if S in _cache:
        return _cache[S]
    if len(S) == 0:
        return ("OPTIMAL", (), 0.0)
    if len(S) == 1:
        i = S[0]
        st = max(ready[i], dist(0, i))
        if st > due[i] or st + svc[i] + dist(i, 0) > depot_due + 1e-9:
            return ("INFEASIBLE", None, None)
        return ("OPTIMAL", (i,), 2*dist(0, i))
    nodes = [0] + list(S)
    m_ = len(nodes)
    model = cp_model.CpModel()
    xv = [[model.NewBoolVar(f"x{a}_{b}") for b in range(m_)] for a in range(m_)]
    model.AddCircuit([(a, b, xv[a][b]) for a in range(m_) for b in range(m_) if a != b])
    T = depot_due * SCALE
    t = [model.NewIntVar(0, T, f"t{a}") for a in range(m_)]
    model.Add(t[0] == 0)
    q = [model.NewIntVar(0, cap, f"q{a}") for a in range(m_)]
    model.Add(q[0] == 0)
    for a in range(m_):
        for b in range(1, m_):
            if a == b:
                continue
            j = nodes[b]
            model.Add(t[b] >= t[a] + svc[nodes[a]]*SCALE + d_scaled[nodes[a]][j] - BIG*(1 - xv[a][b]))
            model.Add(q[b] >= q[a] + dem[j] - BIG*(1 - xv[a][b]))
    for a in range(1, m_):
        i = nodes[a]
        model.Add(t[a] >= ready[i]*SCALE)
        model.Add(t[a] <= due[i]*SCALE)
        model.Add(t[a] + svc[i]*SCALE + d_scaled[i][0] <= depot_due*SCALE)
    model.Minimize(sum(d_scaled[nodes[a]][nodes[b]] * xv[a][b]
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
        for b in range(1, m_):
            if solver.Value(xv[cur][b]):
                nxt = b
                break
        if nxt is None:
            break
        route.append(nodes[nxt])
        cur = nxt
    cost = 0.0
    prev = 0
    for j in route:
        cost += dist(prev, j)
        prev = j
    cost += dist(prev, 0)
    _cache[S] = (status, tuple(route), cost)
    return _cache[S]

def tw_feasible(S):
    """快速 TW+容量可行性 DP（状态记忆 DFS），用于扫描预筛。"""
    S = sorted(S, key=lambda i: due[i])
    if sum(dem[i] for i in S) > cap:
        return False
    memo = {}
    def dfs(last, mask, t):
        key = (last, mask)
        if key in memo:
            if t >= memo[key]:
                return False
        memo[key] = t
        if mask == (1 << len(S)) - 1:
            return t + svc[last] + dist(last, 0) <= depot_due + 1e-9
        for pos, j in enumerate(S):
            if mask & (1 << pos):
                continue
            arr = max(ready[j], t + svc[last] + dist(last, j))
            if arr <= due[j] and dfs(j, mask | (1 << pos), arr):
                return True
        return False
    for pos, j in enumerate(S):
        arr = max(ready[j], 0 + svc[0] + dist(0, j))
        if arr <= due[j] and dfs(j, 1 << pos, arr):
            return True
    return False

def shrink(S):
    """贪心收缩到最小不可行核心。"""
    S = list(S)
    while True:
        for i in list(S):
            if len(S) <= 2:
                return tuple(S)
            S2 = tuple(j for j in S if j != i)
            if solve_route(S2)[0] == "INFEASIBLE":
                S.remove(i)
                break
        else:
            return tuple(S)

def build_master(nogoods, hooker, exact_cuts, fix_v, theta_ub):
    model = cp_model.CpModel()
    a = [[model.NewBoolVar(f"a{i}_{k}") for k in range(K)] for i in range(1, n+1)]
    v = [model.NewBoolVar(f"v{k}") for k in range(K)]
    for i in range(1, n+1):
        model.Add(sum(a[i-1][k] for k in range(K)) == 1)
    for k in range(K):
        for i in range(1, n+1):
            model.Add(a[i-1][k] <= v[k])
        if k + 1 < K:
            model.Add(v[k] >= v[k+1])
        model.Add(sum(dem[i] * a[i-1][k] for i in range(1, n+1)) <= cap)
    for (k, S) in nogoods:
        model.Add(sum(1 - a[i-1][k] for i in S) >= 1)
    model.Add(sum(v) == fix_v)
    theta = model.NewIntVar(0, 10**7, "theta")
    for (cks, deltas) in hooker:
        model.Add(theta >= sum(cks[k] - sum(deltas[k][i] * (1 - a[i-1][k]) for i in deltas[k]) for k in cks))
    for (cc, S_by_k) in exact_cuts:
        model.Add(theta >= cc * (1 - sum(sum(1 - a[i-1][k] for i in S) for k, S in S_by_k.items())))
    if theta_ub is not None:
        model.Add(theta <= theta_ub)
    return model, a, v, theta

def extract_assignment(solver, a):
    S_by_k = {}
    for k in range(K):
        S = tuple(i for i in range(1, n+1) if solver.Value(a[i-1][k]))
        if S:
            S_by_k[k] = S
    return S_by_k

def make_cuts(S_by_k):
    cks = {}
    deltas = {}
    C = 0.0
    for k, S in S_by_k.items():
        ck = solve_route(S)[2]
        C += ck
        cks[k] = int(round(ck * 100))
        deltas[k] = {}
        for i in S:
            S2 = tuple(j for j in S if j != i)
            st3, _, c2 = solve_route(S2)
            dlt = ck - c2 if st3 == "OPTIMAL" else ck
            deltas[k][i] = int(round(dlt * 100))
    return C, cks, deltas

def run_lbbd(verbose=True):
    """完整 LBBD 流程。返回结果 dict。"""
    LB_K = math.ceil(total_dem / cap)
    if verbose:
        print(f"车辆数下界 LB_K = ceil({total_dem}/{cap}) = {LB_K}")
    t0 = time.time()
    ang = sorted(range(1, n+1), key=lambda i: math.atan2(yc[i]-yc[0], xc[i]-xc[0]))
    cum = [0]*(n+1)
    for t_, i in enumerate(ang):
        cum[t_+1] = cum[t_] + dem[i]
    best = None
    nfeas = 0
    nchecked = 0
    for p in range(1, n-1):
        if cum[p] > cap:
            break
        for q in range(p+1, n):
            if cum[q] - cum[p] > cap:
                break
            if total_dem - cum[q] > cap:
                continue
            g1, g2, g3 = tuple(ang[:p]), tuple(ang[p:q]), tuple(ang[q:])
            nchecked += 1
            if not (tw_feasible(g1) and tw_feasible(g2) and tw_feasible(g3)):
                continue
            sts = [solve_route(g) for g in (g1, g2, g3)]
            if all(s[0] == "OPTIMAL" for s in sts):
                c = sum(s[2] for s in sts)
                nfeas += 1
                if best is None or c < best[0]:
                    best = (c, {0: g1, 1: g2, 2: g3})
    t_sweep = time.time() - t0
    if verbose:
        print(f"角度扫描: 容量可行划分 {nchecked} 个, TW预筛通过且子问题可行 {nfeas} 个, 耗时 {round(t_sweep, 1)}s")
    assert best is not None, "扫描未找到可行解"
    UB = best[0]
    UB_cents = int(round(UB * 100))
    best_routes = {k: solve_route(S)[1] for k, S in best[1].items()}
    if verbose:
        print(f"初始可行解: 距离 {round(UB, 4)}（两位小数 {round(UB, 2)}）")
        print("== 阶段 1：LBBD 最优性割迭代 ==")
    nogoods = []
    hooker = []
    exact_cuts = []
    C0, cks0, d0 = make_cuts(best[1])
    exact_cuts.append((int(round(C0*100)), best[1]))
    hooker.append((cks0, d0))
    wall0 = time.time()
    proven = False
    for it in range(MAX_ITER):
        model, a, v, theta = build_master(nogoods, hooker, exact_cuts, fix_v=LB_K, theta_ub=UB_cents - 1)
        for k, S in best[1].items():
            for i in S:
                model.AddHint(a[i-1][k], 1)
        model.Minimize(theta)
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = 20
        solver.parameters.num_search_workers = 8
        st = solver.Solve(model)
        if st == cp_model.INFEASIBLE:
            if verbose:
                print(f"iter {it+1:2d}: 主问题不可行（theta <= {UB_cents-1} 无解）-> UB={UB_cents} 分已证明最优")
            proven = True
            break
        assert st == cp_model.OPTIMAL
        th = solver.Value(theta)
        S_by_k = extract_assignment(solver, a)
        infeas = []
        costs_k = {}
        for k, S in S_by_k.items():
            st2, route, cost = solve_route(S)
            if st2 == "INFEASIBLE":
                infeas.append((k, S))
            else:
                costs_k[k] = cost
        if infeas:
            for (k, S) in infeas:
                core = shrink(S)
                nogoods.append((k, core))
            if verbose:
                print(f"iter {it+1:2d}: theta={th} | {len(infeas)} 车 TW 不可行 -> no-good 核心割 (累计 {len(nogoods)})")
            continue
        C = sum(costs_k.values())
        C_cents = int(round(C * 100))
        if C < UB:
            UB = C
            UB_cents = C_cents
            best_routes = {k: solve_route(S)[1] for k, S in S_by_k.items()}
            if verbose:
                print(f"iter {it+1:2d}: theta={th} | 可行 C={round(C, 4)} < UB -> 更新 UB={round(UB, 4)}")
        else:
            if verbose:
                print(f"iter {it+1:2d}: theta={th} | 可行 C={round(C, 4)} >= UB（不改善）")
        _, cks, dlt = make_cuts(S_by_k)
        hooker.append((cks, dlt))
        exact_cuts.append((C_cents, S_by_k))
        if th >= UB_cents - 1:
            if verbose:
                print("theta 达到 UB -> 收敛")
            break
        if time.time() - wall0 > WALL:
            if verbose:
                print("达到墙钟上限")
            break
    wall = time.time() - wall0
    routes = {k: tuple(r) for k, r in sorted(best_routes.items())}
    if verbose:
        print("== LBBD 结果 ==")
        print("LB_K =", LB_K, "| UB =", round(UB, 6), "| UB_cents =", UB_cents)
        print("no-good 割:", len(nogoods), "| Hooker 割:", len(hooker), "| 精确界割:", len(exact_cuts),
              "| 子问题缓存:", len(_cache))
        print("证明状态:", "主问题不可行 => 已证明最优" if proven else "达到迭代/时间上限（未证明）")
        print("阶段1 墙钟:", round(wall, 2), "s")
    return dict(LB_K=LB_K, sweep_feas=nfeas, sweep_checked=nchecked, sweep_time=t_sweep,
                iterations=it+1, nogoods=len(nogoods), hooker=len(hooker), exact_cuts=len(exact_cuts),
                proven=proven, ub=UB, ub_cents=UB_cents, wall=wall, routes=routes, cache=len(_cache))

if __name__ == "__main__":
    res = run_lbbd()
    exact = 0.0
    for k, r in res["routes"].items():
        seq = [0] + list(r) + [0]
        c = sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1))
        exact += c
        print("  车", k, "路线", r, "距离", round(c, 6))
    print("精确总距离", round(exact, 6), "| 两位小数", round(exact, 2),
          "| BKS 191.81 match:", abs(round(exact, 2) - 191.81) < 1e-9)
