# -*- coding: utf-8 -*-
"""04_lagrangian 原型：松弛覆盖约束的拉格朗日对偶 + 次梯度 + 贪心修复 + MIP 修复。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, math, time
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg_cpsat import build_data, enumerate_pool, solve_rmp
from ortools.sat.python import cp_model

n, xc, yc, dem, ready, due, svc, cap, depot_due, dist, d_scaled = build_data()
K = 25
LB_K = 3
EPS = 1e-7

t0 = time.time()
full_paths, full_costs, full_masks, _ = enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due)
# 候选列集 = 池扫描列生成收敛池（其 LP 松弛值 = 完整池 LP = 191.813620）
# --- 池扫描 CG ---
path_to_idx = {p: i for i, p in enumerate(full_paths)}
sel = [path_to_idx[(i,)] for i in range(1, n+1)]
sel_set = set(sel)
for it in range(200):
    obj, pi, mu = solve_rmp(n, K, full_paths, full_costs, full_masks, sel)
    add = []
    for p in range(len(full_paths)):
        if p in sel_set:
            continue
        s = 0.0
        mm = full_masks[p]
        while mm:
            lb = mm & -mm
            s += pi[lb.bit_length()-1]
            mm -= lb
        rc = full_costs[p] - s - mu
        if rc < -EPS:
            add.append((rc, p))
    add.sort()
    add = [p for _, p in add[:2000]]
    if not add:
        break
    for p in add:
        if p not in sel_set:
            sel_set.add(p)
            sel.append(p)
print("CG 池:", len(sel), "列 | LP 下界", round(obj, 6), "| 池构建耗时", round(time.time()-t0, 1), "s")

paths = [full_paths[p] for p in sel]
costs = [full_costs[p] for p in sel]
masks = [full_masks[p] for p in sel]
P = len(sel)

def rc_of(p, lam):
    s = 0.0
    mm = masks[p]
    while mm:
        lb = mm & -mm
        s += lam[lb.bit_length()-1]
        mm -= lb
    return costs[p] - s

def lagrangian(lam):
    """子问题：取 rc 最小的至多 K 个负列。返回 (L, S)。"""
    neg = [(rc_of(p, lam), p) for p in range(P)]
    neg.sort()
    S = [p for rc, p in neg[:K] if rc < -EPS]
    return sum(lam[1:]) + sum(rc_of(p, lam) for p in S), S

def repair(S):
    """贪心修复：给未覆盖客户补最便宜的覆盖列（≤K 列）。"""
    S = list(S)
    covered = [False]*(n+1)
    for p in S:
        mm = masks[p]
        while mm:
            lb = mm & -mm
            covered[lb.bit_length()-1] = True
            mm -= lb
    for i in range(1, n+1):
        if covered[i]:
            continue
        best = None
        for p in range(P):
            if p in S or not ((masks[p] >> i) & 1):
                continue
            if len(S) >= K:
                break
            if best is None or costs[p] < costs[best]:
                best = p
        if best is None:
            return None
        S.append(best)
        mm = masks[best]
        while mm:
            lb = mm & -mm
            covered[lb.bit_length()-1] = True
            mm -= lb
    return S

lam = [0.0]*(n+1)
rho = 2.0
UB_anchor = sum(2*dist(0, i) for i in range(1, n+1))   # 单客户往返总成本（安全上界）
LAM_MAX = 1000.0     # λ 截断（列成本 ~250，λ 超过上界无意义且防溢出）
STEP_MAX = 100.0     # 步长截断（‖g‖ 趋零时防发散）
best_LB = -1e18
best_UB = UB_anchor
best_routes = None
no_improve = 0
wall0 = time.time()
log = []
for it in range(600):
    L, S = lagrangian(lam)
    if L > best_LB:
        best_LB = L
        no_improve = 0
    else:
        no_improve += 1
        if no_improve >= 30:
            rho /= 2.0
            no_improve = 0
    g = [0.0]*(n+1)
    cnt = [0]*(n+1)
    for p in S:
        mm = masks[p]
        while mm:
            lb = mm & -mm
            cnt[lb.bit_length()-1] += 1
            mm -= lb
    for i in range(1, n+1):
        g[i] = 1.0 - float(cnt[i])      # 正确次梯度：1 - 覆盖次数（而非布尔标志）
    S2 = repair(S)
    if S2 is not None:
        c = sum(costs[p] for p in S2)
        if c < best_UB:
            best_UB = c
            best_routes = [paths[p] for p in S2]
    gnorm2 = sum(g[i]*g[i] for i in range(1, n+1))
    if gnorm2 <= 1e-12:
        step = 0.0
    else:
        step = min(rho * (UB_anchor - L) / gnorm2, STEP_MAX)
    for i in range(1, n+1):
        lam[i] = min(LAM_MAX, max(0.0, lam[i] + step * g[i]))
    if it % 50 == 0 or it == 599:
        gap = (best_UB-best_LB)/best_LB*100 if best_LB > 1e-9 else float('nan')
        print(f"iter {it+1:4d}: L={L:10.4f} best_LB={best_LB:10.4f} best_UB={best_UB:10.4f} gap={gap:.4f}% rho={rho:.4f} lam_max={max(lam):.1f}")
    if rho < 1e-5 or time.time()-wall0 > 110:
        print("停机 at iter", it+1)
        break

# 最终 MIP 修复（CG 池 CP-SAT）
model = cp_model.CpModel()
yv = [model.NewBoolVar(f"y{p}") for p in range(P)]
for i in range(1, n+1):
    model.Add(sum(yv[p] for p in range(P) if (masks[p] >> i) & 1) >= 1)
model.Add(sum(yv) <= K)
model.Minimize(sum(int(round(costs[p]*100))*yv[p] for p in range(P)))
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = 120
solver.parameters.num_search_workers = 8
st = solver.Solve(model)
routes = [paths[p] for p in range(P) if solver.Value(yv[p])]
exact = 0.0
for r in routes:
    seq = [0] + list(r) + [0]
    exact += sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1))
print()
print("== 拉格朗日结果 ==")
print("对偶下界 best_LB =", round(best_LB, 6))
print("贪心修复 best_UB =", round(best_UB, 6))
print("MIP 修复:", solver.StatusName(st), "目标", round(exact, 6), "| 车辆", len(routes))
for r in routes:
    print("  路线", r)
print("对偶 gap =", round((exact - best_LB)/best_LB*100, 6), "%")
print("与 BKS 191.813620 一致:", abs(exact - 191.813620) < 1e-6, "| 两位小数 191.81:", abs(round(exact,2)-191.81) < 1e-9)
