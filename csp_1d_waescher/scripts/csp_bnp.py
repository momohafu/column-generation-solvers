# -*- coding: utf-8 -*-
"""Waescher_TEST0005 分支-定价（Gilmore-Gomory）：
固定 28 辊子树：节点 CG（numpy 背包定价，尊重 forbid pairs）→ pair 分支（禁用 vs 合并）。
节点 LP=28 整数 => 28 可行(最优)；全部节点不可行 => 最优 29。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import numpy as np, math, time, datetime
from ortools.math_opt.python import mathopt

p = "/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt"
lines = open(p, encoding="utf-8").read().splitlines()
m0 = int(lines[0].strip()); L = int(lines[1].strip())
ITEMS0 = []
for line in lines[2:2+m0]:
    a = line.split(); ITEMS0.append((int(a[0]), int(a[1])))

def expand(types):
    out = []
    for idx, (l, d) in enumerate(types):
        for _ in range(d):
            out.append((l, idx))
    return out

def knap_value(types, pi, excl_types=frozenset()):
    """0/1 背包最优值（向量化，不重建）。"""
    exp = expand(types)
    if not exp:
        return 0.0
    W = np.array([e[0] for e in exp], dtype=np.int64)
    vv = np.array([pi[e[1]] for e in exp])
    dp = np.zeros(L+1, dtype=np.float64)
    for pos in range(len(exp)):
        if exp[pos][1] in excl_types:
            continue
        w = W[pos]
        cand = dp[:L+1-w] + vv[pos]
        better = cand > dp[w:]
        dp[w:][better] = cand[better]
    return float(dp.max())

def knap_rebuild(types, pi, excl_types=frozenset()):
    """背包 + 快照重建模式。返回 (值, 各类型计数)。"""
    exp = expand(types)
    W = np.array([e[0] for e in exp], dtype=np.int64)
    vv = np.array([pi[e[1]] for e in exp])
    dp = np.zeros(L+1, dtype=np.float64)
    snaps = []
    order = []
    for pos in range(len(exp)):
        if exp[pos][1] in excl_types:
            continue
        snaps.append(dp.copy())
        order.append(pos)
        w = W[pos]
        cand = dp[:L+1-w] + vv[pos]
        better = cand > dp[w:]
        dp[w:][better] = cand[better]
    if not order:
        return 0.0, [0]*len(types)
    best_cap = int(np.argmax(dp))
    best_val = float(dp[best_cap])
    counts = [0]*len(types)
    cap = best_cap
    cur = dp
    for k in range(len(order)-1, -1, -1):
        pos = order[k]
        w = W[pos]
        if cap >= w and cur[cap] == snaps[k][cap-w] + vv[pos]:
            counts[exp[pos][1]] += 1
            cur = snaps[k]
            cap -= w
    return best_val, counts

def pricing(types, pi, forbids):
    """尊重 forbid pairs 的定价：枚举每对排除哪个端点（2^F）。返回 (最优值, 最优排除集)。"""
    F = list(forbids)
    if len(F) > 8:
        raise RuntimeError("forbids 过多")
    best_val = -1e300
    best_excl = frozenset()
    for bits in range(1 << len(F)):
        excl = set()
        for k, (a, b) in enumerate(F):
            if bits >> k & 1:
                excl.add(a)
            else:
                excl.add(b)
        v = knap_value(types, pi, frozenset(excl))
        if v > best_val:
            best_val = v
            best_excl = frozenset(excl)
    return best_val, best_excl

def solve_rmp(types, patterns, K_lb, K_ub):
    nt = len(types)
    m = mathopt.Model()
    x = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in range(len(patterns))]
    yv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(nt)]
    covers = []
    for i in range(nt):
        covers.append(m.add_linear_constraint(
            mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(patterns))]) + yv[i] >= types[i][1], name=f"c{i}"))
    ubc = m.add_linear_constraint(mathopt.fast_sum(x) <= K_ub, name="ub")
    lbc = m.add_linear_constraint(mathopt.fast_sum(x) >= K_lb, name="lb")
    m.minimize(mathopt.fast_sum(x) + 10**6*mathopt.fast_sum(yv))
    res = mathopt.solve(m, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
    dv = res.dual_values()
    pi = [max(0.0, dv[covers[i]]) for i in range(nt)]
    mu_ub = dv[ubc]
    mu_lb = dv[lbc]
    xvals = {k: res.variable_values()[x[k]] for k in range(len(patterns))}
    y_used = sum(res.variable_values()[yv[i]] for i in range(nt))
    return res.objective_value(), pi, mu_ub, mu_lb, xvals, y_used

def cg_node(types, forbids, K_lb, K_ub, time_limit=120.0):
    nt = len(types)
    patterns = []
    for i, (l, d) in enumerate(types):
        a = [0]*nt
        a[i] = min(d, L//l)
        patterns.append(a)
    t0 = time.time()
    for it in range(400):
        obj, pi, mu_ub, mu_lb, xvals, y_used = solve_rmp(types, patterns, K_lb, K_ub)
        mu_sum = mu_ub + mu_lb
        best_val, best_excl = pricing(types, pi, forbids)
        rc = 1.0 - best_val - mu_sum
        if rc >= -1e-7:
            break
        _, pat = knap_rebuild(types, pi, best_excl)
        patterns.append(pat)
        if time.time()-t0 > time_limit:
            break
    infeasible = (y_used > 1e-6)
    integral = (not infeasible) and all(abs(v-round(v)) < 1e-6 for v in xvals.values())
    return dict(lp=None if infeasible else obj, infeasible=infeasible, integral=integral,
                xvals=xvals, patterns=patterns, iters=it+1, time=time.time()-t0)

def together_flow(types, xvals, patterns):
    """f[i][j] = Σ_p x_p [a_ip>0 ∧ a_jp>0]。"""
    nt = len(types)
    f = {}
    for k, val in xvals.items():
        if val <= 1e-6:
            continue
        pat = patterns[k]
        ones = [i for i in range(nt) if pat[i] > 0]
        for a in range(len(ones)):
            for b in range(a+1, len(ones)):
                key = (ones[a], ones[b])
                f[key] = f.get(key, 0.0) + val
    return f

def branch_and_price(time_limit=280.0, node_limit=200):
    """固定 K=28 的 GG 搜索。返回 (最优辊数, 证明状态, 统计)。"""
    root = (tuple((l, d) for l, d in ITEMS0), frozenset(), 28, 28)
    stack = [root]
    wall0 = time.time()
    nodes = 0
    branches = 0
    found28 = None
    while stack and time.time()-wall0 < time_limit and nodes < node_limit:
        types, forbids, K_lb, K_ub = stack.pop()
        nodes += 1
        res = cg_node(types, forbids, K_lb, K_ub)
        tag = f"node{nodes} F={len(forbids)}"
        if res["infeasible"]:
            print(f"{tag}: LP 不可行 -> 剪枝 | CG {res['iters']} 轮 {res['time']:.1f}s")
            continue
        lp = res["lp"]
        if res["integral"]:
            print(f"{tag}: LP={lp:.6f} 整数解 -> 28 辊可行！")
            found28 = True
            break
        # 分支：挑 together flow 分数部分最接近 0.5 的 pair
        f = together_flow(types, res["xvals"], res["patterns"])
        cands = []
        for (a, b), val in f.items():
            frac = val - math.floor(val)
            if 1e-6 < frac < 1-1e-6:
                cands.append((min(frac, 1-frac), a, b, val))
        if not cands:
            print(f"{tag}: LP={lp:.6f} 分数但无分数 pair（数值）-> 停止")
            break
        cands.sort(reverse=True)
        _, a, b, val = cands[0]
        branches += 1
        print(f"{tag}: LP={lp:.6f} 分数 -> pair 分支 ({a},{b}) 流量 {val:.3f}")
        # 合并子节点：类型 a,b 合并（长度和 <= L，需求 min），残余保留
        types_list = list(types)
        la, da = types_list[a]; lb, db = types_list[b]
        if la + lb <= L:
            dm = min(da, db)
            new_types = []
            for idx, (l, d) in enumerate(types_list):
                if idx == a or idx == b:
                    continue
                new_types.append((l, d))
            new_types.append((la + lb, dm))
            if da > db:
                new_types.append((la, da - db))
            elif db > da:
                new_types.append((lb, db - da))
            merged = (tuple(new_types), forbids, K_lb, K_ub)
            stack.append(merged)
        # 禁用子节点
        forb = (tuple(types), forbids | frozenset([(a, b)]), K_lb, K_ub)
        stack.append(forb)
    wall = time.time()-wall0
    print("== B&P 汇总 ==")
    print("节点:", nodes, "| 分支:", branches, "| 墙钟:", round(wall, 1), "s")
    if found28:
        print("结论: 28 辊可行 -> 最优 28（下界 28）")
        return 28, "feasible", nodes, branches, wall
    if nodes >= node_limit or wall >= time_limit:
        print("结论: 达到上限，未决（需更多时间）")
        return None, "limit", nodes, branches, wall
    print("结论: 28 辊子树全部不可行 -> 最优 29（FFD/CP-SAT 已给 29）")
    return 29, "proved", nodes, branches, wall

if __name__ == "__main__":
    branch_and_price()
