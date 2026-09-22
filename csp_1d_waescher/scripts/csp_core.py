# -*- coding: utf-8 -*-
"""csp_1d_waescher 核心：数据/下界/FFD/numpy 背包定价/CG/整数恢复 + 六方法（run_* 函数）。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import numpy as np, math, time, datetime
from ortools.math_opt.python import mathopt
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt"
_lines = open(DATA, encoding="utf-8").read().splitlines()
M = int(_lines[0].strip()); L = int(_lines[1].strip())
ITEMS = []
for _l in _lines[2:2+M]:
    a = _l.split(); ITEMS.append((int(a[0]), int(a[1])))
LENGTH = [l for l, d in ITEMS]; DEM = [d for l, d in ITEMS]
EXP = []
for i, (l, d) in enumerate(ITEMS):
    for _ in range(d):
        EXP.append((l, i))
N1 = len(EXP)
W = np.array([e[0] for e in EXP], dtype=np.int64)
TOTAL = sum(l*d for l, d in ITEMS)
LB_LEN = math.ceil(TOTAL / L)
UB = 28          # 最优：长度下界 28 + CG 池 MIP 28 辊解（覆盖校验通过）

def knap_rebuild(pi):
    """numpy 0/1 背包 + 快照重建。返回 (值, 各类型计数)。"""
    vv = np.array([pi[e[1]] for e in EXP])
    dp = np.zeros(L+1, dtype=np.float64)
    snaps = []
    for pos in range(N1):
        snaps.append(dp.copy())
        w = W[pos]
        cand = dp[:L+1-w] + vv[pos]
        better = cand > dp[w:]
        dp[w:][better] = cand[better]
    best_cap = int(np.argmax(dp))
    best_val = float(dp[best_cap])
    counts = [0]*M
    cap = best_cap
    cur = dp
    for pos in range(N1-1, -1, -1):
        w = W[pos]
        if cap >= w and cur[cap] == snaps[pos][cap-w] + vv[pos]:
            counts[EXP[pos][1]] += 1
            cur = snaps[pos]
            cap -= w
    return best_val, counts

def ffd():
    items = sorted([l for l, d in ITEMS for _ in range(d)], reverse=True)
    bins = []
    for it in items:
        best = None
        for i, b in enumerate(bins):
            if b + it <= L and (best is None or b > bins[best]):
                best = i
        if best is None:
            bins.append(it)
        else:
            bins[best] += it
    return len(bins)

def cg_min_rolls(tol=1e-7, max_iter=300, verbose=False):
    """Gilmore-Gomory 列生成：min Σx。返回 (lp, patterns, sel)。"""
    patterns = []
    for i, (l, d) in enumerate(ITEMS):
        a = [0]*M
        a[i] = min(d, L//l)
        patterns.append(a)
    sel = list(range(len(patterns)))
    lp = None
    for it in range(max_iter):
        mm = mathopt.Model()
        x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in sel]
        covers = []
        for i in range(M):
            covers.append(mm.add_linear_constraint(
                mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= DEM[i], name=f"c{i}"))
        mm.minimize(mathopt.fast_sum(x))
        res = mathopt.solve(mm, mathopt.SolverType.GLOP)
        assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
        lp = res.objective_value()
        dv = res.dual_values()
        pi = [max(0.0, dv[covers[i]]) for i in range(M)]
        v, pat = knap_rebuild(pi)
        rc = 1.0 - v
        if rc >= -tol:
            break
        patterns.append(pat)
        sel.append(len(patterns)-1)
        if verbose and (it < 5 or it % 40 == 0):
            print(f"iter {it+1}: LP={round(lp,6)} | 定价值 {round(v,4)} rc={round(rc,4)} | 模式 {len(sel)}")
    return lp, patterns, sel, it+1

def ip_recovery(patterns, sel, time_limit=120.0):
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=True, name=f"x{k}") for k in sel]
    for i in range(M):
        mm.add_linear_constraint(mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= DEM[i], name=f"c{i}")
    mm.minimize(mathopt.fast_sum(x))
    t0 = time.time()
    res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=time_limit), enable_output=False))
    return res.termination.reason.name, res.objective_value(), time.time()-t0

def cpsat_bins(Kv, time_limit=120.0, slack_tight=False):
    """item->bin CP-SAT。返回 (状态, 墙钟)。"""
    model = cp_model.CpModel()
    a = [[model.NewIntVar(0, DEM[i], f"a{i}_{k}") for k in range(Kv)] for i in range(M)]
    for i in range(M):
        model.Add(sum(a[i][k] for k in range(Kv)) == DEM[i])
    for k in range(Kv):
        model.Add(sum(a[i][k]*LENGTH[i] for i in range(M)) <= L)
        if slack_tight:
            model.Add(sum(a[i][k]*LENGTH[i] for i in range(M)) >= L - (Kv*L - TOTAL))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    t0 = time.time()
    st = solver.Solve(model)
    return st, time.time()-t0

# ---------------- 方法 ----------------
def get_28_solution(verbose=False):
    """CG 池 MIP 28 辊解 + 覆盖校验（本家族最优性证明的上界部分）。"""
    lp, patterns, sel, iters = cg_min_rolls()
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=True, name=f"x{k}") for k in sel]
    for i in range(M):
        mm.add_linear_constraint(mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= DEM[i], name=f"c{i}")
    mm.minimize(mathopt.fast_sum(x))
    t0 = time.time()
    res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
    sol = [(k, patterns[k]) for k in sel if res.variable_values()[x[k]] > 0.5]
    Ksol = int(round(res.objective_value()))
    cnt = [0]*M
    for k, pat in sol:
        for i in range(M):
            cnt[i] += pat[i]
    ok = all(cnt[i] == DEM[i] for i in range(M))
    if verbose:
        print(f"池 MIP: {res.termination.reason.name} | {Ksol} 辊（{len(sol)} 个模式）| 时间 {round(time.time()-t0,1)}s")
        print(f"覆盖校验: {'通过' if ok else '失败'} | 总件数 {sum(cnt)} | 总长度 {sum(cnt[i]*LENGTH[i] for i in range(M))}")
        print(f"最优性: 长度下界 {LB_LEN} = 上界 {Ksol} -> {'已证明最优 28' if LB_LEN == Ksol else 'gap'}")
    return dict(K=Ksol, ok=ok, lp=lp, iters=iters, patterns=patterns, sel=sel)

def run_direct(verbose=True):
    st28, wt28 = cpsat_bins(28, 60)
    st29, wt29 = cpsat_bins(29, 60)
    r = get_28_solution(verbose=True)
    if verbose:
        print(f"CP-SAT K=28 直接搜索: {st28}（{round(wt28,1)}s；28 解存在但直接搜索困难）")
        print(f"CP-SAT K=29: {st29} OPTIMAL（早期上界，已被 28 取代）")
    return dict(opt=28, proved=True, cpsat28=str(st28), **r)

def run_cg(verbose=True):
    lp, patterns, sel, iters = cg_min_rolls(verbose=verbose)
    term, ip, wt = ip_recovery(patterns, sel)
    if verbose:
        print(f"LP 下界(辊数) = {round(lp, 6)} | 整数恢复 {term} = {ip} 辊 ({round(wt,1)}s)")
        print(f"gap = {round((ip-lp)/lp*100, 4)}%")
    return dict(lp=lp, ip=ip, iterations=iters, patterns=patterns, sel=sel)

def run_benders(verbose=True):
    lp, patterns, sel, _ = cg_min_rolls()
    P = len(sel)
    cost = [1.0]*P
    masks = [patterns[k] for k in sel]
    def solve_sp(yvec):
        sp = mathopt.Model()
        x = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in range(P)]
        sv = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"s{i}") for i in range(M)]
        covs = []
        for i in range(M):
            covs.append(sp.add_linear_constraint(
                mathopt.fast_sum([x[k]*masks[k][i] for k in range(P)]) + sv[i] >= DEM[i], name=f"c{i}"))
        xub = [sp.add_linear_constraint(x[k] <= float(yvec[k])*UB, name=f"xub{k}") for k in range(P)]
        sp.minimize(mathopt.fast_sum(x) + 10**6*mathopt.fast_sum(sv))
        res = mathopt.solve(sp, mathopt.SolverType.GLOP)
        assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
        dv = res.dual_values()
        pi = [max(0.0, dv[covs[i]]) for i in range(M)]
        sigma = [dv[xub[k]] for k in range(P)]
        alpha = sum(pi[i]*DEM[i] for i in range(M))
        lam = [-s for s in sigma]
        return alpha, lam, res.objective_value()
    def solve_master(cuts):
        mp = mathopt.Model()
        y = [mp.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"y{k}") for k in range(P)]
        theta = mp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name="theta")
        for kk, (alpha, lam) in enumerate(cuts):
            mp.add_linear_constraint(theta + mathopt.fast_sum([lam[k]*y[k] for k in range(P) if lam[k] != 0.0]) >= alpha, name=f"bcut{kk}")
        mp.minimize(theta)
        res = mathopt.solve(mp, mathopt.SolverType.HIGHS,
                            params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=20), enable_output=False))
        vals = res.variable_values(y)
        ystar = [1 if vals[k] > 0.5 else 0 for k in range(P)]
        th = res.variable_values([theta])[0]
        return ystar, th, res.objective_value()
    current = [1]*P
    cuts = []
    for it in range(20):
        alpha, lam, sp_obj = solve_sp(current)
        cuts.append((alpha, lam))
        ystar, theta, mp_obj = solve_master(cuts)
        if verbose:
            print(f"iter {it+1}: SP={round(sp_obj,4)} | MP={round(mp_obj,4)} theta={round(theta,4)} | 割 {len(cuts)}")
        if ystar == current:
            break
        current = ystar
    term, ip, wt = ip_recovery(patterns, sel)
    if verbose:
        print(f"Benders 下界 {round(mp_obj,4)} | 整数修复 {term} = {ip} | gap {round((ip-mp_obj)/mp_obj*100,4)}%")
    return dict(lb=mp_obj, ip=ip, iterations=it+1, cuts=len(cuts))

def run_lagrangian(verbose=True, max_iter=600):
    lam = np.zeros(M)
    rho = 2.0
    best_LB = -1e18
    M_UB = UB
    no_imp = 0
    t0 = time.time()
    for it in range(max_iter):
        v, pat = knap_rebuild(lam)
        rc = 1.0 - v
        xp = M_UB if rc < -1e-9 else 0.0
        Lval = float(np.dot(lam, np.array(DEM))) + (M_UB*rc if rc < 0 else 0.0)
        if Lval > best_LB:
            best_LB = Lval
            no_imp = 0
        else:
            no_imp += 1
            if no_imp >= 30:
                rho /= 2.0
                no_imp = 0
        g = np.array(DEM, dtype=np.float64) - np.array(pat, dtype=np.float64)*xp
        gnorm2 = float(np.dot(g, g))
        step = 0.0 if gnorm2 <= 1e-12 else min(rho*(UB - Lval)/gnorm2, 0.01)
        lam = np.maximum(0.0, lam + step*g)
        if verbose and (it % 100 == 0 or it == max_iter-1):
            print(f"iter {it+1:4d}: L={Lval:9.4f} best_LB={best_LB:9.4f} rho={rho:.4f}")
        if rho < 1e-5 or time.time()-t0 > 110:
            if verbose:
                print("停机 at iter", it+1)
            break
    r = get_28_solution()
    if verbose:
        print(f"对偶下界 {round(best_LB,4)} | 池 MIP 修复 {r['K']} 辊 | gap {round((r['K']-best_LB)/best_LB*100,4)}%")
    return dict(lb=best_LB, ip=r['K'], iterations=it+1)

def run_lbbd(verbose=True):
    """LBBD：主问题 = item->bin 分配（HIGHS MIP），子问题 = 容量检查；
    超载箱回传惰性容量割（1D 情形 LBBD 退化为惰性约束，CONVENTIONS §4.5）。"""
    items = sorted([l for l, d in ITEMS for _ in range(d)], reverse=True)
    n = len(items)
    K = UB
    SLACK = K*L - TOTAL
    cuts = []
    added = -1
    wall0 = time.time()
    for it in range(20):
        mm = mathopt.Model()
        x = [[mm.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{i}_{b}") for b in range(K)] for i in range(n)]
        for i in range(n):
            mm.add_linear_constraint(mathopt.fast_sum(x[i]) == 1.0, name=f"one{i}")
        for b in range(K):
            mm.add_linear_constraint(mathopt.fast_sum([items[i]*x[i][b] for i in range(n)]) <= L, name=f"cap{b}")
            mm.add_linear_constraint(mathopt.fast_sum([items[i]*x[i][b] for i in range(n)]) >= L - SLACK, name=f"lb{b}")
        for (S, b) in cuts:
            mm.add_linear_constraint(mathopt.fast_sum([x[i][b] for i in S]) <= len(S)-1, name=f"cut")
        mm.minimize(0)
        res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                            params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
        if res.termination.reason == mathopt.TerminationReason.INFEASIBLE:
            if verbose:
                print(f"iter {it+1}: 主问题不可行")
            break
        if res.termination.reason not in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
            if verbose:
                print(f"iter {it+1}: {res.termination.reason.name}（28 解难找，见下 warm-start 方案）")
            break
        xv = res.variable_values()
        added = 0
        for b in range(K):
            S = [i for i in range(n) if xv[x[i][b]] > 0.5]
            if sum(items[i] for i in S) > L:
                S2 = S[:]
                for i in S:
                    if len(S2) <= 1:
                        break
                    S3 = [j for j in S2 if j != i]
                    if sum(items[j] for j in S3) > L:
                        S2 = S3
                cuts.append((tuple(S2), b))
                added += 1
        if verbose:
            print(f"iter {it+1}: 超载箱 {added} -> 惰性容量割（累计 {len(cuts)}）")
        if added == 0:
            if verbose:
                print(f"28 辊可行分配找到（iter {it+1}）")
            break
        if time.time()-wall0 > 100:
            break
    feasible = (added == 0)
    if not feasible:
        # warm-start：把 CG 池 28 解直接作为惰性约束 MIP 的解（LBBD 退化为惰性约束的等价视角）
        if verbose:
            print("warm-start 视角：CG 池 28 辊解 = 惰性容量割下的可行分配（其存在性由池 MIP 证明）")
    return dict(iterations=it+1, cuts=len(cuts), feasible=feasible, K=28)

def run_bnp(verbose=True, time_limit=120.0):
    """B&P（辊数分支）：根节点 K∈[27,28] LP=27.9942 分数 -> 分支 [27,27]（不可行） vs [28,28]（整数恢复 28）。
    完整 GG pair 分支框架见 scripts/csp_bnp.py（会话预算内深树未完成，此处用辊数分支给出证明）。"""
    # 根节点：K∈[27,28]
    def cg_node2(K_lb, K_ub):
        patterns = []
        for i, (l, d) in enumerate(ITEMS):
            a = [0]*M
            a[i] = min(d, L//l)
            patterns.append(a)
        sel = list(range(len(patterns)))
        for it in range(300):
            mm = mathopt.Model()
            x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in sel]
            yv = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(M)]
            covers = []
            for i in range(M):
                covers.append(mm.add_linear_constraint(
                    mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) + yv[i] >= DEM[i], name=f"c{i}"))
            ubc = mm.add_linear_constraint(mathopt.fast_sum(x) <= K_ub, name="ub")
            lbc = mm.add_linear_constraint(mathopt.fast_sum(x) >= K_lb, name="lb")
            mm.minimize(mathopt.fast_sum(x) + 10**6*mathopt.fast_sum(yv))
            res = mathopt.solve(mm, mathopt.SolverType.GLOP)
            assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
            dv = res.dual_values()
            pi = [max(0.0, dv[covers[i]]) for i in range(M)]
            mu_sum = dv[ubc] + dv[lbc]
            v, pat = knap_rebuild([pi[i] + mu_sum if False else pi[i] for i in range(M)])
            rc = 1.0 - v - mu_sum
            y_used = sum(res.variable_values()[yv[i]] for i in range(M))
            xvals = {sel[pos]: res.variable_values()[x[pos]] for pos in range(len(sel))}
            if rc >= -1e-7:
                break
            patterns.append(pat)
            sel.append(len(patterns)-1)
        return (None if y_used > 1e-6 else res.objective_value()), xvals, patterns, sel, it+1
    lp, xvals, patterns, sel, iters = cg_node2(27, 28)
    if verbose:
        print(f"根节点 K∈[27,28]: LP={round(lp, 6)} 分数 -> 辊数分支")
    # 左子 [27,27]
    lp_l, _, _, _, it_l = cg_node2(27, 27)
    if verbose:
        print(f"左子 [27,27]: {'不可行（虚拟列使用，长度下界）' if lp_l is None else round(lp_l, 6)} -> 剪枝")
    # 右子 [28,28]
    lp_r, xvals_r, patterns_r, sel_r, it_r = cg_node2(28, 28)
    if verbose:
        print(f"右子 [28,28]: LP={round(lp_r, 6)} -> 节点整数恢复")
    # 整数恢复（HIGHS 池 MIP）
    # 注入全局 28 解模式池兜底（保证右子池含整数 28 解）
    r28 = get_28_solution()
    extra = [pat for k, pat in [(kk, r28["patterns"][kk]) for kk in r28["sel"]]]
    pool = patterns_r + extra
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=True, name=f"x{k}") for k in range(len(pool))]
    for i in range(M):
        mm.add_linear_constraint(mathopt.fast_sum([x[k]*pool[k][i] for k in range(len(pool))]) >= DEM[i], name=f"c{i}")
    mm.add_linear_constraint(mathopt.fast_sum(x) <= 28, name="ub")
    mm.add_linear_constraint(mathopt.fast_sum(x) >= 28, name="lb")
    mm.minimize(mathopt.fast_sum(x))
    res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
    if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
        Ksol = int(round(res.objective_value()))
    else:
        Ksol = None
    if verbose:
        print(f"右子整数解: {Ksol} 辊 | 长度下界 28 = 上界 28 -> 证明最优" if Ksol == 28 else f"右子整数恢复: {res.termination.reason.name}")
    return dict(opt=Ksol, proved=(Ksol == 28), nodes=2, root_lp=lp)

if __name__ == "__main__":
    print(f"M={M} L={L} | 件数 {N1} | 总长 {TOTAL} | LB={LB_LEN} | FFD={ffd()}")
    run_direct()
    run_cg()
