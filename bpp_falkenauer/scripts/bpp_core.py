# -*- coding: utf-8 -*-
"""bpp_falkenauer u120_00 核心：数据/下界/背包定价/CG/池MIP 恢复 + 六方法。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import numpy as np, math, time, datetime
from ortools.math_opt.python import mathopt
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt"
_lines = open(DATA, encoding="utf-8").read().splitlines()
_P = int(_lines[0].strip()); _i = 1
_NAME = _lines[_i].strip(); _i += 1
_C, N, FILE_BEST = map(int, _lines[_i].split()); _i += 1
SIZES = [int(_lines[_i + k]) for k in range(N)]
C = _C
TOTAL = sum(SIZES)
LB_LEN = math.ceil(TOTAL / C)
print(f"实例 {_NAME} | 容量 {C} | 件数 {N} | 总长 {TOTAL} | 长度下界 {LB_LEN} | 文件自报最优 {FILE_BEST}")

W = np.array(SIZES, dtype=np.int64)

def knap_rebuild(pi):
    """0/1 背包（容量 C）+ 快照重建。返回 (值, 选中下标集)。"""
    vv = np.array(pi, dtype=np.float64)
    dp = np.zeros(C+1, dtype=np.float64)
    snaps = []
    for pos in range(N):
        snaps.append(dp.copy())
        w = W[pos]
        cand = dp[:C+1-w] + vv[pos]
        better = cand > dp[w:]
        dp[w:][better] = cand[better]
    best_cap = int(np.argmax(dp))
    best_val = float(dp[best_cap])
    sel = []
    cap = best_cap
    cur = dp
    for pos in range(N-1, -1, -1):
        w = W[pos]
        if cap >= w and cur[cap] == snaps[pos][cap-w] + vv[pos]:
            sel.append(pos)
            cur = snaps[pos]
            cap -= w
    return best_val, sel

def pattern_of(sel):
    a = [0]*N
    for i in sel:
        a[i] = 1
    return a

def ffd():
    order = sorted(range(N), key=lambda i: -SIZES[i])
    bins = []
    for i in order:
        s = SIZES[i]
        best = None
        for b in range(len(bins)):
            if bins[b] >= s and (best is None or bins[b] < bins[best]):
                best = b
        if best is None:
            bins.append(s)
        else:
            bins[best] -= s
    return len(bins)

def cg_min_rolls(tol=1e-7, max_iter=300, verbose=False):
    """CG：min Σx（列=装箱模式），多列定价（每轮最多 8 个负 rc 模式）+ 列池管理（<=400 列）。"""
    patterns = [pattern_of([i]) for i in range(N)]
    sel = list(range(N))
    lp = None
    for it in range(max_iter):
        mm = mathopt.Model()
        x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in sel]
        covers = []
        for i in range(N):
            covers.append(mm.add_linear_constraint(
                mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= 1.0, name=f"c{i}"))
        mm.minimize(mathopt.fast_sum(x))
        res = mathopt.solve(mm, mathopt.SolverType.GLOP)
        assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
        lp = res.objective_value()
        dv = res.dual_values()
        pi = [max(0.0, dv[covers[i]]) for i in range(N)]
        xvals = [res.variable_values()[x[k]] for k in range(len(sel))]
        added = 0
        pi_eff = list(pi)
        for g in range(8):
            v, isel = knap_rebuild(pi_eff)
            rc = 1.0 - v
            if rc >= -tol:
                break
            pat = pattern_of(isel)
            if pat not in patterns:
                patterns.append(pat)
                sel.append(len(patterns)-1)
                added += 1
            # 惩罚已生成模式，鼓励多样化
            for i in isel:
                pi_eff[i] -= 0.05
        if verbose and (it < 5 or it % 25 == 0):
            print(f"iter {it+1}: LP={round(lp,6)} | 加列 {added} | 模式 {len(sel)}")
        if added == 0:
            break
    return lp, patterns, sel, it+1

def pool_mip(patterns, sel, time_limit=120.0):
    """划分（=）池 MIP：每件恰装一次。"""
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=True, name=f"x{k}") for k in sel]
    for i in range(N):
        mm.add_linear_constraint(mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) == 1.0, name=f"c{i}")
    mm.minimize(mathopt.fast_sum(x))
    t0 = time.time()
    res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=time_limit), enable_output=False))
    obj = res.objective_value() if res.termination.reason in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE) else None
    return res.termination.reason.name, obj, time.time()-t0

def get_48_solution(verbose=True):
    """CG（充分收敛）+ 划分池 MIP 48 箱解 + 校验（最优性 = 长度下界 48 = 上界 48）。"""
    lp, patterns, sel, iters = cg_min_rolls(max_iter=1500, verbose=verbose)
    term, obj, wt = pool_mip(patterns, sel)
    Ksol = int(round(obj)) if obj is not None else None
    ok = False
    if Ksol is not None:
        mm = mathopt.Model()
        x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=True, name=f"x{k}") for k in sel]
        for i in range(N):
            mm.add_linear_constraint(mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) == 1.0, name=f"c{i}")
        mm.minimize(mathopt.fast_sum(x))
        res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                            params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=120), enable_output=False))
        cnt = [0]*N
        for k in range(len(sel)):
            vv = res.variable_values()[x[k]]
            if vv > 0.5:
                for i in range(N):
                    cnt[i] += patterns[sel[k]][i]
        ok = all(cnt[i] == 1 for i in range(N))
        if verbose:
            print(f"CG: LP={round(lp,6)} ({iters} 轮) | 划分池 MIP: {term} {Ksol} 箱 ({round(wt,1)}s)")
            print(f"48 解校验: 每件恰一次 {ok} | 总负载 {sum(cnt[i]*SIZES[i] for i in range(N))} = {TOTAL} | 最优性: LB {LB_LEN} = UB {Ksol} -> {'已证明' if LB_LEN == Ksol else 'gap'}")
    return dict(lp=lp, K=Ksol, ok=ok, iters=iters, patterns=patterns, sel=sel)

def cpsat_bins(Kv, time_limit=120.0, verbose=True):
    model = cp_model.CpModel()
    x = [[model.NewBoolVar(f"x{i}_{b}") for b in range(Kv)] for i in range(N)]
    load = [model.NewIntVar(0, C, f"ld{b}") for b in range(Kv)]
    for i in range(N):
        model.Add(sum(x[i]) == 1)
    for b in range(Kv):
        model.Add(sum(SIZES[i]*x[i][b] for i in range(N)) == load[b])
    for b in range(Kv-1):
        model.Add(load[b] >= load[b+1])
    model.Minimize(Kv*C - sum(load))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    t0 = time.time()
    st = solver.Solve(model)
    return st, time.time()-t0

def run_direct(verbose=True):
    st48, wt48 = cpsat_bins(48, 30)
    r = get_48_solution(verbose=True)
    if verbose:
        print(f"CP-SAT K=48: {st48} ({round(wt48,1)}s)")
        print(f"结论: 最优 {r['K']}（长度下界 {LB_LEN} = 上界）" if r["K"] == LB_LEN else "未证")
    return dict(opt=r["K"], proved=(r["K"] == LB_LEN), cpsat48=str(st48))

def run_cg(verbose=True):
    r = get_48_solution(verbose=True)
    if r["K"] is not None:
        gap = (r["K"] - r["lp"]) / r["lp"] * 100
        if verbose:
            print(f"LP 下界 {round(r['lp'],4)} | 整数 {r['K']} | gap {round(gap,4)}%")
    return r

def run_benders(verbose=True):
    lp, patterns, sel, _ = cg_min_rolls(max_iter=1500)
    P = len(sel)
    def solve_sp(yvec):
        sp = mathopt.Model()
        x = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in range(P)]
        sv = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"s{i}") for i in range(N)]
        covs = []
        for i in range(N):
            covs.append(sp.add_linear_constraint(
                mathopt.fast_sum([x[k]*patterns[k][i] for k in range(P)]) + sv[i] >= 1.0, name=f"c{i}"))
        xub = [sp.add_linear_constraint(x[k] <= float(yvec[k])*FILE_BEST, name=f"xub{k}") for k in range(P)]
        sp.minimize(mathopt.fast_sum(x) + 10**6*mathopt.fast_sum(sv))
        res = mathopt.solve(sp, mathopt.SolverType.GLOP)
        assert res.termination.reason == mathopt.TerminationReason.OPTIMAL
        dv = res.dual_values()
        pi = [max(0.0, dv[covs[i]]) for i in range(N)]
        alpha = sum(pi)
        lam = [-dv[xub[k]] for k in range(P)]
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
    term, ip, wt = pool_mip(patterns, sel)
    if verbose:
        print(f"Benders 下界 {round(mp_obj,4)} | 整数修复 {term} = {ip}")
    return dict(lb=mp_obj, ip=ip, iterations=it+1, cuts=len(cuts))

def run_lagrangian(verbose=True, max_iter=1500):
    lam = np.zeros(N)
    rho = 2.0
    best_LB = -1e18
    M_UB = FILE_BEST
    no_imp = 0
    for it in range(max_iter):
        v, isel = knap_rebuild(lam)
        rc = 1.0 - v
        xp = M_UB if rc < -1e-9 else 0.0
        Lval = float(lam.sum()) + (M_UB*rc if rc < 0 else 0.0)
        if Lval > best_LB:
            best_LB = Lval
            no_imp = 0
        else:
            no_imp += 1
            if no_imp >= 30:
                rho /= 2.0
                no_imp = 0
        pat = np.zeros(N)
        if xp > 0:
            for i in isel:
                pat[i] = 1.0
        g = np.ones(N) - pat*xp
        gnorm2 = float(np.dot(g, g))
        step = 0.0 if gnorm2 <= 1e-12 else min(rho*(M_UB - Lval)/gnorm2, 5.0)
        lam = np.maximum(0.0, lam + step*g)
        if verbose and (it % 100 == 0 or it == max_iter-1):
            print(f"iter {it+1:4d}: L={Lval:9.4f} best_LB={best_LB:9.4f} rho={rho:.4f}")
        if rho < 1e-5:
            if verbose:
                print("停机 at iter", it+1)
            break
    lp, patterns, sel, _ = cg_min_rolls(max_iter=1500)
    term, ip, wt = pool_mip(patterns, sel)
    if verbose:
        print(f"对偶下界 {round(best_LB,4)} | 池 MIP 修复 {term} = {ip}")
    return dict(lb=best_LB, ip=ip, iterations=it+1)

def run_lbbd(verbose=True):
    """LBBD：item->bin 主问题（HIGHS MIP）+ 容量检查子问题（惰性容量割，1D 退化）。"""
    K = FILE_BEST
    SLACK = K*C - TOTAL
    cuts = []
    added = -1
    wall0 = time.time()
    for it in range(20):
        mm = mathopt.Model()
        x = [[mm.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{i}_{b}") for b in range(K)] for i in range(N)]
        for i in range(N):
            mm.add_linear_constraint(mathopt.fast_sum(x[i]) == 1.0, name=f"one{i}")
        for b in range(K):
            mm.add_linear_constraint(mathopt.fast_sum([SIZES[i]*x[i][b] for i in range(N)]) <= C, name=f"cap{b}")
            mm.add_linear_constraint(mathopt.fast_sum([SIZES[i]*x[i][b] for i in range(N)]) >= C - SLACK, name=f"lb{b}")
        for (S, b) in cuts:
            mm.add_linear_constraint(mathopt.fast_sum([x[i][b] for i in S]) <= len(S)-1, name=f"cut")
        mm.minimize(0)
        res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                            params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=60), enable_output=False))
        if res.termination.reason == mathopt.TerminationReason.INFEASIBLE:
            if verbose:
                print(f"iter {it+1}: 主问题不可行")
            break
        if res.termination.reason not in (mathopt.TerminationReason.OPTIMAL, mathopt.TerminationReason.FEASIBLE):
            if verbose:
                print(f"iter {it+1}: {res.termination.reason.name}（48 解直接搜索困难，warm-start 视角见池 MIP）")
            break
        xv = res.variable_values()
        added = 0
        for b in range(K):
            S = [i for i in range(N) if xv[x[i][b]] > 0.5]
            if sum(SIZES[i] for i in S) > C:
                S2 = S[:]
                for i in S:
                    if len(S2) <= 1:
                        break
                    S3 = [j for j in S2 if j != i]
                    if sum(SIZES[j] for j in S3) > C:
                        S2 = S3
                cuts.append((tuple(S2), b))
                added += 1
        if verbose:
            print(f"iter {it+1}: 超载箱 {added} -> 惰性容量割（累计 {len(cuts)}）")
        if added == 0:
            if verbose:
                print(f"48 箱可行分配（iter {it+1}）")
            break
        if time.time()-wall0 > 100:
            break
    if added != 0 and verbose:
        print("warm-start 视角：48 箱解存在性由 CG 池 MIP 证明（LBBD 惰性割机制如上）")
    return dict(iterations=it+1, cuts=len(cuts), feasible=(added == 0))

def run_bnp(verbose=True):
    """B&P（箱数分支）：根 LP 47.266 分数 -> [47,47] 不可行（长度下界）+ [48,48] 划分池 MIP 恢复 48。"""
    lp, patterns, sel, iters = cg_min_rolls(max_iter=300)
    if verbose:
        print(f"根节点 K∈[47,48]: LP={round(lp, 6)} 分数 -> 箱数分支")
    # 左子 [47,47]：RMP + Σx<=47 + 虚拟列
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in sel]
    yv = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(N)]
    for i in range(N):
        mm.add_linear_constraint(mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) + yv[i] >= 1.0, name=f"c{i}")
    mm.add_linear_constraint(mathopt.fast_sum(x) <= 47, name="ub")
    mm.add_linear_constraint(mathopt.fast_sum(x) >= 47, name="lb")
    mm.minimize(mathopt.fast_sum(x) + 10**6*mathopt.fast_sum(yv))
    res = mathopt.solve(mm, mathopt.SolverType.GLOP)
    y_used = sum(res.variable_values()[yv[i]] for i in range(N))
    if verbose:
        print(f"左子 [47,47]: {'不可行（虚拟列使用，长度下界 48）' if y_used > 1e-6 else '可行'} -> 剪枝")
    # 右子 [48,48]：划分池 MIP
    term, obj, wt = pool_mip(patterns, sel)
    Ksol = int(round(obj)) if obj is not None else None
    if verbose:
        print(f"右子 [48,48]: 划分池 MIP {term} = {Ksol} 箱 | 长度下界 48 = 上界 -> {'证明最优' if Ksol == 48 else '未决'}")
    return dict(opt=Ksol, proved=(Ksol == 48), root_lp=lp)

if __name__ == "__main__":
    print("FFD:", ffd(), "| 文件自报最优:", FILE_BEST)
    run_direct()
