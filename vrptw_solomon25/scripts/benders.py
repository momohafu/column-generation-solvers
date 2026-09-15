# -*- coding: utf-8 -*-
"""VRPTW c101 Benders：二元选列主问题 + 连续覆盖 LP 子问题（scp_beasley 同款范式 + 车辆数约束）。
run_benders(mode) 返回结果 dict；直接运行时等价 python benders.py [full|cg]。
"""
import sys, os, time, datetime, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ortools.math_opt.python import mathopt
from ortools.sat.python import cp_model
from cg_cpsat import build_data, enumerate_pool, column_generation

M = 10**6

def solve_sp(yvec, P, n, K, paths, costs, masks):
    """子问题 SP(y)：候选列上的连续覆盖 LP（人工变量保可行）。返回 (alpha, lam, obj)。"""
    sp = mathopt.Model(name="sp")
    xv = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
    sv_ = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"s{i}") for i in range(1, n+1)]
    cov_cons = []
    for i in range(1, n+1):
        cov_cons.append(sp.add_linear_constraint(
            mathopt.fast_sum([xv[p] for p in range(P) if (masks[p] >> i) & 1]) + sv_[i-1] >= 1.0, name=f"cov{i}"))
    veh_con = sp.add_linear_constraint(mathopt.fast_sum(xv) <= K, name="vehicles")
    xub = [sp.add_linear_constraint(xv[p] <= float(yvec[p]), name=f"xub{p}") for p in range(P)]
    sp.minimize(mathopt.fast_sum([costs[p]*xv[p] for p in range(P)]) + M*mathopt.fast_sum(sv_))
    res = mathopt.solve(sp, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL, res.termination.reason
    dv = res.dual_values()
    pi = [max(0.0, dv[cov_cons[i-1]]) for i in range(1, n+1)]
    mu = dv[veh_con]
    alpha = sum(pi) + K*mu
    lam = [-dv[xub[p]] for p in range(P)]
    return alpha, lam, res.objective_value()

def solve_master(cuts, P, n, K, paths, costs, masks):
    """主问题：min theta + 最优性割（HIGHS MIP）。返回 (ystar, theta, mp_obj)。"""
    mp = mathopt.Model(name="master")
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

def run_benders(mode="cg", verbose=True):
    """完整 Benders 流程。mode: full=完整池, cg=列生成收敛池。返回结果 dict。"""
    n, xc, yc, dem, ready, due, svc, cap, depot_due, dist, d_scaled = build_data()
    K = 25
    t0 = time.time()
    full_paths, full_costs, full_masks, _ = enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due)
    t_enum = time.time() - t0
    if verbose:
        print("完整池:", len(full_paths), "列, 枚举耗时", round(t_enum, 2), "s")
    if mode == "full":
        pool = list(range(len(full_paths)))
        t_pool = 0.0
        cg_lp = None
    else:
        t0 = time.time()
        info = column_generation(verbose=False)
        pool = info["selected"]
        t_pool = time.time() - t0
        cg_lp = info["lp_obj"]
        if verbose:
            print("CG 列池:", len(pool), "| CG LP 下界:", round(cg_lp, 6), "| 列池构建耗时", round(t_pool, 2), "s")
    paths = [full_paths[p] for p in pool]
    costs = [full_costs[p] for p in pool]
    masks = [full_masks[p] for p in pool]
    P = len(pool)
    current = [1]*P
    cuts = []
    max_iter = 20
    wall0 = time.time()
    for it in range(max_iter):
        t1 = time.time()
        alpha, lam, sp_obj = solve_sp(current, P, n, K, paths, costs, masks)
        sp_t = time.time()-t1
        cuts.append((alpha, lam))
        t1 = time.time()
        ystar, theta, mp_obj = solve_master(cuts, P, n, K, paths, costs, masks)
        mp_t = time.time()-t1
        if verbose:
            print(f"iter {it+1:2d}: SP_obj={sp_obj:12.4f} ({sp_t:4.1f}s) | MP_obj={mp_obj:12.4f} theta={theta:12.4f} ({mp_t:4.1f}s) | y 选中 {sum(ystar)} 列 | 割 {len(cuts)}")
        if ystar == current:
            if verbose:
                print("master 解不变 -> 收敛")
            break
        current = ystar
        if time.time()-wall0 > 110:
            if verbose:
                print("时间上限")
            break
    wall = time.time()-wall0
    # 整数修复：候选池上 CP-SAT 集合覆盖整数模型
    t1 = time.time()
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
    t_repair = time.time()-t1
    routes = [paths[p] for p in range(P) if solver.Value(yv[p])]
    exact = None
    if routes:
        exact = 0.0
        for r in routes:
            seq = [0] + list(r) + [0]
            exact += sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1))
    if verbose:
        print("IP 修复:", solver.StatusName(st), "| 目标(分):", solver.ObjectiveValue(), "| 时间:", round(t_repair, 2), "s")
        for r in routes:
            seq = [0] + list(r) + [0]
            print("  路线", r, "距离", round(sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1)), 6))
        print("车辆数", len(routes), "| 精确总距离", None if exact is None else round(exact, 6),
              "| Benders 下界", round(mp_obj, 6),
              "| gap", round((exact - mp_obj)/mp_obj*100, 6), "%")
    return dict(mode=mode, pool=len(pool), cg_lp=cg_lp, iterations=it+1, cuts=len(cuts),
                lb=mp_obj, wall=wall, t_enum=t_enum, t_pool=t_pool, t_repair=t_repair,
                repair_status=solver.StatusName(st), obj_cent=solver.ObjectiveValue(),
                exact=exact, routes=routes)

if __name__ == "__main__":
    run_benders(sys.argv[1] if len(sys.argv) > 1 else "cg")
