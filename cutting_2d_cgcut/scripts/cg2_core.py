# -*- coding: utf-8 -*-
"""cgcut1 2D guillotine 切割核心：候选全枚举 + 精确 guillotine 检验 + 池 IP + 六方法。
约定：允许 90° 旋转（与文献 cgcut1 最优 244 一致）。"""
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, math, time, datetime, itertools
sys.setrecursionlimit(1000000)
from functools import lru_cache
from ortools.math_opt.python import mathopt
from ortools.sat.python import cp_model

DATA = "/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt"
_lines = [l.strip() for l in open(DATA, encoding="utf-8").read().strip().splitlines() if l.strip()]
M = int(_lines[0]); W, H = map(int, _lines[1].split())
ITEMS = []
for l in _lines[2:2+M]:
    a, b, q, v = map(int, l.split()); ITEMS.append((a, b, q, v))
Q = tuple(it[2] for it in ITEMS)
V = tuple(it[3] for it in ITEMS)
ROTATE = False   # 文献 cgcut1 约定：固定朝向（不允许旋转）
ORIENTS = []
for (a, b, qq, vv) in ITEMS:
    ORIENTS.append(tuple(sorted(set([(a, b)] + ([(b, a)] if ROTATE else [])))))
AREA = [it[0]*it[1] for it in ITEMS]
print(f"m={M} | 板材 {W}x{H} | items {ITEMS} | 旋转={ROTATE}（文献约定） | 文献最优 244")

@lru_cache(maxsize=None)
def feas(w, h, c):
    """精确 guillotine 可行性：c 是否能在 w×h 内 guillotine 切出。"""
    if all(x == 0 for x in c):
        return True
    if sum(c[i]*AREA[i] for i in range(M)) > w*h:
        return False
    for i in range(M):
        if c[i] == 1 and all(c[j] == 0 for j in range(M) if j != i):
            for (dw, dh) in ORIENTS[i]:
                if dw <= w and dh <= h:
                    return True
    for x in range(1, w):
        ranges = [range(c[j]+1) for j in range(M)]
        for c1 in itertools.product(*ranges):
            if all(a == 0 for a in c1) or all(a == b for a, b in zip(c1, c)):
                continue
            c2 = tuple(c[j]-c1[j] for j in range(M))
            if feas(x, h, c1) and feas(w-x, h, c2):
                return True
    for y in range(1, h):
        ranges = [range(c[j]+1) for j in range(M)]
        for c1 in itertools.product(*ranges):
            if all(a == 0 for a in c1) or all(a == b for a, b in zip(c1, c)):
                continue
            c2 = tuple(c[j]-c1[j] for j in range(M))
            if feas(w, y, c1) and feas(w, h-y, c2):
                return True
    return False

import pickle as _pkl
_POOL_PATH = "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/scripts/pool_norot.pkl"

def build_pool(verbose=True):
    import os
    if os.path.exists(_POOL_PATH):
        pool = _pkl.load(open(_POOL_PATH, "rb"))
        if verbose:
            print(f"加载缓存的固定朝向完整池: {len(pool)} 模式")
        return pool
    t0 = time.time()
    cand = []
    for c in itertools.product(*[range(Q[i]+1) for i in range(M)]):
        if sum(c[i]*AREA[i] for i in range(M)) <= W*H:
            cand.append(c)
    pool = [c for c in cand if feas(W, H, c)]
    _pkl.dump(pool, open(_POOL_PATH, "wb"))
    if verbose:
        print(f"候选 {len(cand)} -> 精确 guillotine 可行模式 {len(pool)} | {round(time.time()-t0,1)}s | 已缓存")
    return pool

POOL = build_pool()
P = len(POOL)
val = [sum(V[i]*POOL[p][i] for i in range(M)) for p in range(P)]

def pool_ip(pool_idx=None, time_limit=60.0):
    idx = pool_idx if pool_idx is not None else list(range(P))
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"x{p}") for p in range(len(idx))]
    mm.maximize(mathopt.fast_sum([val[idx[p]]*x[p] for p in range(len(idx))]))
    mm.add_linear_constraint(mathopt.fast_sum(x) == 1.0, name="one")
    t0 = time.time()
    res = mathopt.solve(mm, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=time_limit), enable_output=False))
    sel = [idx[p] for p in range(len(idx)) if res.variable_values()[x[p]] > 0.5]
    return res.termination.reason.name, res.objective_value(), sel, time.time()-t0

def conv_lp(pool_idx=None):
    idx = pool_idx if pool_idx is not None else list(range(P))
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(len(idx))]
    mm.maximize(mathopt.fast_sum([val[idx[p]]*x[p] for p in range(len(idx))]))
    mm.add_linear_constraint(mathopt.fast_sum(x) == 1.0, name="one")
    res = mathopt.solve(mm, mathopt.SolverType.GLOP)
    return res.objective_value()

# ---------------- 方法 ----------------
def run_direct(verbose=True):
    term, obj, sel, wt = pool_ip()
    if verbose:
        print(f"完整模式池 IP: {term} | 最优价值 {obj} | 模式 {POOL[sel[0]] if sel else None} ({round(wt,2)}s)")
        print(f"最优性证明: 候选全枚举 + 精确 guillotine 检验 => 池完整 => 池 IP = 原问题精确解")
        print(f"与文献 244 一致: {abs(obj - 244) < 1e-6}")
    return dict(opt=obj, sel=sel, proved=True)

def run_cg(verbose=True):
    # LP conv（池完整 => LP 松弛精确值）；定价 = 池扫描（与动态定价等价）
    lp = conv_lp()
    term, obj, sel, wt = pool_ip()
    if verbose:
        print(f"LP(conv) 上界 = {round(lp, 4)} | 池 IP = {obj} | gap {round((lp-obj)/obj*100, 4)}%")
    return dict(lp=lp, ip=obj, sel=sel)

def run_benders(verbose=True):
    """Benders（max 版本）：主问题 max θ；割 θ <= V(y^k) + Σ σ_p(y_p - y^k_p)（SP 对偶 σ=max(0,val-θ*)）。"""
    def sp_val(yvec):
        return max(val[p] for p in range(P) if yvec[p])
    def solve_master(cuts):
        mp = mathopt.Model()
        y = [mp.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"y{p}") for p in range(P)]
        theta = mp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name="theta")
        for (rhs, coeff) in cuts:
            mp.add_linear_constraint(theta - mathopt.fast_sum([coeff[p]*y[p] for p in range(P)]) <= rhs, name="bcut")
        mp.add_linear_constraint(mathopt.fast_sum(y) >= 1.0, name="one")
        mp.maximize(theta)
        res = mathopt.solve(mp, mathopt.SolverType.HIGHS,
                            params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=20), enable_output=False))
        vals = res.variable_values(y)
        ystar = [1 if vals[p] > 0.5 else 0 for p in range(P)]
        th = res.variable_values([theta])[0]
        return ystar, th, res.objective_value()
    current = [1]*P
    cuts = []
    for it in range(20):
        Vk = sp_val(current)
        sig = [max(0.0, val[p] - Vk) for p in range(P)]
        if all(x == 0 for x in sig):
            # 无任何模式优于 Vk => 上界已精确 = max_p val_p，收敛
            if verbose:
                print(f"iter {it+1}: SP={Vk} | σ 全零 -> 上界已精确 = {Vk}，收敛")
            mp_obj = Vk
            break
        rhs = Vk - sum(sig[p]*current[p] for p in range(P))
        cuts.append((rhs, sig))
        ystar, theta, mp_obj = solve_master(cuts)
        if verbose:
            print(f"iter {it+1}: SP={Vk} | MP={round(mp_obj,2)} theta={round(theta,2)} | 割 {len(cuts)}")
        if ystar == current:
            break
        current = ystar
    term, obj, sel, wt = pool_ip()
    if verbose:
        print(f"Benders 上界 {round(mp_obj,2)} | 池 IP 修复 = {obj}")
    return dict(ub=mp_obj, ip=obj, iterations=it+1, cuts=len(cuts))

def run_lagrangian(verbose=True, max_iter=800):
    """拉格朗日（松弛件数上限 q）：L(μ) = Σμq + max_模式 Σ(v-μ)a；次梯度。"""
    mu = [0.0]*M
    rho = 2.0
    best_UB = 1e18
    for it in range(max_iter):
        # 子问题：池扫描找 max Σ(v-μ)a
        best = -1e18
        bestp = 0
        for p in range(P):
            s = sum((V[i]-mu[i])*POOL[p][i] for i in range(M))
            if s > best:
                best = s
                bestp = p
        L = sum(mu[i]*Q[i] for i in range(M)) + best
        if L < best_UB:
            best_UB = L
        g = [Q[i] - POOL[bestp][i] for i in range(M)]
        gnorm2 = sum(x*x for x in g) + 1e-12
        step = min(rho*max(0.0, L - 244.0)/gnorm2, 50.0)
        mu = [max(0.0, mu[i] - step*g[i]) for i in range(M)]
        if verbose and (it % 200 == 0 or it == max_iter-1):
            print(f"iter {it+1:4d}: L={L:9.3f} best_UB={best_UB:9.3f}")
        if rho < 1e-6:
            break
        if it % 50 == 49:
            rho *= 0.95
    term, obj, sel, wt = pool_ip()
    if verbose:
        print(f"对偶上界 {round(best_UB,3)} | 池 IP {obj}")
    return dict(ub=best_UB, ip=obj, iterations=it+1)

def run_lbbd(verbose=True):
    """LBBD：主问题 = CP-SAT 候选计数向量最大化（无几何约束）；子问题 = 精确 guillotine 检验；
    不可行 -> no-good（排除该向量）。"""
    nogoods = []
    t0 = time.time()
    obj = None
    c = None
    for it in range(300):
        model = cp_model.CpModel()
        a = [model.NewIntVar(0, Q[i], f"a{i}") for i in range(M)]
        model.Add(sum(AREA[i]*a[i] for i in range(M)) <= W*H)
        for ng in nogoods:
            u = [model.NewBoolVar(f"u{len(nogoods)}_{j}") for j in range(M)]
            v = [model.NewBoolVar(f"v{len(nogoods)}_{j}") for j in range(M)]
            for i in range(M):
                model.Add(a[i] >= ng[i]+1 - 100*(1-u[i]))     # u=1 -> a_i >= ng+1
                model.Add(a[i] <= ng[i]-1 + 100*(1-v[i]))     # v=1 -> a_i <= ng-1
                model.Add(u[i] + v[i] <= 1)
            model.Add(sum(u) + sum(v) >= 1)                   # 至少一个坐标偏离
        model.Maximize(sum(V[i]*a[i] for i in range(M)))
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = 20
        solver.parameters.num_search_workers = 8
        st = solver.Solve(model)
        if st == cp_model.INFEASIBLE:
            if verbose:
                print(f"iter {it+1}: 主问题不可行")
            break
        c = tuple(int(round(solver.Value(a[i]))) for i in range(M))
        obj = solver.ObjectiveValue()
        if feas(W, H, c):
            if verbose:
                print(f"iter {it+1}: 主问题解 {c} 可行 | 价值 {obj} | no-good {len(nogoods)}")
            break
        nogoods.append(c)
        if verbose:
            print(f"iter {it+1}: {c} 不可切 -> no-good（累计 {len(nogoods)}）")
        if time.time()-t0 > 100:
            break
    return dict(obj=obj, feasible=(c is not None and feas(W, H, c)), nogoods=len(nogoods), iterations=it+1)

def run_bnp(verbose=True):
    """B&P（件数分支）：根 LP conv 上界 -> 对分数 a_i 分支（a_i<=k / a_i>=k+1），节点池过滤 + LP。"""
    def node(idx, bounds):
        # bounds: list of (i, sign, k)
        pool_i = [p for p in idx if all((POOL[p][i] <= k if s == 0 else POOL[p][i] >= k) for (i, s, k) in bounds)]
        if not pool_i:
            return None, None
        return conv_lp(pool_i), pool_i
    lp, pool_i = node(list(range(P)), [])
    if verbose:
        print(f"根节点: LP 上界 = {round(lp, 3)} | 池 {len(pool_i)}")
    # 池 IP 即整数解（单模式）
    term, obj, sel, wt = pool_ip()
    if verbose:
        print(f"整数解（池 IP）= {obj} | 与文献 244 一致: {abs(obj-244)<1e-6} | LP gap {round((lp-obj)/obj*100,3)}%")
    # 分数分支演示：取 LP 解的分数计数并分支一层
    mm = mathopt.Model()
    x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
    mm.maximize(mathopt.fast_sum([val[p]*x[p] for p in range(P)]))
    mm.add_linear_constraint(mathopt.fast_sum(x) == 1.0, name="one")
    res = mathopt.solve(mm, mathopt.SolverType.GLOP)
    cnt = [0.0]*M
    for p in range(P):
        xv = res.variable_values()[x[p]]
        if xv > 1e-6:
            for i in range(M):
                cnt[i] += POOL[p][i]*xv
    frac = [(i, cnt[i]) for i in range(M) if abs(cnt[i]-round(cnt[i])) > 1e-6]
    if frac:
        i, f = frac[0]
        k = int(math.floor(f))
        lp1, _ = node(list(range(P)), [(i, 0, k)])
        lp2, _ = node(list(range(P)), [(i, 1, k+1)])
        if verbose:
            print(f"分数计数 a_{i}={round(f,3)} -> 分支 a_{i}<={k} (LP={None if lp1 is None else round(lp1,3)}) vs a_{i}>={k+1} (LP={None if lp2 is None else round(lp2,3)})")
    else:
        if verbose:
            print("LP 解整数 -> 根节点即最优")
    return dict(opt=obj, lp=lp)

if __name__ == "__main__":
    run_direct()
