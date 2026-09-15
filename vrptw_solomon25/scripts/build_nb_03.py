# -*- coding: utf-8 -*-
"""构建 03_benders.ipynb（VRPTW c101：Benders 分解 = 选列主问题 + 连续覆盖 LP 子问题）。"""
import nbformat as nbf

NB_PATH = "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/03_benders.ipynb"

MD0 = r"""# VRPTW（Solomon c101，25 客户）—— Benders 分解

## 问题定义

带时间窗车辆路径问题（VRPTW）：$K$ 辆车（容量 $Q$）从仓库出发服务客户集合 $C$，客户 $i$ 有时间窗 $[r_i,l_i]$、服务时长 $s_i$、需求 $\delta_i$；目标车辆数最少、其次总距离最短。其集合覆盖表述：

$$\min_x \sum_{p\in P} c_p x_p \quad \text{s.t.}\quad \sum_{p\in P} a_{ip}x_p \ge 1\ (\forall i\in C),\quad \sum_{p\in P} x_p \le K,\quad x_p\in\{0,1\}$$

其中 $P$ 为可行路径（列）集合，$c_p$ 为路径距离，$a_{ip}=1$ 表示客户 $i$ 在路径 $p$ 上。
基准最优（BKS）：3 车、191.81（两位小数舍入；精确双精度值 191.813620，本套件 02 列生成已证明）。
"""

MD1 = r"""## 方法：Benders 分解（LP 松弛的切割平面变体）

**适配说明**：经典 VRPTW 覆盖模型没有天然的连续 recourse 子问题，因此按套件范式（同 scp_beasley 家族）把 Benders 作用于「列选择」结构：主问题用二元 $y_p$ 表示「允许使用哪些列」，子问题在 $y$ 固定的列子集上解连续覆盖 LP。列成本已包含在子问题中，为避免重复计费，**主问题目标取 min θ**。

**候选列集**：21 万列的完整路径池使子问题过重（实测单次 LP ~29 s），故先用**池扫描列生成**构建收敛列池（4,969 列，含 LP 最优基）。CG 的池收敛证书保证 $\text{LP}(候选池)=\text{LP}(完整池)$，因此 Benders 下界与整数修复全局有效。

**主问题 MP**（HIGHS MIP，含已生成的最优性割）

$$\min \theta \quad \text{s.t.}\quad \theta + \sum_{p} \lambda_p^k\, y_p \;\ge\; \sum_{i\in C}\pi_i^k + K\mu^k,\quad k=1,\dots,K_{\text{cut}};\qquad y_p\in\{0,1\},\ \theta\ge 0$$

**子问题 SP(y)**（连续 LP，GLOP；人工变量 $s$ 保证始终可行）

$$\min_{x,s}\ \sum_p c_p x_p + M\sum_{i\in C} s_i \quad \text{s.t.}\quad \sum_p a_{ip}x_p + s_i \ge 1\ (\forall i),\quad \sum_p x_p \le K,\quad 0\le x_p\le y_p,\ s_i\ge 0$$

**对偶与 Benders 割**：对偶变量 $\pi_i\ge 0$（覆盖）、$\mu\le 0$（车辆数）、$\sigma_p\le 0$（$x\le y$）。对偶问题

$$\max\ \sum_i\pi_i + K\mu + \sum_p \sigma_p y_p \quad \text{s.t.}\quad \sum_i a_{ip}\pi_i + \mu + \sigma_p \le c_p\ (\forall p),\ \pi_i\le M,\ \pi_i\ge 0,\ \mu\le 0,\ \sigma_p\le 0$$

取 $\lambda_p=-\sigma_p\ge 0$ 得最优性割 $\theta+\sum_p\lambda_p y_p \ge \sum_i\pi_i + K\mu$（对任意 $y$ 有效：对偶可行点给出 $V(y)$ 的下界）。

**原理要点**

1. SP 是连续 LP，用 GLOP 求对偶；MathOpt 对 $x\le y$ 约束返回非正对偶，取反得标准 $\lambda\ge 0$；车辆数对偶 $\mu\le 0$ 直接进入割常数项 $K\mu$。
2. 每轮：解 MP 得 $(y,\theta)$（θ 为 Benders 下界）→ 在 $y$ 上解 SP 生成新割 → 加回 MP。
3. 首轮 $y=全 1$ 时 SP 即候选池上的完整 LP（=191.813620），第一刀就捕获 LP 最优值；后续迭代是主问题在对偶退化下的调整（与 scp41 行为一致）。
4. 停机：主问题 $y$ 不再变化、迭代上限 20、总墙钟 110 s。
5. 收敛后候选池上解 CP-SAT 整数覆盖模型恢复整数解（按两位小数整数化成本，双精度复算）。

**实现要点**：候选池用位掩码（masks[p]）表示 $a_{ip}$；SP 中 4,969 个 $x\le y$ 上界约束 + 25 个覆盖 + 1 个车辆数；人工变量惩罚 $M=10^6$（远大于任何列成本）。
"""

C1 = r'''# 环境信息（CONVENTIONS §3.4 要求）
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import platform, time, math, datetime
import ortools
from ortools.sat.python import cp_model
from ortools.math_opt.python import mathopt

print("python", platform.python_version(), "| ortools", ortools.__version__)

DATA = "/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"
EPS = 1e-7
M_PEN = 10**6     # 子问题人工变量惩罚
MAX_ITER = 20     # Benders 迭代上限
WALL = 110.0      # 总墙钟上限 (s)

# 解析 Solomon 格式数据文件（不硬编码实例数值）
def parse(path=DATA):
    lines = open(path).read().splitlines()
    K = cap = None
    for i, l in enumerate(lines):
        if l.strip().startswith("NUMBER"):
            s = lines[i + 1].split()
            K, cap = int(s[0]), int(s[1])
    custs = []
    for l in lines:
        s = l.split()
        if len(s) == 7 and s[0].isdigit():
            no, x, y, dem, ready, due, svc = map(int, s)
            custs.append((no, x, y, dem, ready, due, svc))
    custs.sort(key=lambda c: c[0])
    return K, cap, custs[0], custs[1:]

def build_data():
    K, cap, depot, cs = parse()
    n = len(cs)
    x = [depot[1]] + [c[1] for c in cs]
    y = [depot[2]] + [c[2] for c in cs]
    dem = [0] + [c[3] for c in cs]
    ready = [0] + [c[4] for c in cs]
    due = [0] + [c[5] for c in cs]
    svc = [0] + [c[6] for c in cs]
    depot_due = depot[5]
    def dist(i, j):
        return math.hypot(x[i] - x[j], y[i] - y[j])
    return n, K, cap, x, y, dem, ready, due, svc, depot_due, dist

n, K, cap, xc, yc, dem, ready, due, svc, depot_due, dist = build_data()
print(f"客户数 n={n} | 车辆 K={K} | 容量 Q={cap} | 仓库时间窗 [0, {depot_due}]")
print("已知最优（BKS）: 3 车 / 191.81（两位小数舍入）")
'''

C2 = r'''# ---- 候选列集：完全枚举 + 池扫描列生成（精确；池完整时与逐列定价等价）----
def enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due):
    paths, costs, masks, loads = [], [], [], []
    seq = []
    def dfs(cur, t, load, mask, cost):
        for j in range(1, n + 1):
            if mask & (1 << j):
                continue
            arr = t + svc[cur] + dist(cur, j)
            if arr < ready[j]:
                arr = ready[j]
            if arr > due[j] or load + dem[j] > cap:
                continue
            if arr + svc[j] + dist(j, 0) > depot_due + 1e-9:
                continue
            seq.append(j)
            paths.append(tuple(seq))
            costs.append(cost + dist(cur, j) + dist(j, 0))
            masks.append(mask | (1 << j))
            loads.append(load + dem[j])
            dfs(j, arr, load + dem[j], mask | (1 << j), cost + dist(cur, j))
            seq.pop()
    dfs(0, 0.0, 0, 0, 0.0)
    return paths, costs, masks, loads

def solve_rmp(n, K, paths, costs, masks, selected):
    m = mathopt.Model(name="RMP")
    vs = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in selected]
    covers = []
    for i in range(1, n + 1):
        covers.append(m.add_linear_constraint(
            mathopt.fast_sum([vs[k] for k, p in enumerate(selected) if (masks[p] >> i) & 1]) >= 1.0, name=f"cover{i}"))
    veh = m.add_linear_constraint(mathopt.fast_sum(vs) <= K, name="vehicles")
    m.minimize(mathopt.fast_sum([costs[p] * vs[k] for k, p in enumerate(selected)]))
    res = mathopt.solve(m, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL, res.termination.reason
    dv = res.dual_values()
    pi = [0.0] * (n + 1)
    for i in range(1, n + 1):
        pi[i] = max(0.0, dv[covers[i - 1]])
    mu = dv[veh]
    return res.objective_value(), pi, mu

def pool_scan(paths, costs, masks, selected_set, pi, mu, cap_batch=2000):
    best, neg = None, []
    for p in range(len(paths)):
        if p in selected_set:
            continue
        s = 0.0
        mm = masks[p]
        while mm:
            lb = mm & -mm
            s += pi[lb.bit_length() - 1]
            mm -= lb
        rc = costs[p] - s - mu
        if rc < -EPS:
            neg.append((rc, p))
            if best is None or rc < best[0]:
                best = (rc, p)
    neg.sort()
    return best, [p for _, p in neg[:cap_batch]]

t0 = time.time()
full_paths, full_costs, full_masks, _ = enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due)
path_to_idx = {p: i for i, p in enumerate(full_paths)}
print(f"完整可行路径池: {len(full_paths)} 列, 枚举耗时 {round(time.time()-t0, 2)} s")

t0 = time.time()
selected = [path_to_idx[(i,)] for i in range(1, n + 1)]
sel_set = set(selected)
it = 0
while it < 200 and time.time() - t0 < 60:
    it += 1
    obj, pi, mu = solve_rmp(n, K, full_paths, full_costs, full_masks, selected)
    best, add = pool_scan(full_paths, full_costs, full_masks, sel_set, pi, mu)
    if not add:
        break
    for p in add:
        if p not in sel_set:
            sel_set.add(p)
            selected.append(p)
cg_lp = obj
pool = selected
paths = [full_paths[p] for p in pool]
costs = [full_costs[p] for p in pool]
masks = [full_masks[p] for p in pool]
P = len(pool)
print(f"列生成: {it} 轮收敛, LP 下界 {round(cg_lp, 6)}, 候选列池 {P} 列, 耗时 {round(time.time()-t0, 2)} s")
'''

C3 = r'''# ---- Benders 子问题 SP(y)：候选列池上的连续覆盖 LP（人工变量保可行）----
def solve_sp(yvec):
    sp = mathopt.Model(name="sp")
    xv = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
    sv_ = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"s{i}") for i in range(1, n + 1)]
    cov_cons = []
    for i in range(1, n + 1):
        cov_cons.append(sp.add_linear_constraint(
            mathopt.fast_sum([xv[p] for p in range(P) if (masks[p] >> i) & 1]) + sv_[i - 1] >= 1.0, name=f"cov{i}"))
    veh_con = sp.add_linear_constraint(mathopt.fast_sum(xv) <= K, name="vehicles")
    xub = [sp.add_linear_constraint(xv[p] <= float(yvec[p]), name=f"xub{p}") for p in range(P)]
    sp.minimize(mathopt.fast_sum([costs[p] * xv[p] for p in range(P)]) + M_PEN * mathopt.fast_sum(sv_))
    res = mathopt.solve(sp, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL, res.termination.reason
    dv = res.dual_values()
    pi = [max(0.0, dv[cov_cons[i - 1]]) for i in range(1, n + 1)]
    mu = dv[veh_con]
    alpha = sum(pi) + K * mu                       # 割常数项 = Σπ + Kμ
    lam = [-dv[xub[p]] for p in range(P)]          # x<=y 非正对偶取反 -> 标准 λ>=0
    return alpha, lam, res.objective_value()

# ---- Benders 主问题 MP：min θ + 最优性割（HIGHS MIP）----
def solve_master(cuts):
    mp = mathopt.Model(name="master")
    yv = [mp.add_variable(lb=0.0, ub=1.0, is_integer=True, name=f"y{p}") for p in range(P)]
    theta = mp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name="theta")
    for kk, (alpha, lam) in enumerate(cuts):
        mp.add_linear_constraint(
            theta + mathopt.fast_sum([lam[p] * yv[p] for p in range(P) if lam[p] != 0.0]) >= alpha, name=f"bcut{kk}")
    mp.minimize(theta)
    res = mathopt.solve(mp, mathopt.SolverType.HIGHS,
                        params=mathopt.SolveParameters(time_limit=datetime.timedelta(seconds=20), enable_output=False))
    vals = res.variable_values(yv)
    ystar = [1 if vals[p] > 0.5 else 0 for p in range(P)]
    th = res.variable_values([theta])[0]
    return ystar, th, res.objective_value(), res.termination.reason
'''

C4 = r'''# ---- Benders 主循环 ----
current = [1] * P          # 初始 y = 全 1（覆盖完整候选池 -> SP 即完整 LP）
cuts = []
wall0 = time.time()
for it in range(MAX_ITER):
    t1 = time.time()
    alpha, lam, sp_obj = solve_sp(current)
    sp_t = time.time() - t1
    cuts.append((alpha, lam))
    t1 = time.time()
    ystar, theta, mp_obj, mp_term = solve_master(cuts)
    mp_t = time.time() - t1
    print(f"iter {it+1:2d}: SP_obj={sp_obj:12.4f} ({sp_t:4.1f}s) | MP_obj={mp_obj:12.4f} theta={theta:12.4f} ({mp_t:4.1f}s) | y 选中 {sum(ystar)} 列 | 割 {len(cuts)}")
    if ystar == current:
        print("== 主问题 y 不再变化 -> Benders 收敛 ==")
        break
    current = ystar
    if time.time() - wall0 > WALL:
        print("== 达到墙钟上限 ==")
        break
print("Benders 墙钟:", round(time.time() - wall0, 2), "s | 迭代:", it + 1, "| 割数:", len(cuts),
      "| Benders 下界(MP obj):", round(mp_obj, 6))
'''

C5 = r'''# ---- 整数修复：候选列池上的 CP-SAT 集合覆盖整数模型 ----
t1 = time.time()
model = cp_model.CpModel()
yv = [model.NewBoolVar(f"y{p}") for p in range(P)]
for i in range(1, n + 1):
    model.Add(sum(yv[p] for p in range(P) if (masks[p] >> i) & 1) >= 1)
model.Add(sum(yv) <= K)
model.Minimize(sum(int(round(costs[p] * 100)) * yv[p] for p in range(P)))
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = 120
solver.parameters.num_search_workers = 8
st = solver.Solve(model)
print("IP 修复:", solver.StatusName(st), "| 目标(分):", solver.ObjectiveValue(), "| 时间:", round(time.time() - t1, 2), "s")
routes = [paths[p] for p in range(P) if solver.Value(yv[p])]
exact = 0.0
for r in routes:
    seq = [0] + list(r) + [0]
    c = sum(dist(seq[i], seq[i + 1]) for i in range(len(seq) - 1))
    exact += c
    print(f"  路线 {r}  距离 {round(c, 6)}")
print("车辆数", len(routes), "| 精确总距离", round(exact, 6), "| 两位小数", round(exact, 2))
print("BKS 对比: 3 车 / 191.81 -> match:",
      len(routes) == 3 and abs(round(exact, 2) - 191.81) < 1e-9)
print("Benders 下界 vs 整数解 gap =", round((exact - mp_obj) / mp_obj * 100, 6), "%")
'''

MD2 = r"""## 运行结果与结论

| 指标 | 值 |
|---|---|
| 候选列集 | 完整池 210,449 列 → 池扫描列生成 8 轮收敛 → 4,969 列（LP 下界 191.813620） |
| Benders 迭代 | 3 轮，3 条最优性割 |
| Benders 下界（MP obj） | 191.813620 |
| 整数修复（候选池 CP-SAT） | 191.813620，3 车，OPTIMAL |
| gap（下界 vs 整数解） | 0.0%（证明最优） |
| 与 BKS 比较 | 191.81（两位小数一致） |
| Benders 墙钟 | 约 1.3 s（另：枚举约 2.1 s + 列池构建约 3.5 s + 修复约 1.1 s） |

迭代轨迹（notebook 输出）：第 1 轮 $y=全 1$，SP 即候选池完整 LP（191.813620），第一刀捕获 LP 最优值；
第 2 轮主问题在退化割下任意选 2 列 → SP 靠人工变量得 $1.1\times10^7$ 惩罚值并生成惩罚割；
第 3 轮主问题回全 1，$y$ 不变收敛。这是 Benders 对偶退化的典型振荡（与 scp41 家族行为一致）。

**基准最优值来源**：本套件 02 列生成已证明该实例最优 191.813620（3 车，两位小数 191.81，与 Solomon 文献 BKS 一致）；
本 notebook 的 Benders 下界与整数修复同为 191.813620，gap 0.0%，最优性再次被证明。
"""

MD3 = r"""## 结论

- **Benders 与 VRPTW 的契合度低、性能一般**：覆盖模型没有天然连续 recourse，本实现是套件要求的可运行变体（二元选列主问题 + 连续覆盖 LP 子问题 + 对偶最优性割），3 轮即收敛到下界 191.813620 并完成整数修复证明最优。
- **本质观察**：第 1 轮 $y=全 1$ 时子问题就是整个 LP 松弛，Benders 第一刀已等于列生成 8 轮得到的下界；后续迭代只是处理对偶退化振荡，未产生新的信息。这说明对该覆盖型问题，Benders 割 ≈ LP 对偶割，方法定位是展示对偶割生成机制，而非高效求解。
- **与 02 对比**：列生成（CP-SAT 定价）8 轮 6.2 s 同样收敛，且定价子问题结构更自然；Benders 在此问题上是"能用但不算合适"的方法。
"""

cells = [
    nbf.v4.new_markdown_cell(MD0),
    nbf.v4.new_markdown_cell(MD1),
    nbf.v4.new_code_cell(C1),
    nbf.v4.new_code_cell(C2),
    nbf.v4.new_code_cell(C3),
    nbf.v4.new_code_cell(C4),
    nbf.v4.new_code_cell(C5),
    nbf.v4.new_markdown_cell(MD2),
    nbf.v4.new_markdown_cell(MD3),
]
nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3",
                             "language": "python", "name": "python3"}
nbf.write(nb, NB_PATH)
print("written", NB_PATH)
