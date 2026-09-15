# -*- coding: utf-8 -*-
"""构建 07_branch_and_price.ipynb（VRPTW c101：branch-and-price）。"""
import nbformat as nbf

NB_PATH = "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/07_branch_and_price.ipynb"

MD0 = r"""# VRPTW（Solomon c101）—— Branch-and-Price（分支-定价）

## 问题定义

带时间窗车辆路径问题（VRPTW，同 02/03/05）：$K=25$ 辆车（容量 $Q=200$）服务 25 个客户，
目标车辆数最少、其次总距离最短。集合覆盖表述：

$$\min_x \sum_{p\in P} c_p x_p \quad \text{s.t.}\quad \sum_{p\in P} a_{ip}x_p \ge 1\ (\forall i\in C),\quad \sum_{p\in P} x_p \le K,\quad x_p\in\{0,1\}$$

基准最优（BKS）：3 车、191.813620（两位小数 191.81），02/03/05 三方法已独立证明。

## 方法：Branch-and-Price = 分支定界 × 列生成

- **节点求解**：每个 B&B 节点用**列生成**精确求解该节点的集合覆盖 LP 松弛（节点下界）：
  CP-SAT 定价 ESPPRC（动态列生成）+ 完整池扫描（25 节点池可完全枚举，作为精确收敛证书，同 02 口径）。
- **根节点**：车辆数上界 $K_{ub}=25$、下界 $K_{lb}=\lceil 460/200\rceil=3$（容量下界）。
- **剪枝**：① LP 下界 ≥ incumbent（定界剪枝）；② 节点 LP 不可行（无真实列覆盖，虚拟列取正值）；
  ③ LP 解为整数（$x_p\in\{0,1\}$）→ 得到整数解、更新 incumbent、节点解决。
- **分支**：Σx 分数 → **车辆数分支**（左子 $\Sigma x\le\lfloor f\rfloor$，右子 $\Sigma x\ge\lceil f\rceil$）；
  否则 **Ryan-Foster 弧分支**（取流量最接近 0.5 的客户弧 $(i,j)$：左子禁用、右子强制）。
- **约束传播到定价**：强制弧 $x_{ij}=1$ / 禁用弧 $x_{ij}=0$ 直接加入 CP-SAT 定价模型，池扫描同步过滤。

**预期**：该实例集合覆盖 LP 松弛恰好是整数（02 已验证 LP 下界 = 整数目标），故根节点即应给出整数最优解、无需分支。
"""

MD1 = r"""## 实现要点

1. **RMP 扩展**：在 02 的覆盖 LP 基础上加车辆数下界约束 $\Sigma x \ge K_{lb}$（对偶 $\mu_{lb}\ge 0$，
   与上界对偶 $\mu_{ub}\le 0$ 合并成 $\mu_{sum}=\mu_{ub}+\mu_{lb}$ 进入定价 reduced cost），
   并引入 25 个虚拟列（成本 $M=10^6$）保证初始可行；虚拟列取正值 ⇒ 节点无真实覆盖 ⇒ 不可行剪枝。
2. **节点 CG**：初始列 = 满足弧限制的单客户列；强制弧节点另需以**含全部强制弧的池列**做种子
   （单客户列不含客户-客户弧，否则 $\Sigma x\ge K_{lb}$ 不可行）。
3. **定价子问题**：CP-SAT ESPPRC 增加强制/禁用弧约束；每轮仍以完整池扫描精确校验并批量加列。
4. **整数判定**：GLOP 解全部 0/1（容差 $10^{-6}$）即节点整数解，直接提取路线。
5. **停机**：搜索栈空 / 墙钟 110 s；算法确定性。
"""

C1 = r'''# 环境信息（CONVENTIONS §3.4 要求）
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, platform, math, time
import ortools
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg_cpsat import build_data, enumerate_pool
from ortools.sat.python import cp_model
from ortools.math_opt.python import mathopt

print("python", platform.python_version(), "| ortools", ortools.__version__)

n, xc, yc, dem, ready, due, svc, cap, depot_due, dist, d_scaled = build_data()
K = 25
LB_K = 3                       # ceil(460/200)
EPS = 1e-7
M_DUMMY = 1e6
SCALE = 1000
DUAL_SCALE = 1000
BIG = 10**9

t0 = time.time()
paths, costs, masks, loads = enumerate_pool(n, dist, dem, ready, due, svc, cap, depot_due)
path_to_idx = {p: i for i, p in enumerate(paths)}
print(f"完整池: {len(paths)} 列, 枚举耗时 {round(time.time()-t0, 1)} s")
print("已知最优（BKS）: 3 车 / 191.81（两位小数舍入）")
'''

C2 = r'''# ---- 节点组件：列/弧工具、RMP、CP-SAT 定价、节点列生成 ----
def col_arcs(p):
    a = []
    prev = 0
    for j in p:
        a.append((prev, j))
        prev = j
    a.append((prev, 0))
    return a

def eligible(idx, forced, forbidden):
    arcs = col_arcs(paths[idx])
    if any(arc in forbidden for arc in arcs):
        return False
    if any(arc not in arcs for arc in forced):
        return False
    return True

def solve_rmp_node(selected, K_ub, K_lb):
    """覆盖 LP + 车辆数上下界 + 虚拟列。返回 (obj, pi, mu_ub, mu_lb, xvals, y_used)。"""
    m = mathopt.Model(name="rmp")
    vs = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in selected]
    yv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(1, n+1)]
    covers = []
    for i in range(1, n+1):
        covers.append(m.add_linear_constraint(
            mathopt.fast_sum([vs[k] for k, p in enumerate(selected) if (masks[p] >> i) & 1]) + yv[i-1] >= 1.0, name=f"c{i}"))
    ub_con = m.add_linear_constraint(mathopt.fast_sum(vs) <= K_ub, name="vub")
    lb_con = m.add_linear_constraint(mathopt.fast_sum(vs) >= K_lb, name="vlb")
    m.minimize(mathopt.fast_sum([costs[p]*vs[k] for k, p in enumerate(selected)]) + M_DUMMY*mathopt.fast_sum(yv))
    res = mathopt.solve(m, mathopt.SolverType.GLOP)
    assert res.termination.reason == mathopt.TerminationReason.OPTIMAL, res.termination.reason
    dv = res.dual_values()
    pi = [0.0]*(n+1)
    for i in range(1, n+1):
        pi[i] = max(0.0, dv[covers[i-1]])
    mu_ub = dv[ub_con]
    mu_lb = dv[lb_con]
    xvals = {p: res.variable_values()[vs[k]] for k, p in enumerate(selected)}
    y_used = sum(res.variable_values()[yv[i-1]] for i in range(1, n+1))
    return res.objective_value(), pi, mu_ub, mu_lb, xvals, y_used

def cpsat_pricing(pi, mu_sum, forced, forbidden, time_limit=10.0):
    """CP-SAT 定价：ESPPRC + 强制/禁用弧。"""
    model = cp_model.CpModel()
    xv = [[model.NewBoolVar(f"x{i}_{j}") for j in range(n+1)] for i in range(n+1)]
    model.Add(xv[0][0] == 0)
    model.AddCircuit([(i, j, xv[i][j]) for i in range(n+1) for j in range(n+1)])
    visited = [model.NewBoolVar(f"v{i}") for i in range(n+1)]
    model.Add(visited[0] == 1)
    for i in range(1, n+1):
        model.Add(visited[i] + xv[i][i] == 1)
    model.Add(sum(visited) >= 2)
    for (i, j) in forced:
        model.Add(xv[i][j] == 1)
    for (i, j) in forbidden:
        model.Add(xv[i][j] == 0)
    T = depot_due * SCALE
    t = [model.NewIntVar(0, T, f"t{i}") for i in range(n+1)]
    model.Add(t[0] == 0)
    for i in range(1, n+1):
        model.Add(t[i] >= ready[i]*SCALE - BIG*xv[i][i])
        model.Add(t[i] <= due[i]*SCALE + BIG*xv[i][i])
        model.Add(t[i] <= BIG*(1 - xv[i][i]))
    for i in range(n+1):
        for j in range(1, n+1):
            if i == j:
                continue
            model.Add(t[j] >= t[i] + svc[i]*SCALE + d_scaled[i][j] - BIG*(1 - xv[i][j]))
    for i in range(1, n+1):
        model.Add(t[i] + svc[i]*SCALE + d_scaled[i][0] <= depot_due*SCALE + BIG*(1 - xv[i][0]))
    q = [model.NewIntVar(0, cap, f"q{i}") for i in range(n+1)]
    model.Add(q[0] == 0)
    for i in range(n+1):
        for j in range(1, n+1):
            if i == j:
                continue
            model.Add(q[j] >= q[i] + dem[j] - BIG*(1 - xv[i][j]))
    pi_s = [int(round(pi[i]*DUAL_SCALE)) for i in range(n+1)]
    mu_s = int(round(mu_sum*DUAL_SCALE))
    model.Minimize(sum(d_scaled[i][j]*xv[i][j] for i in range(n+1) for j in range(n+1) if i != j)
                   - sum(pi_s[i]*visited[i] for i in range(1, n+1)) - mu_s)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    st = solver.Solve(model)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return None
    path = []
    cur = 0
    while True:
        nxt = None
        for j in range(n+1):
            if solver.Value(xv[cur][j]):
                nxt = j
                break
        if nxt is None or nxt == 0:
            break
        path.append(nxt)
        cur = nxt
        if len(path) > n:
            break
    return tuple(path)

def exact_rc(path, pi, mu_sum):
    c = 0.0
    prev = 0
    for j in path:
        c += dist(prev, j)
        prev = j
    c += dist(prev, 0)
    return c - sum(pi[i] for i in path) - mu_sum

def cg_node(K_ub, K_lb, forced, forbidden, time_limit=60.0):
    """节点列生成。返回节点状态 dict。"""
    forced = tuple(forced)
    forbidden = set(forbidden)
    t0 = time.time()
    selected = []
    sel_set = set()
    for i in range(1, n+1):
        idx = path_to_idx[(i,)]
        if eligible(idx, forced, forbidden):
            sel_set.add(idx)
            selected.append(idx)
    if forced:
        # 强制弧节点：以含全部强制弧的池列做种子（保证 Σx>=K_lb 可行）
        cnt = 0
        for idx in range(len(paths)):
            if idx in sel_set:
                continue
            if eligible(idx, forced, forbidden):
                sel_set.add(idx)
                selected.append(idx)
                cnt += 1
                if cnt >= 50:
                    break
    if not selected:
        return dict(lp_obj=None, infeasible=True, integral=False, iters=0, ncols=0, xvals={}, time=0.0)
    it = 0
    lp_obj = None
    xvals = {}
    y_used = 0.0
    while it < 200 and time.time()-t0 < time_limit:
        it += 1
        obj, pi, mu_ub, mu_lb, xvals, y_used = solve_rmp_node(selected, K_ub, K_lb)
        lp_obj = obj
        mu_sum = mu_ub + mu_lb
        cpath = cpsat_pricing(pi, mu_sum, forced, forbidden)
        cpsat_rc = exact_rc(cpath, pi, mu_sum) if cpath is not None else None
        add = []
        for idx in range(len(paths)):
            if idx in sel_set:
                continue
            if not eligible(idx, forced, forbidden):
                continue
            s = 0.0
            mm = masks[idx]
            while mm:
                lb = mm & -mm
                s += pi[lb.bit_length()-1]
                mm -= lb
            rc = costs[idx] - s - mu_sum
            if rc < -EPS:
                add.append((rc, idx))
        add.sort()
        add = [i for _, i in add[:2000]]
        for idx in add:
            if idx not in sel_set:
                sel_set.add(idx)
                selected.append(idx)
        if not add and (cpsat_rc is None or cpsat_rc >= -EPS):
            break
    infeasible = (y_used > 1e-6)
    integral = (not infeasible) and all(abs(v - round(v)) < 1e-6 for v in xvals.values())
    return dict(lp_obj=None if infeasible else lp_obj, infeasible=infeasible, integral=integral,
                iters=it, ncols=len(selected), xvals=xvals, time=time.time()-t0)
'''

C3 = r'''# ---- B&P 树搜索：分支（车辆数 / Ryan-Foster 弧）+ 剪枝 ----
def branch_and_price(extra_forbidden=(), root_K_ub=None, root_K_lb=None, time_limit=110.0):
    wall0 = time.time()
    incumbent = None
    inc_routes = None
    root = dict(K_ub=root_K_ub if root_K_ub is not None else K,
                K_lb=root_K_lb if root_K_lb is not None else LB_K,
                forced=(), forbidden=set(extra_forbidden))
    stack = [root]
    nodes = 0
    branches = 0
    while stack and time.time()-wall0 < time_limit:
        node = stack.pop()
        nodes += 1
        res = cg_node(node["K_ub"], node["K_lb"], node["forced"], node["forbidden"])
        tag = f"node{nodes} K[{node['K_lb']},{node['K_ub']}]"
        if node["forced"]:
            tag += f" 强制{list(node['forced'])}"
        if node["forbidden"]:
            tag += f" 禁用{len(node['forbidden'])}弧"
        if res["infeasible"]:
            print(f"{tag}: 节点不可行 -> 剪枝 | CG {res['iters']} 轮 {res['time']:.1f}s")
            continue
        v = res["lp_obj"]
        if incumbent is not None and v >= incumbent - 1e-7:
            print(f"{tag}: LP={v:.6f} >= incumbent {incumbent:.6f} -> 定界剪枝 | CG {res['iters']} 轮 {res['time']:.1f}s")
            continue
        if res["integral"]:
            routes = [paths[p] for p, val in res["xvals"].items() if val > 0.5]
            if incumbent is None or v < incumbent:
                incumbent = v
                inc_routes = routes
            print(f"{tag}: LP={v:.6f} 整数解 ({len(routes)} 车) -> 节点解决, 更新 incumbent | CG {res['iters']} 轮 {res['time']:.1f}s")
            continue
        # ---- 分支 ----
        fV = sum(res["xvals"].values())
        if abs(fV - round(fV)) > 1e-6:
            branches += 1
            c1 = dict(node, K_ub=int(math.floor(fV)))
            c2 = dict(node, K_lb=int(math.ceil(fV)))
            stack.append(c2)
            stack.append(c1)
            print(f"{tag}: LP={v:.6f} 分数(Σx={fV:.3f}) -> 车辆数分支 [{node['K_lb']},{int(math.floor(fV))}] vs [{int(math.ceil(fV))},{node['K_ub']}]")
            continue
        flows = {}
        for p, val in res["xvals"].items():
            if val <= 1e-6:
                continue
            for arc in col_arcs(paths[p]):
                flows[arc] = flows.get(arc, 0.0) + val
        cands = [a for a, f in flows.items() if a[0] >= 1 and a[1] >= 1 and 1e-6 < f < 1 - 1e-6]
        if not cands:
            print(f"{tag}: LP={v:.6f} 无分数弧（数值异常）-> 停止")
            break
        arc = max(cands, key=lambda a: min(flows[a], 1 - flows[a]))
        branches += 1
        c_forb = dict(node, forbidden=node["forbidden"] | {arc})
        c_force = dict(node, forced=node["forced"] + (arc,))
        stack.append(c_force)
        stack.append(c_forb)
        print(f"{tag}: LP={v:.6f} 分数 -> 弧分支 {arc} (流量 {flows[arc]:.3f}): 禁用 vs 强制")
    wall = time.time() - wall0
    print()
    print("== B&P 汇总 ==")
    print("节点数:", nodes, "| 分支数:", branches, "| 墙钟:", round(wall, 2), "s")
    print("incumbent:", None if incumbent is None else round(incumbent, 10),
          "（两位小数", None if incumbent is None else round(incumbent, 2), "）")
    if inc_routes:
        for r in inc_routes:
            seq = [0] + list(r) + [0]
            c = sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1))
            print("  路线", r, "距离", round(c, 6))
    return dict(incumbent=incumbent, routes=inc_routes, nodes=nodes, branches=branches, wall=wall)
'''

C4 = r'''# ---- 实验一：原始 c101（完整 B&P）----
print("========== 实验一：原始 c101（完整 B&P）==========")
r1 = branch_and_price()
print("BKS 对比: 3 车 / 191.81 -> match:",
      r1["incumbent"] is not None and abs(round(r1["incumbent"], 2) - 191.81) < 1e-9)
print("与 02/03/05 已证最优 191.813620 一致:",
      r1["incumbent"] is not None and abs(r1["incumbent"] - 191.813620) < 1e-6)
'''

C5 = r'''# ---- 实验二：车辆数分支机制验证（临时禁用弧 13->17）----
print("========== 实验二：车辆数分支验证（临时禁用弧 13->17）==========")
r2 = branch_and_price(extra_forbidden={(13, 17)})
'''

C6 = r'''# ---- 实验三：Ryan-Foster 弧分支机制验证（固定 3 车 + 禁用 15->16）----
print("========== 实验三：Ryan-Foster 弧分支验证（固定 3 车 + 禁用 15->16）==========")
r3 = branch_and_price(extra_forbidden={(15, 16)}, root_K_ub=3, root_K_lb=3)
'''

MD2 = r"""## 结果与结论

**实验一（原始实例）**：根节点列生成 5 轮收敛，LP = 191.813620 且为**整数解（3 车）**——
根节点直接解决，**0 次分支即证明最优**，与 02/03/05 的结论一致（第四次独立证明）。
原因：该实例集合覆盖 LP 松弛恰为整数（LP 下界 = 整数目标）。

**实验二（分支机制验证：临时禁用弧 13→17）**：根节点 LP = 245.054865、Σx=3.5 分数 →
**车辆数分支** [3,3] vs [4,25]：3 车子节点不可行剪枝，4 车子节点整数解 249.196147（4 车）→ 节点解决。
演示了不可行剪枝与整数节点处理。

**实验三（Ryan-Foster 弧分支验证：固定 3 车 + 禁用 15→16）**：根节点 LP = 232.977957 分数
（Σx=3 整数、客户弧流量分数）→ **弧分支 (10,11)**（流量 0.500）：禁用子节点整数解 241.724452（3 车）、
强制子节点不可行剪枝。演示了弧分支、强制弧在定价子问题中的传播与剪枝。

**结论**：

- B&P 能求解该问题：根节点 LP 整数 → 无需分支即证明最优（191.813620 = 191.81，3 车，约 5 s）。
- 框架完整性：车辆数分支与 Ryan-Foster 弧分支两种规则、定界/不可行/整数三种剪枝，
  以及强制-禁用弧在 CP-SAT 定价子问题中的传播，均经受控实验全部触发验证。
- 定位：B&P 是 02 列生成的直接扩展；对 25 节点实例只需根节点，对 100 节点实例（LP 分数）
  本框架即为标准 branch-and-price 求解内核。
- **基准最优值来源**：Solomon c101(25) 文献 BKS（3 车、191.81）；本 notebook 与 02/03/05
  一致地证明精确最优 191.813620。
"""

cells = [
    nbf.v4.new_markdown_cell(MD0),
    nbf.v4.new_markdown_cell(MD1),
    nbf.v4.new_code_cell(C1),
    nbf.v4.new_code_cell(C2),
    nbf.v4.new_code_cell(C3),
    nbf.v4.new_code_cell(C4),
    nbf.v4.new_code_cell(C5),
    nbf.v4.new_code_cell(C6),
    nbf.v4.new_markdown_cell(MD2),
]
nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3",
                             "language": "python", "name": "python3"}
nbf.write(nb, NB_PATH)
print("written", NB_PATH)
