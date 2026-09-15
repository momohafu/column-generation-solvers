# -*- coding: utf-8 -*-
"""构建 crew_scheduling 的 06 统一报告 + theory 总览 notebook。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/crew_scheduling"

MD06A = r"""# 机组排班 csp50 —— 三方法统一报告（列生成 / Benders / LBBD）

问题与基准同各方法 notebook（50 任务、T=480、弧成本；最少 crew 数→最小总成本；
**基准：27 crew / 3139，01 直接模型证明**）。本 notebook 复用 scripts 的
crew_cg / crew_benders / crew_lbbd 运行三种方法，汇总对比并验证一致性。
"""

C06A = r'''# 环境信息
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, platform, time
import ortools
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts")
import crew_cg, crew_benders, crew_lbbd
print("python", platform.python_version(), "| ortools", ortools.__version__)
'''

C06B = r'''# ---- 方法一：列生成 ----
print("=" * 60)
print("方法一：列生成（覆盖 LP + CP-SAT 定价）")
print("=" * 60)
res_cg = crew_cg.run_cg(verbose=True)

# ---- 方法二：Benders ----
print()
print("=" * 60)
print("方法二：Benders（选列主问题 + 覆盖 LP 子问题）")
print("=" * 60)
res_bd = crew_benders.run_benders(verbose=True)

# ---- 方法三：LBBD（主方法）----
print()
print("=" * 60)
print("方法三：LBBD（选列主问题 + 覆盖逻辑割）")
print("=" * 60)
res_lb = crew_lbbd.run_lbbd(verbose=True)
'''

C06C = r'''# ---- 统一汇总与一致性验证 ----
def route_set(routes):
    if isinstance(routes, dict):
        routes = list(routes.values())
    return set(frozenset(r) for r in routes)

s_cg = route_set(res_cg["routes"])
s_bd = route_set(res_bd["routes"])
s_lb = route_set(res_lb["routes"])
print("三方法 crew 集合一致:", s_cg == s_bd == s_lb, "| crew 数:", len(s_cg))
print()
print("=" * 90)
print("统一结果对比（基准最优：27 crew / 3139）")
print("=" * 90)
rows = [
    ("02 列生成", res_cg["ip"], f"{res_cg['iterations']} 轮", f"LP 下界 {round(res_cg['lp'],4)}"),
    ("03 Benders", res_bd["ip"], f"{res_bd['iterations']} 轮 / {res_bd['cuts']} 割", f"下界 {round(res_bd['lb'],4)}"),
    ("05 LBBD", res_lb["obj"], f"{res_lb['iterations']} 轮 / {res_lb['cuts']} 逻辑割", "主问题=完整池 IP"),
]
print(f"{'方法':<12}{'目标':>10}{'迭代/割':>22}  证明机制")
for name, obj, itc, mech in rows:
    print(f"{name:<12}{obj:>10}{itc:>22}  {mech}")
print()
print("与基准 3139 一致:", res_cg["ip"] == res_bd["ip"] == res_lb["obj"] == 3139)
print("最优 crew 路线（三方法一致，示例前 8 条）:")
for r in sorted(s_cg, key=lambda s: -len(s))[:8]:
    print("  ", tuple(r))
'''

MD06B = r"""## 统一报告与结论

三种方法各自独立运行并**一致得到 27 crew / 3139**：

- **02 列生成**：2 轮收敛，LP 下界 = 3139 = 整数恢复（该实例覆盖 LP 松弛恰为整数，与 vrptw 同现象）。
- **03 Benders**：3 轮 / 3 条对偶最优性割，下界 = 整数修复 = 3139。
- **05 LBBD**：2 轮 / 50 条覆盖逻辑割，主问题等价完整池 IP ⇒ 证明最优。

最优性由 01 直接模型（K=26 不可行 + K=27 最优）与三种分解方法多重确认。
**基准最优值来源**：本家族 01_direct 自证（OR-Library csp50 文献实例，此前草稿 pool_opt 亦为 27/3139）。
"""

# -*- coding: utf-8 -*-
"""构建 theory/crew_对偶推导总览.ipynb。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/crew_scheduling/theory"
os.makedirs(BASE, exist_ok=True)

MD = r"""# 机组排班 csp50 —— 六方法对偶/割推导总览（与 vrptw 同构的覆盖结构）

问题：集合覆盖 $\min\{\sum_p c_px_p:\ \sum_p a_{ip}x_p\ge 1,\ \sum_p x_p\le 27,\ x_p\in\{0,1\}\}$，
列 = 可行 crew 调度（完整池 266 条）。基准最优 27 crew / 3139，且 **LP 松弛 = 整数最优**。

## 01 直接建模 —— 无对偶
完整池集合覆盖 IP（HIGHS）直接分支定界；最优性证书 = K=26 不可行 + K=27 目标=界。

## 02 列生成 —— LP 对偶
对偶规则：覆盖（≥）→ $\pi_i\ge0$；crew 数（≤）→ $\mu\le0$；$x_p\ge0$ → $\sum_i a_{ip}\pi_i+\mu\le c_p$。
$$\max\ \sum_i\pi_i+27\mu\ \ \text{s.t.}\ \sum_i a_{ip}\pi_i+\mu\le c_p\ (\forall p),\ \pi_i\ge0,\ \mu\le0$$
reduced cost $rc_p=c_p-\sum_{i\in p}\pi_i-\mu$ = 对偶约束松弛量；定价 = 找最违反的对偶约束
（单 crew 最短路径 + 跨度 ≤480，CP-SAT）；无负 rc 列 ⇒ LP 最优。数值：LP=3139=整数。

## 03 Benders —— 子问题对偶 → 割
SP(y) 对偶：$\max\sum_i\pi_i+27\mu+\sum_p\sigma_py_p$，s.t. $\sum_i a_{ip}\pi_i+\mu+\sigma_p\le c_p,\ \pi_i\le M$；
弱对偶 ⇒ 割 $\theta+\sum_p\lambda_py_p\ge\sum_i\pi_i+27\mu$（$\lambda=-\sigma\ge0$）对任意 y 有效。

## 04 拉格朗日 —— 乘子对偶 + 次梯度
$L(\lambda)=\sum_i\lambda_i+\min_{\Sigma x\le27}\sum_p(c_p-\sum_i a_{ip}\lambda_i)x_p$（子问题=排序取前 27 负列）；
integrality property ⇒ 对偶 = LP 松弛 = 3139；次梯度 $g_i=1-\sum_p a_{ip}x_p$（覆盖次数），
支撑不等式由子问题最优性直接导出。

## 05 LBBD —— 逻辑割（无需对偶）
主方法：选列主问题 + 未覆盖任务检查 → 覆盖逻辑割 $\sum_{p\ni i}y_p\ge1$（= 原覆盖约束，惰性约束 MIP）；
分配变体：调度子问题不可行 → no-good $\sum_{i\in S_c}(1-a_{ic})\ge1$；可行 → Hooker 割
$\theta\ge\sum_c[c_c-\sum\delta_{ic}(1-a_{ic})]$（δ=精确单任务移除节省）。

## 07 Branch-and-Price —— 节点对偶
节点 RMP 加 crew 数上下界：下界（≥）对偶 $\mu_{lb}\ge0$、上界（≤）对偶 $\mu_{ub}\le0$；
$rc_p=c_p-\sum\pi-\mu_{ub}-\mu_{lb}$；弧流量 $f_{ij}=\sum_px_p\mathbb1[(i,j)\in p]$ 分数 ⇒ 分支。
本实例所有节点 LP 均为整数/不可行（无分数节点），根节点即证明最优。
"""

CODE = r'''# 数值验证：RMP 对偶 + 强对偶 + 拉格朗日强对偶点（复用 crew 脚本）
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, platform, math
import ortools
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/crew_scheduling/scripts")
import crew_cg
from crew_cg import N, K_MIN, P, pseq, pcost, pmask
from ortools.math_opt.python import mathopt
print("python", platform.python_version(), "| ortools", ortools.__version__)

sel = list(range(P))
obj, pi, mu, xvals, y_used = crew_cg.solve_rmp(sel)
dual_obj = sum(pi[1:]) + K_MIN*mu
print(f"LP 目标 = {round(obj,4)} | 对偶目标 Σπ+27μ = {round(dual_obj,4)} | 强对偶: {abs(obj-dual_obj)<1e-6}")
print(f"mu = {round(mu,4)} | 非零 π 数 = {sum(1 for v in pi[1:] if v > 1e-6)}")
nz = 0
for p in range(P):
    if xvals[p] > 1e-6:
        s2 = 0.0
        mm = pmask[p]
        while mm:
            lb = mm & -mm
            s2 += pi[lb.bit_length()-1]
            mm -= lb
        rc = pcost[p] - s2 - mu
        assert abs(rc) < 1e-5, (p, rc)
        nz += 1
print(f"绑定列数 {nz}，全部 rc≈0（互补松弛）✓")
lam = pi
neg = sorted((pcost[p]-sum(lam[i] for i in pseq[p]), p) for p in range(P))
S = [p for rc, p in neg[:K_MIN] if rc < -1e-7]
L = sum(lam[1:]) + sum(pcost[p]-sum(lam[i] for i in pseq[p]) for p in S)
print(f"λ=LP对偶π: L(λ) = {round(L,4)} = LP 值 {round(obj,4)} -> {abs(L-obj)<1e-6}")
cnt = [0]*(N+1)
for p in S:
    mm = pmask[p]
    while mm:
        lb = mm & -mm
        cnt[lb.bit_length()-1] += 1
        mm -= lb
g = [0.0]*(N+1)
for i in range(1, N+1):
    g[i] = 1.0 - cnt[i]
rhs = L + sum(g[i]*(0.0-lam[i]) for i in range(1, N+1))
print(f"支撑性质 L(0)=0 >= L(λ)+Σg(0-λ) = {round(rhs,4)} -> {rhs <= 1e-9}")
'''

nb6 = nbf.v4.new_notebook()
nb6.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
for cell in [nbf.v4.new_markdown_cell(MD06A), nbf.v4.new_code_cell(C06A), nbf.v4.new_code_cell(C06B),
             nbf.v4.new_code_cell(C06C), nbf.v4.new_markdown_cell(MD06B)]:
    nb6.cells.append(cell)
nbf.write(nb6, "/mnt/d/exactTest/column-generation-solvers/crew_scheduling/06_unified_report.ipynb")
print("written 06_unified_report.ipynb")

nbt = nbf.v4.new_notebook()
nbt.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nbt.cells.append(nbf.v4.new_markdown_cell(MD))
nbt.cells.append(nbf.v4.new_code_cell(CODE))
nbf.write(nbt, "/mnt/d/exactTest/column-generation-solvers/crew_scheduling/theory/crew_对偶推导总览.ipynb")
print("written theory/crew_对偶推导总览.ipynb")
