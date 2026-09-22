# -*- coding: utf-8 -*-
"""构建 csp 06 统一报告 + theory 总览。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/csp_1d_waescher"
CORE = open(BASE + "/scripts/csp_core.py", encoding="utf-8").read()
CORE = CORE[:CORE.index('if __name__ == "__main__":')].rstrip() + "\n"

MD06A = r"""# 1D 切割下料 TEST0005 —— 三方法统一报告（列生成 / Benders / 拉格朗日）

基准：**28 辊，已证明**（长度下界 28 + 池 MIP 28 辊解）。本 notebook 复用 csp_core
运行三种方法并汇总对比。
"""
C06 = r'''# 核心实现 + 三方法
'''
C06 += CORE
C06 += r'''
print("=" * 60)
print("方法一：列生成（背包 DP 定价）")
print("=" * 60)
r1 = run_cg(verbose=True)
print()
print("=" * 60)
print("方法二：Benders")
print("=" * 60)
r2 = run_benders(verbose=True)
print()
print("=" * 60)
print("方法三：拉格朗日松弛")
print("=" * 60)
r3 = run_lagrangian(verbose=True, max_iter=800)
print()
print("=" * 80)
print("统一对比（基准最优 28 辊）")
print("=" * 80)
print(f"{'方法':<14}{'下界':>12}{'整数解':>10}  gap")
print(f"{'列生成':<14}{round(r1['lp'],4):>12}{r1['ip']:>10}  {round((r1['ip']-r1['lp'])/r1['lp']*100,4)}%")
print(f"{'Benders':<14}{round(r2['lb'],4):>12}{r2['ip']:>10}  {round((r2['ip']-r2['lb'])/r2['lb']*100,4)}%")
print(f"{'拉格朗日':<14}{round(r3['lb'],4):>12}{r3['ip']:>10}  {round((r3['ip']-r3['lb'])/r3['lb']*100,4)}%")
print()
print("三方法整数解一致(28):", r1["ip"] == r2["ip"] == r3["ip"] == 28)
'''

MD06B = r"""## 统一报告与结论

- **列生成**：LP 下界 27.9942，整数恢复 28（gap 0.02%）——与长度下界 28 共同证明最优。
- **Benders**：下界 27.9942，修复 28。
- **拉格朗日**：对偶下界 →27.87+（600~800 轮），修复 28。
- 该实例 LP 松弛 = 27.9942 < 28（有积分性间隙 0.02%），最优 28 由「长度下界 + 28 辊解」双重锁定；
  28 解极难直接构造（CP-SAT 560s 未果），但列生成生成的模式池使其自然浮现——本家族是
  **列生成对 1D 切割下料的经典主场**。
**基准最优值来源**：本家族 01_direct（长度下界 + 池 MIP 28 + 覆盖校验）。
"""

# -*- coding: utf-8 -*-
"""构建 theory/csp_对偶推导总览.ipynb。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/csp_1d_waescher/theory"
os.makedirs(BASE, exist_ok=True)

MD = r"""# 1D 切割下料 TEST0005 —— 对偶/割推导总览（Gilmore-Gomory 结构）

问题：$\min\{\sum_p x_p:\ \sum_p a_{ip}x_p\ge d_i,\ x_p\in\mathbb{Z}_+\}$，列=切割模式
（$\sum l_i a_i\le 10000$），最优 28（长度下界 28 = 28 辊解）。

## 01 直接建模 —— 无对偶
证明链：总长度下界 28 + 池 MIP 28 辊解（覆盖校验）；CP-SAT item→bin 直接模型 K=28 难构造（UNKNOWN 560s）。

## 02 列生成 —— LP 对偶
RMP（min Σx）：对偶 $\max\sum_i d_i\pi_i$ s.t. $\sum_i a_{ip}\pi_i\le 1\ \forall$模式, $\pi_i\ge0$。
reduced cost $rc_p=1-\sum_i a_{ip}\pi_i$；定价 = 0/1 背包 $\max\sum\pi_i a_i$（numpy DP）；
无负 rc ⇒ LP 最优（27.9942，对偶=原值强对偶）。

## 03 Benders —— 子问题对偶 → 割
SP(y) 对偶：$\max\sum_i d_i\pi_i+\sum_p\sigma_p(M y_p)$，s.t. $\sum_i a_{ip}\pi_i+\sigma_p\le 1$；
弱对偶 ⇒ 割 $\theta+\sum\lambda_p y_p\ge\sum d_i\pi_i$（λ=−σ）。

## 04 拉格朗日 —— 乘子对偶 + 次梯度
$L(\lambda)=\sum_i d_i\lambda_i+M\cdot\min(0,\ 1-v(\lambda))$，$v(\lambda)$=背包最优值；
子问题=背包（无需求解器）；次梯度 $g_i=d_i-\sum_p a_{ip}x_p$；对偶 → LP 27.9942。
（注意 x 需有界 M=28，否则 rc<0 时对偶函数 −∞。）

## 05 LBBD —— 逻辑割（1D 退化）
主问题 item→bin；子问题=容量检查（线性）；超载箱回传惰性容量割
$\sum_{i\in S}x_{ib}\le|S|-1$（最小超载核心）——1D 下 LBBD 退化为惰性约束主问题（CONVENTIONS §4.5）。

## 07 Branch-and-Price —— 辊数分支
根 LP 27.9942 分数 → 分支 Σx≤27（长度下界不可行）vs Σx≥28（节点整数恢复 28）；
GG pair 分支框架（forbid=双 DP 排除端点 / merge=超件合并）已实现于 scripts/csp_bnp.py。
"""

CODE = r'''# 数值验证：RMP 对偶 + 强对偶 + 定价值 = 1（对偶可行）
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, platform
import ortools
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/csp_1d_waescher/scripts")
import csp_core as cc
from ortools.math_opt.python import mathopt
print("python", platform.python_version(), "| ortools", ortools.__version__)

lp, patterns, sel, iters = cc.cg_min_rolls()
# 最终对偶：重解 RMP 取 π
mm = mathopt.Model()
x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in sel]
covers = []
for i in range(cc.M):
    covers.append(mm.add_linear_constraint(
        mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= cc.DEM[i], name=f"c{i}"))
mm.minimize(mathopt.fast_sum(x))
res = mathopt.solve(mm, mathopt.SolverType.GLOP)
dv = res.dual_values()
pi = [max(0.0, dv[covers[i]]) for i in range(cc.M)]
dual_obj = sum(pi[i]*cc.DEM[i] for i in range(cc.M))
print(f"LP 目标 = {round(lp,6)} | 对偶目标 Σdπ = {round(dual_obj,6)} | 强对偶: {abs(lp-dual_obj)<1e-6}")
v, pat = cc.knap_rebuild(pi)
print(f"定价背包最优值 v(π) = {round(v,6)} | 对偶约束 Σaπ<=1 满足: {v <= 1 + 1e-9}")
rc = 1.0 - v
print(f"最负 rc = {round(rc,8)}（≈0 ⇒ 无改进列 ⇒ LP 最优）")
'''

nb6 = nbf.v4.new_notebook()
nb6.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb6.cells.append(nbf.v4.new_markdown_cell(MD06A))
nb6.cells.append(nbf.v4.new_code_cell(C06))
nb6.cells.append(nbf.v4.new_markdown_cell(MD06B))
nbf.write(nb6, "/mnt/d/exactTest/column-generation-solvers/csp_1d_waescher/06_unified_report.ipynb")
print("written 06_unified_report.ipynb")

nbt = nbf.v4.new_notebook()
nbt.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nbt.cells.append(nbf.v4.new_markdown_cell(MD))
nbt.cells.append(nbf.v4.new_code_cell(CODE))
nbf.write(nbt, "/mnt/d/exactTest/column-generation-solvers/csp_1d_waescher/theory/csp_对偶推导总览.ipynb")
print("written theory/csp_对偶推导总览.ipynb")
