# -*- coding: utf-8 -*-
"""构建 bpp 06 统一报告 + theory 总览。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer"
CORE = open(BASE + "/scripts/bpp_core.py", encoding="utf-8").read()
CORE = CORE[:CORE.index('if __name__ == "__main__":')].rstrip() + "\n"

MD06A = r"""# 1D 装箱 u120_00 —— 三方法统一报告（列生成 / Benders / 拉格朗日）

基准：**48 箱，已证明**（长度下界 48 + 划分池 MIP 48 箱解）。LP 松弛 47.266。
"""
C06 = "# 核心实现 + 三方法\n" + CORE + r'''
print("=" * 60)
print("方法一：列生成（背包 DP 定价 + 划分池 MIP）")
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
r3 = run_lagrangian(verbose=True, max_iter=1500)
print()
print("=" * 80)
print("统一对比（基准最优 48 箱）")
print("=" * 80)
print(f"{'方法':<14}{'下界':>12}{'整数解':>10}  gap")
print(f"{'列生成':<14}{round(r1['lp'],4):>12}{r1['K']:>10}  {round((r1['K']-r1['lp'])/r1['lp']*100,4)}%")
print(f"{'Benders':<14}{round(r2['lb'],4):>12}{int(round(r2['ip'])):>10}  {round((r2['ip']-r2['lb'])/r2['lb']*100,4)}%")
print(f"{'拉格朗日':<14}{round(r3['lb'],4):>12}{int(round(r3['ip'])):>10}  {round((r3['ip']-r3['lb'])/r3['lb']*100,4)}%")
print()
print("三方法整数解一致(48):", r1["K"] == int(round(r2["ip"])) == int(round(r3["ip"])) == 48)
'''

MD06B = r"""## 统一报告与结论

- **列生成**：CG 362 轮，LP=47.265957，划分池 MIP 48（gap ~1.5%）——与长度下界 48 共同证明最优。
- **Benders**：下界 47.266，修复 48。
- **拉格朗日**：对偶下界 46.58+（1500 轮），修复 48。
- 该实例 LP 与整数最优有真实间隙（47.27 vs 48），Falkenauer 经典硬例；
  48 解由列生成模式池自然构造（CP-SAT 直接搜索 120s 未果）。
**基准最优值来源**：本家族 01_direct（长度下界 + 划分池 MIP 48 + 校验；与文件自报 48 一致）。
"""

# -*- coding: utf-8 -*-
"""构建 theory/bpp_对偶推导总览.ipynb。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer/theory"
os.makedirs(BASE, exist_ok=True)

MD = r"""# 1D 装箱 u120_00 —— 对偶/割推导总览（模式覆盖结构）

问题：$\min\{\sum_p x_p:\ \sum_p a_{ip}x_p=1,\ x_p\in\mathbb{Z}_+\}$，列=装箱模式
（$\sum s_i a_i\le150$）。最优 48（长度下界 48 = 48 箱解）；LP 松弛 47.265957。

## 01 直接建模 —— 无对偶
证明链：长度下界 48 + 划分池 MIP 48 箱解（每件恰一次、总负载校验）；CP-SAT K=48 直接搜索 UNKNOWN。

## 02 列生成 —— LP 对偶
RMP 对偶：$\max\sum_i\pi_i$ s.t. $\sum_i a_{ip}\pi_i\le 1\ \forall$模式（划分约束 ⇒ 对偶变量自由，
本实现用覆盖+划分 MIP 恢复）；$rc_p=1-\sum a_{ip}\pi_i$；定价 = 0/1 背包（容量 150）。

## 03 Benders —— 子问题对偶 → 割
SP(y) 对偶：$\max\sum_i\pi_i+\sum_p\sigma_p(M y_p)$，s.t. $\sum a_{ip}\pi_i+\sigma_p\le1$；
弱对偶 ⇒ 割 $\theta+\sum\lambda_p y_p\ge\sum\pi_i$。

## 04 拉格朗日 —— 乘子对偶 + 次梯度
$L(\lambda)=\sum_i\lambda_i+M\min(0,\ 1-v(\lambda))$；子问题=背包；$g_i=1-\sum_p a_{ip}x_p$；
对偶 = LP 松弛 47.266。

## 05 LBBD —— 逻辑割（1D 退化）
主问题 item→bin；子问题容量检查；超载箱 → 惰性容量割（最小超载核心）——
1D 下 LBBD 退化为惰性约束主问题。

## 07 Branch-and-Price —— 箱数分支
根 LP 47.35 分数 → 分支 Σx≤47（长度下界不可行）vs Σx≥48（节点划分池 MIP 恢复 48）。
"""

CODE = r'''# 数值验证：RMP 对偶 + 强对偶 + 定价值=1
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, platform
import ortools
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer/scripts")
import bpp_core as bc
from ortools.math_opt.python import mathopt
print("python", platform.python_version(), "| ortools", ortools.__version__)

lp, patterns, sel, iters = bc.cg_min_rolls(max_iter=1500)
mm = mathopt.Model()
x = [mm.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{k}") for k in sel]
covers = []
for i in range(bc.N):
    covers.append(mm.add_linear_constraint(
        mathopt.fast_sum([x[k]*patterns[k][i] for k in range(len(sel))]) >= 1.0, name=f"c{i}"))
mm.minimize(mathopt.fast_sum(x))
res = mathopt.solve(mm, mathopt.SolverType.GLOP)
dv = res.dual_values()
pi = [max(0.0, dv[covers[i]]) for i in range(bc.N)]
dual_obj = sum(pi)
print(f"LP 目标 = {round(lp,6)} | 对偶目标 Σπ = {round(dual_obj,6)} | 强对偶: {abs(lp-dual_obj)<1e-6}")
v, isel = bc.knap_rebuild(pi)
print(f"定价背包最优值 v(π) = {round(v,6)} | 对偶约束 Σaπ<=1 满足: {v <= 1 + 1e-9}")
print(f"最负 rc = {round(1.0 - v, 8)}（≈0 ⇒ LP 最优）")
'''

nb6 = nbf.v4.new_notebook()
nb6.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb6.cells.append(nbf.v4.new_markdown_cell(MD06A))
nb6.cells.append(nbf.v4.new_code_cell(C06))
nb6.cells.append(nbf.v4.new_markdown_cell(MD06B))
nbf.write(nb6, "/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer/06_unified_report.ipynb")
print("written 06_unified_report.ipynb")

nbt = nbf.v4.new_notebook()
nbt.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nbt.cells.append(nbf.v4.new_markdown_cell(MD))
nbt.cells.append(nbf.v4.new_code_cell(CODE))
nbf.write(nbt, "/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer/theory/bpp_对偶推导总览.ipynb")
print("written theory/bpp_对偶推导总览.ipynb")
