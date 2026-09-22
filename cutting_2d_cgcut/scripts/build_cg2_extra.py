# -*- coding: utf-8 -*-
"""构建 cg2 06 统一报告 + theory 总览。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut"
CORE = open(BASE + "/scripts/cg2_core.py", encoding="utf-8").read()
CORE = CORE[:CORE.index('if __name__ == "__main__":')].rstrip() + "\n"

MD06A = r"""# 2D guillotine 切割 cgcut1 —— 三方法统一报告（列生成 / Benders / 拉格朗日）

基准：**244，已证明**（完整池 IP = LP 上界）。本 notebook 复用 cg2_core 运行三种方法并汇总。
"""
C06 = "# 核心实现 + 三方法\n" + CORE + r'''
print("=" * 60)
print("方法一：列生成（池定价）")
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
r3 = run_lagrangian(verbose=True)
print()
print("=" * 80)
print("统一对比（基准最优 244）")
print("=" * 80)
print(f"{'方法':<14}{'LP/上界':>12}{'整数解':>10}  gap")
print(f"{'列生成':<14}{round(r1['lp'],2):>12}{r1['ip']:>10}  {round((r1['lp']-r1['ip'])/r1['ip']*100,4)}%")
print(f"{'Benders':<14}{round(r2['ub'],2):>12}{r2['ip']:>10}  {round((r2['ub']-r2['ip'])/r2['ip']*100,4)}%")
print(f"{'拉格朗日':<14}{round(r3['ub'],2):>12}{r3['ip']:>10}  {round((r3['ub']-r3['ip'])/r3['ip']*100,4)}%")
print()
print("三方法一致(244):", r1["ip"] == r2["ip"] == r3["ip"] == 244)
'''

MD06B = r"""## 统一报告与结论

- **列生成**：LP(conv) = 244 = 整数解（该单张板材结构下 LP 恰为最大模式价值）。
- **Benders**：割收敛上界 244，修复 244。
- **拉格朗日**：对偶上界 244（第 1 轮即达），修复 244。
- 三种方法一致 244；最优性由 01 的完整池枚举 + 精确 guillotine 检验 + 池 IP 证明。
**基准最优值来源**：文献 cgcut1 最优 244；本家族 01_direct 自证（完整池）。
"""

# -*- coding: utf-8 -*-
"""构建 theory/cg2_对偶推导总览.ipynb。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/theory"
os.makedirs(BASE, exist_ok=True)

MD = r"""# 2D guillotine 切割 cgcut1 —— 对偶/割推导总览（单张板材模式结构）

问题：$\max\{\sum v_i a_i:\ a\in P\}$，$P$=guillotine 可行模式集（固定朝向，完整池 1748 模式）。
最优 244（文献一致；LP = IP）。

## 01 直接建模 —— 无对偶
候选全枚举（2119 计数向量）+ 精确 guillotine 检验递归（切分归纳 + 记忆化）⇒ 池完整 ⇒ 池 IP = 精确解 244。

## 02 列生成 —— LP 即最大模式价值
主问题 = 模式凸包 LP：$\max\sum_p val_p\lambda_p$ s.t. $\sum\lambda_p=1$——线性目标在单纯形上的
最大 = 最大系数 = 244 = IP。定价 = 池扫描（池完整 ⇒ 等价动态定价）。

## 03 Benders（max 版本）—— SP 对偶 → 上界割
SP(y)：$\max\sum\lambda_p val_p$，$\sum\lambda=1$，$\lambda\le y$；对偶 $\min\theta+\sum\sigma_p y_p$
s.t. $\theta+\sigma_p\ge val_p$。取 $\sigma_p=\max(0,val_p-\theta^*)$ ⇒ 割
$\theta\le V(y^k)+\sum\sigma_p(y_p-y^k_p)$（对任意 y 有效）；主问题 $\max\theta$ 收敛到 $\max_p val_p=244$。

## 04 拉格朗日 —— 松弛件数上限
$L(\mu)=\sum\mu_i q_i+\max_{a\in P}\sum(v_i-\mu_i)a_i$（子问题=池扫描）；次梯度 $g=q-a^*$，
下降 $\mu\leftarrow\mu-\alpha g$；$\min L$ = LP = 244（strong duality）。

## 05 LBBD —— 逻辑割
主问题=候选最大化（无几何）；子问题=精确 guillotine 检验；不可切 → no-good 排除该向量
（$|a-c|\ge1$ 线性化）。首轮排除 (2,1,3,0,1,1,1)，收敛到可行 244。

## 07 Branch-and-Price —— 件数分支
根 LP = 244 整数（单模式）→ 0 分支证明最优；分支规则（$a_i\le k$ / $a_i\ge k+1$ + 节点池过滤）已实现。
"""

CODE = r'''# 数值验证：池 IP = LP(conv) = 244
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, platform
import ortools
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/scripts")
import cg2_core as cc
from ortools.math_opt.python import mathopt
print("python", platform.python_version(), "| ortools", ortools.__version__)

term, obj, sel, wt = cc.pool_ip()
lp = cc.conv_lp()
print(f"完整池 {cc.P} 模式 | 池 IP = {obj}（模式 {cc.POOL[sel[0]] if sel else None}）")
print(f"LP(conv) = {round(lp,4)} | LP = IP = 244: {abs(lp-244)<1e-9 and abs(obj-244)<1e-9}")
print(f"最优模式价值复算: {sum(cc.V[i]*cc.POOL[sel[0]][i] for i in range(cc.M))}")
'''

nb6 = nbf.v4.new_notebook()
nb6.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nb6.cells.append(nbf.v4.new_markdown_cell(MD06A))
nb6.cells.append(nbf.v4.new_code_cell(C06))
nb6.cells.append(nbf.v4.new_markdown_cell(MD06B))
nbf.write(nb6, "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/06_unified_report.ipynb")
print("written 06_unified_report.ipynb")

nbt = nbf.v4.new_notebook()
nbt.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
nbt.cells.append(nbf.v4.new_markdown_cell(MD))
nbt.cells.append(nbf.v4.new_code_cell(CODE))
nbf.write(nbt, "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut/theory/cg2_对偶推导总览.ipynb")
print("written theory/cg2_对偶推导总览.ipynb")
