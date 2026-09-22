# -*- coding: utf-8 -*-
"""构建 cutting_2d_cgcut 家族全部 notebook。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/cutting_2d_cgcut"
os.makedirs(BASE + "/theory", exist_ok=True)
CORE = open(BASE + "/scripts/cg2_core.py", encoding="utf-8").read()
CORE = CORE[:CORE.index('if __name__ == "__main__":')].rstrip() + "\n"

MD0 = r"""# 2D guillotine 切割（cgcut1）—— {METHOD}

## 问题定义

单张板材 $15\times10$，7 种件型 $(l,w,q,v)$：$(8,4,2,66),(3,7,1,35),(8,2,3,24),(3,4,5,17),(3,3,2,11),(3,2,2,8),(2,1,1,2)$。
**固定朝向（不允许 90° 旋转，文献 cgcut1 约定）**；切割必须为 **guillotine**（贯通的横/竖切）。
目标：**最大化切出件总价值**。模式模型：

$$\max_a \sum_i v_i a_i \quad \text{s.t.}\quad a\in P,\quad a_i\le q_i$$

$P$ = 全部 guillotine 可行模式（$a_i$ = 件型 $i$ 的件数）。
**基准最优：244，已证明**（候选全枚举 2119 个计数向量 + 精确 guillotine 检验递归 → 完整池 1748 模式 → 池 IP 244），
与文献 cgcut1（Christofides-Whitlock 家族）一致。
注意：若允许旋转，可另得 260（本家族固定朝向与文献一致，旋转敏感性已在 01 中说明）。
"""

MD2 = r"""## 运行结果与结论

见上方输出。结论：**最优 244 已证明**（完整池 IP = LP 上界 = 244）。
基准最优值来源：文献 cgcut1 最优 244；本家族 01_direct（完整池枚举 + 精确检验 + 池 IP）自证。
"""

def build(name, method_title, md1, call):
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb.cells.append(nbf.v4.new_markdown_cell(MD0.replace("{METHOD}", method_title)))
    nb.cells.append(nbf.v4.new_code_cell("# 核心实现（数据/精确 guillotine 检验/完整池/六方法，scripts/cg2_core.py 同源）\n" + CORE))
    nb.cells.append(nbf.v4.new_markdown_cell(md1))
    nb.cells.append(nbf.v4.new_code_cell(call))
    nb.cells.append(nbf.v4.new_markdown_cell(MD2))
    nbf.write(nb, os.path.join(BASE, name))
    print("written", name)

MD_DR = r"""## 方法：直接建模基准

- 完整模式枚举：候选计数向量 2119 个（面积 ≤150）→ 精确 guillotine 检验递归（递归切分 + 记忆化）
  → 1748 个可行模式（池**完整**）。
- 池 IP（HIGHS，选一个模式）= **244**（模式 (2,1,1,2,1,1,0)）；池完整性 ⇒ 精确最优。
- 旋转敏感性：允许旋转时最优为 260（文档约定：本家族固定朝向，与文献 244 一致）。
"""
build("01_direct.ipynb", "直接建模基准", MD_DR, "\nr = run_direct(verbose=True)\n")

MD_CG = r"""## 方法：列生成（池定价）

- 主问题 = 模式凸包 LP：$\max \sum_p val_p \lambda_p$ s.t. $\sum\lambda_p=1$——
  线性目标在单纯形上的最大 = 最大模式价值 = **244 = IP**（LP=IP，与 vrptw/crew 同现象）。
- 定价 = 池扫描（池完整 ⇒ 与动态定价等价）；整数恢复 = 池 IP。
"""
build("02_column_generation.ipynb", "列生成（池定价）", MD_CG, "\nr = run_cg(verbose=True)\n")

MD_BD = r"""## 方法：Benders（max 版本）

- 主问题：$\max\theta$，$y_p\in\{0,1\}$，$\sum y\ge1$；最优性割（SP 对偶 $\sigma=\max(0,val-\theta^*)$）：
  $$\theta \le V(y^k)+\sum_p\sigma_p(y_p-y_p^k)$$
- 子问题 SP(y)：受限池上 $\max\sum\lambda_p val_p$，$\sum\lambda=1$，$\lambda\le y$。
- 收敛：$\theta\to\max_p val_p=244$；整数修复 = 池 IP 244。
"""
build("03_benders.ipynb", "Benders 分解", MD_BD, "\nr = run_benders(verbose=True)\n")

MD_LG = r"""## 方法：拉格朗日松弛（松弛件数上限）

- 松弛 $a_i\le q_i$（乘子 $\mu_i\ge0$），对偶函数（子问题 = 池扫描找最大修正价值模式）：
  $$L(\mu)=\sum_i\mu_i q_i+\max_{a\in P}\sum_i(v_i-\mu_i)a_i$$
- 次梯度 $g_i=q_i-a_i^*$（下降方向 $\mu\leftarrow\mu-\alpha g$，锚点 244）；
  $\min L(\mu)$ = LP 值 = 244（强对偶）。
"""
build("04_lagrangian.ipynb", "拉格朗日松弛", MD_LG, "\nr = run_lagrangian(verbose=True)\n")

MD_LB = r"""## 方法：LBBD（逻辑 Benders）

- 主问题：候选计数向量上最大化 $\sum v_i a_i$（面积 ≤150、$a_i\le q_i$，HIGHS，无几何约束）。
- 子问题：**精确 guillotine 检验**（递归切分 + 记忆化）——不可切 → no-good（排除该向量）：
  $$\sum_i z_i\ge1,\quad a_i\le c_i-1+M z_i,\ a_i\ge c_i+1-M(1-z_i)$$
- 迭代至主问题解可行 → 244（首个不可切候选 (2,1,3,0,1,1,1) 被排除，与池检验一致）。
"""
build("05_lbbd.ipynb", "LBBD（逻辑 Benders）", MD_LB, "\nr = run_lbbd(verbose=True)\n")

MD_BP = r"""## 方法：Branch-and-Price（件数分支）

- 根节点 LP(conv) = 244 且为**整数**（单模式）→ **0 分支即证明最优**（与 vrptw 同现象）。
- 分支规则（件数 $a_i\le k$ vs $a_i\ge k+1$，节点池过滤 + LP）已实现；本实例 LP 整数无需触发。
- 结论：LP 上界 244 = 整数解 244 ⇒ 证明最优。
"""
build("07_branch_and_price.ipynb", "Branch-and-Price", MD_BP, "\nr = run_bnp(verbose=True)\n")
