# -*- coding: utf-8 -*-
"""构建 csp_1d_waescher 家族全部 notebook。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/csp_1d_waescher"
os.makedirs(BASE + "/theory", exist_ok=True)
CORE = open(BASE + "/scripts/csp_core.py", encoding="utf-8").read()
# 去掉 __main__ 尾
CORE = CORE[:CORE.index('if __name__ == "__main__":')].rstrip() + "\n"

MD0 = r"""# 1D 切割下料（Wäscher & Gau，TEST0005）—— {METHOD}

## 问题定义

57 种件型、共 114 件，标准辊长 $L=10000$；需求 $d_i$。目标：**用最少标准辊数**满足全部需求。
集合覆盖（Gilmore-Gomory）模型：列 = 一个切割模式（$\sum l_i a_i \le 10000$），$c_p=1$（每辊）：

$$\min_x \sum_p x_p \quad \text{s.t.}\quad \sum_p a_{ip}x_p \ge d_i\ (\forall i),\quad x_p\in\mathbb{Z}_+$$

**基准最优：28 辊，已证明**（总长度下界 $\lceil 279935/10000\rceil=28$ + CG 池 MIP 28 辊解，
覆盖/件数/长度校验通过）。模式总数 2.8 亿 ⇒ 不可枚举，定价必须用背包 DP。
（FFD=29；CP-SAT K=28 直接搜索 560s 未能找到 28 解——28 解存在但极难直接构造，见 02/07。）
"""

MD2 = r"""## 运行结果与结论

见上方输出。结论：**最优 28 辊已证明**（长度下界 28 = 池 MIP 28 解）；
LP 松弛 27.9942（分数），整数最优 28。
基准最优值来源：本家族 01_direct（长度下界 + 池 MIP 28 解 + 覆盖校验）。
"""

def build(name, method_title, md1, call, extra_md=""):
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb.cells.append(nbf.v4.new_markdown_cell(MD0.replace("{METHOD}", method_title)))
    nb.cells.append(nbf.v4.new_code_cell("# 核心实现（数据/下界/背包定价/CG/池MIP 等，本家族 scripts/csp_core.py 同源）\n" + CORE))
    nb.cells.append(nbf.v4.new_markdown_cell(md1))
    nb.cells.append(nbf.v4.new_code_cell(call))
    if extra_md:
        nb.cells.append(nbf.v4.new_markdown_cell(extra_md))
    nb.cells.append(nbf.v4.new_markdown_cell(MD2))
    nbf.write(nb, os.path.join(BASE, name))
    print("written", name)

MD_CG = r"""## 方法：列生成（Gilmore-Gomory，背包 DP 定价）

- **RMP（GLOP LP）**：初始列=57 条单件满复制模式；对偶 $\pi_i\ge 0$；
  reduced cost：$rc_p = 1 - \sum_i a_{ip}\pi_i$。
- **定价子问题 = 0/1 背包**：$\max \sum_i \pi_i a_i$ s.t. $\sum_i l_i a_i\le 10000,\ a_i\le d_i$；
  numpy 向量化 DP（10ms/次，快照重建保证模式合法）。
- **整数恢复**：生成模式池 MIP（HIGHS）→ 28 辊（覆盖校验通过）。
"""
build("02_column_generation.ipynb", "列生成（背包 DP 定价）", MD_CG,
      "\nr = run_cg(verbose=True)\n")

MD_BD = r"""## 方法：Benders（选列主问题 + 覆盖 LP 子问题）

- **主问题**：$\min\theta$ + 最优性割 $\theta+\sum_p\lambda_p^k y_p\ge\sum_i\pi_i^kd_i$（y 选模式）。
- **子问题 SP(y)**：$\min\sum_p x_p+M\sum_i s_i$ s.t. 覆盖$+s\ge d$、$x_p\le M\cdot y_p$（$M=28$）；
  对偶 $\lambda=-\sigma$ 生成割。列成本(=1)含在子问题中。
"""
build("03_benders.ipynb", "Benders 分解", MD_BD,
      "\nr = run_benders(verbose=True)\n")

MD_LG = r"""## 方法：拉格朗日松弛（松弛需求约束）

- 对偶函数（子问题 = 背包，$x_p$ 有界 $M=28$）：
  $$L(\lambda)=\sum_i d_i\lambda_i + M\cdot\min\Big(0,\ 1-\max_{a:\ \sum l a\le L}\sum_i\lambda_i a_i\Big)$$
- 次梯度：$g_i = d_i - \sum_p a_{ip}x_p(\lambda)$（$x=M$ 于最负 rc 模式）；步长
  $\alpha=\rho(UB-L)/\|g\|^2$（步长/λ 截断；锚点 UB=28）。
- 对偶收敛到 LP 松弛 27.9942；池 MIP 修复 28。
"""
build("04_lagrangian.ipynb", "拉格朗日松弛", MD_LG,
      "\nr = run_lagrangian(verbose=True, max_iter=800)\n")

MD_LB = r"""## 方法：LBBD（逻辑 Benders）

- **主问题**：item→bin 分配 MIP（HIGHS，28 箱，每箱负载∈[9935,10000]——总空隙 65 的必然推论）。
- **子问题**：容量检查（线性）；超载箱回传**惰性容量割**（最小超载核心）：
  $\sum_{i\in S} x_{ib} \le |S|-1$。1D 情形 LBBD 退化为带惰性约束的主问题（CONVENTIONS §4.5 的合理形式）。
- 说明：28 分配极难直接搜索（主问题在限时内未自行找到）；其存在性由 CG 池 MIP 证明
  （warm-start 视角），LBBD 割机制如上展示。
"""
build("05_lbbd.ipynb", "LBBD（惰性容量割）", MD_LB,
      "\nr = run_lbbd(verbose=True)\n")

MD_BP = r"""## 方法：Branch-and-Price（辊数分支）

- **根节点** K∈[27,28]：LP=27.9942 分数 → **辊数分支**：左子 Σx≤27（不可行：总长度
  279935>27×10000，虚拟列剪枝）、右子 Σx≥28。
- **右子** K=28：LP=28.0 → 节点整数恢复（池 MIP + 全局 28 模式注入）→ 28 辊整数解。
- 结论：**左子不可行 + 右子整数 28 ⇒ 证明最优**（长度下界 28）。
- 完整 Gilmore-Gomory pair 分支框架见 scripts/csp_bnp.py（本会话深树探索 8 节点未完成，
  辊数分支已给出决定性证明）。
"""
build("07_branch_and_price.ipynb", "Branch-and-Price", MD_BP,
      "\nr = run_bnp(verbose=True)\n")

MD_DR = r"""## 方法：直接建模基准

- **最优性证明链**：①总长度下界 $\lceil279935/10000\rceil=28$；②CG 生成模式池上的
  集合覆盖整数规划（HIGHS，0.02~20s）得 28 辊解；③覆盖/件数(114)/总长度(279935)三重校验通过
  ⇒ 28 为最优。
- CP-SAT item→bin 直接模型：K=28 为 UNKNOWN（560s 未找到——28 解存在但难直接构造）、
  K=29 OPTIMAL（早期上界，已被 28 取代）；如实记录。
"""
build("01_direct.ipynb", "直接建模基准", MD_DR,
      "\nr = run_direct(verbose=True)\n")
