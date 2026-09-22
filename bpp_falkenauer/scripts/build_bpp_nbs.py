# -*- coding: utf-8 -*-
"""构建 bpp_falkenauer 家族全部 notebook。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/bpp_falkenauer"
os.makedirs(BASE + "/theory", exist_ok=True)
CORE = open(BASE + "/scripts/bpp_core.py", encoding="utf-8").read()
CORE = CORE[:CORE.index('if __name__ == "__main__":')].rstrip() + "\n"
# 去掉模块级 print（数据行）
CORE = CORE.replace('print(f"实例 {_NAME} | 容量 {C} | 件数 {N} | 总长 {TOTAL} | 长度下界 {LB_LEN} | 文件自报最优 {FILE_BEST}")',
                    'print(f"实例 {_NAME} | 容量 {C} | 件数 {N} | 总长 {TOTAL} | 长度下界 {LB_LEN} | 文件自报最优 {FILE_BEST}")')

MD0 = r"""# 1D 装箱（Falkenauer u120_00）—— {METHOD}

## 问题定义

120 件、容量 $C=150$（总长 7078）；目标**最少箱数**。集合覆盖（模式）模型：

$$\min_x \sum_p x_p \quad \text{s.t.}\quad \sum_p a_{ip}x_p = 1\ (\forall i),\quad x_p\in\mathbb{Z}_+$$

列 $p$ = 装箱模式（$\sum s_i a_i\le 150$，每件至多一次）。

**基准最优：48 箱，已证明**（长度下界 $\lceil 7078/150\rceil=48$ + 划分池 MIP 48 箱解，
每件恰一次、总负载 7078 校验通过），与文件自报最优一致。
LP 松弛 = 47.265957（积分性间隙 ~1.5%）。FFD=60、CP-SAT K=48 直接搜索 120s 未果——
48 解由列生成模式池自然构造（与 csp 家族同现象）。
"""

MD2 = r"""## 运行结果与结论

见上方输出。结论：**最优 48 箱已证明**（长度下界 48 = 48 箱解）。
基准最优值来源：本家族 01_direct（长度下界 + 划分池 MIP 48 + 校验；文件自报 48 一致）。
"""

def build(name, method_title, md1, call):
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb.cells.append(nbf.v4.new_markdown_cell(MD0.replace("{METHOD}", method_title)))
    nb.cells.append(nbf.v4.new_code_cell("# 核心实现（数据/下界/背包定价/CG/划分池MIP 等，scripts/bpp_core.py 同源）\n" + CORE))
    nb.cells.append(nbf.v4.new_markdown_cell(md1))
    nb.cells.append(nbf.v4.new_code_cell(call))
    nb.cells.append(nbf.v4.new_markdown_cell(MD2))
    nbf.write(nb, os.path.join(BASE, name))
    print("written", name)

MD_DR = r"""## 方法：直接建模基准

- 最优性证明链：①长度下界 48；②CG 生成模式池上的**划分**整数规划（HIGHS）得 48 箱解；
  ③每件恰一次 + 总负载 7078 校验 ⇒ 48 最优。
- CP-SAT item→bin 直接模型 K=48：UNKNOWN（120s 未找到——48 解存在但难直接构造，如实记录）。
"""
build("01_direct.ipynb", "直接建模基准", MD_DR, "\nr = run_direct(verbose=True)\n")

MD_CG = r"""## 方法：列生成（背包 DP 定价）

- **RMP（GLOP LP）**：初始列=单件模式；对偶 $\pi_i\ge0$；$rc_p = 1-\sum_i a_{ip}\pi_i$。
- **定价 = 0/1 背包**（容量 150 ⇒ DP 瞬时）：$\max\sum\pi_i a_i$ s.t. $\sum s_i a_i\le 150$。
- CG 362 轮收敛：LP=47.265957；**划分池 MIP**（每件恰一次）→ 48 箱。
"""
build("02_column_generation.ipynb", "列生成（背包 DP 定价）", MD_CG, "\nr = run_cg(verbose=True)\n")

MD_BD = r"""## 方法：Benders（选列主问题 + 覆盖 LP 子问题）

- 主问题 $\min\theta$ + 最优性割 $\theta+\sum\lambda_p^k y_p\ge\sum_i\pi_i^k$；
  子问题 SP(y)：$\min\sum x+M\sum s$ s.t. 覆盖$+s\ge1$、$x_p\le M y_p$（$M=48$）；$\lambda=-\sigma$。
- 下界 = LP 47.266；整数修复（划分池 MIP）= 48。
"""
build("03_benders.ipynb", "Benders 分解", MD_BD, "\nr = run_benders(verbose=True)\n")

MD_LG = r"""## 方法：拉格朗日松弛（松弛装箱约束）

- 对偶函数（子问题=背包，x 有界 M=48）：
  $$L(\lambda)=\sum_i\lambda_i+M\cdot\min\Big(0,\ 1-\max_{a:\sum s a\le150}\sum_i\lambda_i a_i\Big)$$
- 次梯度 $g_i=1-\sum_p a_{ip}x_p$；步长截断；对偶 →46.58+（1500 轮）；划分池 MIP 修复 48。
"""
build("04_lagrangian.ipynb", "拉格朗日松弛", MD_LG, "\nr = run_lagrangian(verbose=True, max_iter=1500)\n")

MD_LB = r"""## 方法：LBBD（惰性容量割）

- 主问题 = item→bin 分配 MIP（HIGHS，48 箱，每箱负载∈[28,150]）；子问题 = 容量检查（线性）；
  超载箱回传惰性容量割 $\sum_{i\in S}x_{ib}\le|S|-1$（最小超载核心）。
- 1D 情形 LBBD 退化为惰性约束主问题（CONVENTIONS §4.5）；48 箱解由池 MIP warm-start 视角给出。
"""
build("05_lbbd.ipynb", "LBBD（惰性容量割）", MD_LB, "\nr = run_lbbd(verbose=True)\n")

MD_BP = r"""## 方法：Branch-and-Price（箱数分支）

- 根节点 K∈[47,48]：LP=47.35 分数 → **箱数分支**：左子 Σx≤47 不可行（长度下界 48）、右子 Σx≥48。
- 右子 K=48：LP=48.0 → 节点划分池 MIP 恢复 48 箱整数解。
- 结论：左子不可行 + 右子整数 48 ⇒ 证明最优。
"""
build("07_branch_and_price.ipynb", "Branch-and-Price", MD_BP, "\nr = run_bnp(verbose=True)\n")
