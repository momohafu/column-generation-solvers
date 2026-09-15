# -*- coding: utf-8 -*-
"""构建 crew_scheduling 家族全部 notebook（01-07 + theory 总览）。"""
import nbformat as nbf, os, re

BASE = "/mnt/d/exactTest/column-generation-solvers/crew_scheduling"
os.makedirs(BASE + "/theory", exist_ok=True)
SB = BASE + "/scripts"

def read(p):
    return open(p, encoding="utf-8").read()

# 脚本头（数据+池）与主体拆分
def split_script(fname):
    src = read(os.path.join(SB, fname))
    i = src.index("def col_arcs")
    header = src[:i]
    body = src[i:]
    # 去掉 __main__ 尾
    j = body.index('if __name__ == "__main__":')
    body = body[:j].rstrip() + "\n"
    return header, body

HDR_CG, BODY_CG = split_script("crew_cg.py")
HDR_BD, BODY_BD = split_script("crew_benders.py")
HDR_LG, BODY_LG = split_script("crew_lagrangian.py")
HDR_LB, BODY_LB = split_script("crew_lbbd.py")
HDR_BP, BODY_BP = split_script("crew_bnp.py")
HDR_DR, BODY_DR = split_script("crew_direct.py")

POOL_PRINT = '''
print(f"完整池: {P} 条可行路径（单任务 {sum(1 for s in pseq if len(s)==1)} / 双 {sum(1 for s in pseq if len(s)==2)} / 三 {sum(1 for s in pseq if len(s)==3)}）")
print(f"下界: 时长={math.ceil(sum(f-s for s,f in tasks[1:])/T)} | 已知最优: 27 crew / 3139（K=26 不可行）")
'''

MD0 = r"""# 机组排班（Beasley-Cao csp50）—— {METHOD}

## 问题定义

50 个任务 $i$（固定起止时间 $s_i,f_i$）；时间上限 $T=480$；173 条转移弧 $(i,j,c_{ij})$。
一个 crew 的任务序列须逐对由弧连接（弧列表已编码时间兼容性）且**跨度** $f_{last}-s_{first}\le 480$。
目标：**最少 crew 数，其次最小总转移成本**。集合覆盖模型：

$$\min_x \sum_{p\in P} c_p x_p \quad \text{s.t.}\quad \sum_{p\in P} a_{ip}x_p \ge 1\ (\forall i),\quad \sum_{p\in P} x_p \le K,\quad x_p\in\{0,1\}$$

列 $p$ = 一条可行 crew 调度（任务序列），$c_p$ = 序列转移成本之和。
**基准最优**：27 crew、成本 3139（本家族 01 直接模型证明，K=26 不可行）。
"""

MD1_CG = r"""## 方法：列生成（主问题覆盖 LP + 定价子问题）

- **RMP（GLOP LP）**：初始列=50 条单任务列 + 虚拟列（大 M=1e6，因 50>K=27 需 Phase-I）；
  对偶 $\pi_i\ge 0$（覆盖）、$\mu\le 0$（crew 数）。列 $p$ 的 reduced cost：
  $$rc_p = c_p - \sum_{i\in p}\pi_i - \mu$$
- **定价子问题（CP-SAT）**：单 crew 可选任务路径（0→任务…→N+1→0 闭合回路 + 任务自环），
  跨度 $\le 480$，目标 $\min\ \sum c_{ij}x_{ij}-\sum\pi_i v_i-\mu$；时间单调性自动消子环。
- **精确收敛证书**：完整池仅 266 列可完全枚举，每轮池扫描精确校验（与逐列动态定价等价）。
- **整数恢复**：CG 列池 CP-SAT 集合覆盖整数模型。
"""

MD1_BD = r"""## 方法：Benders（选列主问题 + 覆盖 LP 子问题）

- **主问题 MP（HIGHS MIP）**：$\min\theta$，$y_p\in\{0,1\}$，最优性割
  $\theta+\sum_p\lambda_p^k y_p\ge\sum_i\pi_i^k+K\mu^k$（列成本含在子问题中）。
- **子问题 SP(y)（GLOP LP）**：固定 y：$\min\sum c_px_p+M\sum s_i$，
  s.t. 覆盖 $+s_i\ge1$（$\pi_i\ge0$）、$\sum x_p\le K$（$\mu\le0$）、$0\le x_p\le y_p$（$\sigma_p\le0$）；
  对偶 $\lambda=-\sigma$ 生成割（弱对偶：对偶可行点对任意 y 有效）。
- 候选集 = 完整池（266 列）；整数修复 = 完整池 HIGHS。
"""

MD1_LG = r"""## 方法：拉格朗日松弛（松弛覆盖约束）

- 对偶函数（子问题 = 按 rc 排序取前 K 个负列，无需求解器）：
  $$L(\lambda)=\sum_i\lambda_i+\min_{\Sigma x\le K}\sum_p\Big(c_p-\sum_i a_{ip}\lambda_i\Big)x_p$$
- integrality property ⇒ $\max L(\lambda)$ = 覆盖 LP 松弛值 = **3139**（=整数最优，gap 0）。
- 次梯度 $g_i=1-\sum_p a_{ip}x_p(\lambda)$（覆盖次数）；步长
  $\alpha=\rho(UB-L)/\|g\|^2$（ρ=2 起 30 轮无改进减半；步长/λ 截断）；锚点 UB=3139（已知上界）。
- 修复：贪心补列 + 完整池 CP-SAT MIP 修复。
"""

MD1_LB = r"""## 方法：LBBD（逻辑 Benders 分解）

**主方法（scp41 同款，2 轮收敛证明最优）**：主问题 = 选列 $y_p\in\{0,1\}$，$\min\sum c_py_p$，$\sum y\le K$，
初始无覆盖约束；子问题 = 检查未覆盖任务集合 $U$；逻辑割：对每个未覆盖任务
$$\sum_{p: i\in p} y_p \ge 1$$
（恰为覆盖约束本身，收敛时主问题等价于完整池覆盖 IP ⇒ 最优 3139）。

**机制展示（分配主问题 + 调度子问题）**：主问题 = CP-SAT 任务→crew 分配（每 crew ≤3 任务）+ min θ；
子问题 = CP-SAT 单 crew 调度（跨度 ≤480）→ 不可行回传 no-good
$\sum_{i\in S_c}(1-a_{ic})\ge1$，可行回传 Hooker 割 $\theta\ge\sum_c[c_c-\sum\delta_{ic}(1-a_{ic})]$。
如实说明：50 任务/27 crew 组合空间巨大、no-good 割弱，该变体在迭代上限内难以自收敛——
故收敛证明由主方法给出，分配变体用于展示割的构造机制。
"""

MD1_BP = r"""## 方法：Branch-and-Price

- **节点 RMP（GLOP）**：覆盖 LP + 车辆数上下界（对偶 $\mu_{ub}\le0,\ \mu_{lb}\ge0$）+ 虚拟列；
  定价 rc 用 $\mu_{sum}=\mu_{ub}+\mu_{lb}$；虚拟列取正 ⇒ 节点不可行剪枝。
- **定价（CP-SAT）**：单 crew 可选路径 + 强制/禁用弧约束；池扫描精确校验。
- **分支**：Σx 分数 → 车辆数分支；否则 Ryan-Foster 弧分支（流量最接近 0.5 的转移弧：禁用 vs 强制）。
- 剪枝：定界 / 不可行 / LP 整数解。
"""

MD2_COMMON = r"""## 运行结果与结论

见上方输出。结论：

- **{CLAIM}**
- 基准最优值来源：本家族 01_direct 证明（K=26 不可行 + K=27 最优 3139）；
  各方法结果交叉一致。
"""

def build(name, md0_method, md1, hdr, body, call, md2_claim, extra_cells=None):
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb.cells.append(nbf.v4.new_markdown_cell(MD0.replace("{METHOD}", md0_method)))
    nb.cells.append(nbf.v4.new_code_cell(hdr + POOL_PRINT))
    nb.cells.append(nbf.v4.new_markdown_cell(md1))
    nb.cells.append(nbf.v4.new_code_cell(body + call))
    for kind, content in (extra_cells or []):
        nb.cells.append(nbf.v4.new_markdown_cell(content) if kind == "MD" else nbf.v4.new_code_cell(content))
    nb.cells.append(nbf.v4.new_markdown_cell(MD2_COMMON.replace("{CLAIM}", md2_claim)))
    nbf.write(nb, os.path.join(BASE, name))
    print("written", name)

# 01 direct
build("01_direct.ipynb", "直接建模（完整池集合覆盖 IP）",
r"""## 方法：直接建模基准

- **模型**：完整池上的集合覆盖整数规划（HIGHS）：池=全部 266 条可行 crew 调度（DFS 完全枚举，
  跨度≤480 与弧兼容性检查 ⇒ 与原问题**等价**），故这是原问题的精确整数模型：
  $$\min\sum_p c_px_p\ \ \text{s.t.}\ \sum_p a_{ip}x_p\ge 1,\ \sum_p x_p\le K,\ x_p\in\{0,1\}$$
- **最少 crew 数证明**：K=26 不可行 + K=27 最优（0.02s）；容量/重叠/入度下界（14/10/10）仅作参考。
- CP-SAT 交叉验证同解；说明：紧凑 3-index CP-SAT 草稿（direct_solve.py）在 K=16 即 30s 未决，
  故基准采用等价池模型。
""",
HDR_DR, BODY_DR, "\nres = run_direct(verbose=True)\n",
"K=26 不可行 ⇒ 最少 27 crew；K=27 最优成本 3139（HIGHS 与 CP-SAT 交叉验证一致），作为家族基准。")

# 02 CG
build("02_column_generation.ipynb", "列生成（覆盖 LP + CP-SAT 定价）",
MD1_CG, HDR_CG, BODY_CG, "\nres = run_cg(verbose=True)\n",
"CG 2 轮收敛：LP 下界 = 3139 = 整数恢复目标（27 crew）⇒ 证明最优；CP-SAT 定价返回列与池定价一致。")

# 03 Benders
build("03_benders.ipynb", "Benders（选列主问题 + 覆盖 LP 子问题）",
MD1_BD, HDR_BD, BODY_BD, "\nres = run_benders(verbose=True)\n",
"Benders 3 轮 / 3 条对偶割收敛，下界 = 整数修复 = 3139（27 crew）⇒ 证明最优。")

# 04 Lagrangian
build("04_lagrangian.ipynb", "拉格朗日松弛",
MD1_LG, HDR_LG, BODY_LG, "\nres = run_lagrangian(verbose=True)\n",
"次梯度 ~300 轮收敛：对偶下界 = 3139，贪心修复与 MIP 修复均为 3139 ⇒ gap 0.0%，证明最优。")

# 05 LBBD
build("05_lbbd.ipynb", "LBBD（逻辑 Benders 分解）",
MD1_LB, HDR_LB, BODY_LB,
"\nprint('== 主方法：选列主问题 + 覆盖逻辑割 ==')\nres = run_lbbd(verbose=True)\nprint()\nprint('== 机制展示：分配主问题 + 调度子问题 ==')\nresA = run_lbbd_assignment(verbose=True)\n",
"主方法 2 轮 / 50 条覆盖逻辑割收敛，主问题等价完整池 IP ⇒ 证明最优 3139；分配变体如实展示 Hooker/no-good 割机制及其收敛局限。")

# 07 B&P
build("07_branch_and_price.ipynb", "Branch-and-Price（分支-定价）",
MD1_BP, HDR_BP, BODY_BP,
"\nprint('== 实验一：原始 csp50 ==')\nr1 = run_bnp()\nprint()\nprint('== 实验二：禁用最优弧 (1,10) ==')\nr2 = run_bnp(extra_forbidden={(1, 10)})\nprint()\nprint('== 实验三：不可行剪枝演示（禁用弧 (7,20)）==')\nr3 = run_bnp(extra_forbidden={(7, 20)})\n",
"根节点 LP 即整数（3139，27 crew）⇒ 0 分支证明最优；该实例所有受限节点 LP 均为整数/不可行（无分数节点，弧分支规则虽实现但无需触发），分支树机制经受限实验验证。")

print("notebooks built")
