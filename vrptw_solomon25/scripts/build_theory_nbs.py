# -*- coding: utf-8 -*-
"""构建 theory/ 理论详解 notebook 系列（每个方法：完整数学公式 + 对偶推导 + 数值验证）。"""
import nbformat as nbf, os

BASE = "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/theory"
os.makedirs(BASE, exist_ok=True)

PREAMBLE = r'''# 环境与演示数据（28 列小池 = 25 条单客户路径 + 3 条最优路线）
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import sys, platform, math, time
import ortools
sys.path.insert(0, "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/scripts")
from cg_cpsat import build_data
from ortools.math_opt.python import mathopt
print("python", platform.python_version(), "| ortools", ortools.__version__)

n, xc, yc, dem, ready, due, svc, cap, depot_due, dist, d_scaled = build_data()
K = 25
ROUTES = [(13,17,18,19,15,16,14,12), (20,24,25,23,22,21), (5,3,7,8,10,11,9,6,4,2,1)]
pool = [(i,) for i in range(1, n+1)] + ROUTES

def col_cost(p):
    c = 0.0
    prev = 0
    for j in p:
        c += dist(prev, j)
        prev = j
    return c + dist(prev, 0)

def col_mask(p):
    m = 0
    for j in p:
        m |= (1 << j)
    return m

costs = [col_cost(p) for p in pool]
masks = [col_mask(p) for p in pool]
P = len(pool)
print(f"演示池: {P} 列（25 单客户 + 3 最优路线）；最优值 191.813620（=191.81）")
'''

# ================= 01 =================
NB01 = [
 ("MD", r"""# 01 直接建模（3-index CP-SAT）—— 完整数学公式与最优性证书

## 完整模型

**变量**

- $x_{ijk}\in\{0,1\}$：车辆 $k$ 走弧 $(i,j)$，$i,j\in\{0\}\cup C$（$0$=仓库）
- $t_{ik}\in[0,L_0\cdot S]$：车辆 $k$ 在客户 $i$ 的**开始服务时刻**（整数缩放 $S=10000$）
- $q_{ik}\in[0,Q]$：车辆 $k$ 服务完客户 $i$ 后的累计负载

**约束**

1. 每车一个回路：$\text{AddCircuit}(x_{k})$，即
$$\sum_j x_{ijk}=\sum_j x_{jik}=1,\ \forall i,k\qquad(\text{自环允许})$$
2. 每客户恰被一车访问（其余车在该客户自环）：
$$\sum_k x_{iik}=K-1,\ \forall i\in C$$
3. 对称破缺（未使用车辆=后缀，仓库自环 $x_{00k}=1$ 表示车 $k$ 未用）：
$$x_{00k}\le x_{00,k+1}$$
4. 时间窗（允许等待；$M=10^{12}$ 大数，自环 $x_{iik}=1$ 时解除约束）：
$$t_{jk}\ge t_{ik}+s_i S+d_{ij}S-M(1-x_{ijk}),\qquad r_iS\le t_{ik}\le l_iS$$
5. 容量：$q_{jk}\ge q_{ik}+\delta_j-M(1-x_{ijk})$，$0\le q_{ik}\le Q$
6. 返回仓库：$t_{ik}+s_iS+d_{i0}S\le L_0S+M(1-x_{i0k})$
7. 子回路自动消除：弧距离为正 ⇒ 客户环上 $t$ 严格递增矛盾（时间单调性）

**目标**

$$\min\ \sum_{k}\sum_{i\ne j}d_{ij}S\cdot x_{ijk}$$

## 为什么没有"对偶"

- 该模型的回路（AddCircuit）与时间窗约束**没有紧凑的线性化**，LP 松弛/对偶理论不直接适用；
  CP-SAT 用**约束传播 + 分支定界**直接证明最优，**不需要显式对偶**。
- **最优性证书**：求解状态 OPTIMAL ⟺ 可行解（上界）与搜索界（下界）相等：
  $objective = best\_bound = 1918136$（缩放值）。
- 正是这种"无对偶可用"的结构，催生了 02–07 的分解方法（列生成把困难留在定价子问题、
  拉格朗日/LBBD 用乘子与逻辑割构造对偶信息）。
"""),
 ("CODE", r'''from direct_solve import solve
st, obj, routes, el, bound = solve(3, 120)
print(f"K=3: {st} | 缩放目标 {obj} | 下界(best bound) {bound} | 耗时 {round(el,1)}s")
print("最优性证书: 目标 == 下界 ->", obj == bound and st.name == "OPTIMAL")
exact = 0.0
for r in routes:
    seq = [0] + list(r) + [0]
    c = sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1))
    exact += c
    print("  路线", r, "距离", round(c, 6))
print("精确总距离", round(exact, 10), "= 191.81 两位小数（BKS）")
'''),
]

# ================= 02 =================
NB02 = [
 ("MD", r"""# 02 列生成 —— LP 对偶与 reduced cost 的完整推导

## 主问题（RMP，LP 松弛）

$$\min_x \sum_{p\in P'} c_p x_p \quad\text{s.t.}\quad \sum_{p\in P'} a_{ip}x_p\ge 1\ (\forall i\in C),\qquad \sum_{p\in P'} x_p\le K,\qquad x_p\ge 0$$

## 对偶是怎么来的：两步推导

**第一步（机械规则）**。min 问题的 LP 对偶规则：

| 原问题（min） | 对偶问题（max） |
|---|---|
| 约束 $\sum a x \ge b$ | 变量 $y\ge 0$，目标系数 $b$ |
| 约束 $\sum a x \le b$ | 变量 $y\le 0$，目标系数 $b$ |
| 变量 $x\ge 0$ | 约束 $\sum a y \le c$ |

应用：覆盖约束（≥）→ $\pi_i\ge 0$（目标系数 1）；车辆数约束（≤）→ $\mu\le 0$（目标系数 $K$）；
$x_p\ge 0$ → 对偶约束 $\sum_i a_{ip}\pi_i+\mu\le c_p$。得：

$$\max_{\pi,\mu}\ \sum_{i\in C}\pi_i+K\mu \quad\text{s.t.}\quad \sum_{i\in C}a_{ip}\pi_i+\mu\le c_p\ (\forall p),\quad \pi_i\ge 0,\ \mu\le 0$$

**第二步（乘子法，理解"为什么"）**。把约束写成 $\le 0$ 形式并配非负乘子：
$g_i(x)=1-\sum_p a_{ip}x_p\le 0$（乘子 $\pi_i\ge 0$）、$h(x)=\sum_p x_p-K\le 0$（乘子 $\nu\ge 0$），
拉格朗日函数（对 $x\ge 0$）：

$$L(x;\pi,\nu)=\sum_p c_px_p+\sum_i\pi_i\Big(1-\sum_p a_{ip}x_p\Big)+\nu\Big(\sum_p x_p-K\Big)
=\underbrace{\sum_i\pi_i-K\nu}_{\text{常数}}+\sum_p\Big(\underbrace{c_p-\sum_i a_{ip}\pi_i+\nu}_{\text{列 }p\text{ 的系数}}\Big)x_p$$

对 $x\ge 0$ 求 min：若某列系数 $<0$ 则 $x_p\to\infty$ 使 $L\to-\infty$，故有限值要求
$c_p-\sum_i a_{ip}\pi_i+\nu\ge 0\ \forall p$。令 $\mu=-\nu\le 0$ 即得对偶约束
$\sum_i a_{ip}\pi_i+\mu\le c_p$，对偶目标 $\max\sum_i\pi_i+K\mu$。两种推导完全一致。
"""),
 ("MD", r"""## reduced cost = 对偶约束的松弛量

$$rc_p \;=\; c_p-\sum_{i\in C}a_{ip}\pi_i-\mu \;=\; c_p-\sum_{i\in p}\pi_i-\mu$$

三个核心命题：

1. **$rc_p<0$ ⟺ 对偶约束被违反** ⟺ 当前 $(\pi,\mu)$ 对偶不可行 ⟺ 列 $p$ 入基可使 RMP 目标下降（单纯形判据）。
2. **所有列 $rc_p\ge 0$ ⟺ 对偶可行**。由弱对偶（原目标 ≥ 对偶目标）与 RMP 自身的强对偶，
   RMP 解即完整主问题 LP 松弛的最优解 —— 这就是列生成的下界。
3. **互补松弛**：$x_p>0\Rightarrow rc_p=0$；$\pi_i>0\Rightarrow$ 覆盖约束绑定；$\mu<0\Rightarrow$ 车辆数绑定。

**定价子问题** = 找最违反的对偶约束：$\min_p rc_p=\min_p[c_p-\sum_{i\in p}\pi_i-\mu]$，
即带节点收益 $\pi_i$ 的 ESPPRC（CP-SAT 建模见 02 notebook）。
"""),
 ("CODE", r'''# 28 列演示池上的 RMP 与对偶
m = mathopt.Model()
xv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
covers = []
for i in range(1, n+1):
    covers.append(m.add_linear_constraint(
        mathopt.fast_sum([xv[p] for p in range(P) if (masks[p] >> i) & 1]) >= 1.0, name=f"c{i}"))
veh = m.add_linear_constraint(mathopt.fast_sum(xv) <= K, name="veh")
m.minimize(mathopt.fast_sum([costs[p]*xv[p] for p in range(P)]))
res = mathopt.solve(m, mathopt.SolverType.GLOP)
dv = res.dual_values()
pi = [0.0]*(n+1)
for i in range(1, n+1):
    pi[i] = max(0.0, dv[covers[i-1]])
mu = dv[veh]
xvals = {p: res.variable_values()[xv[p]] for p in range(P)}
print("RMP 目标 =", round(res.objective_value(), 6), "| 车辆数对偶 mu =", round(mu, 6))
print("覆盖对偶 pi =", [round(pi[i], 3) for i in range(1, n+1)])
dual_obj = sum(pi[1:]) + K*mu
print("对偶目标 Σπ+Kμ =", round(dual_obj, 6),
      "| 强对偶(原=对偶):", abs(res.objective_value()-dual_obj) < 1e-6)
print()
print("互补松弛（绑定列 rc = 0）:")
for p in range(P):
    if xvals[p] > 1e-6:
        s = 0.0
        mm = masks[p]
        while mm:
            lb = mm & -mm
            s += pi[lb.bit_length()-1]
            mm -= lb
        rc = costs[p] - s - mu
        print(f"  列 {pool[p]}  x={round(xvals[p],3)}  rc={rc:.2e}（应≈0）")
p0 = pool.index((20,))
print()
print(f"非绑定单客户列 (20,): c={round(costs[p0],4)}  rc = c - pi_20 - mu = "
      f"{costs[p0]-pi[20]-mu:.6f}（≥0，对偶约束满足）")
'''),
]

# ================= 03 =================
NB03 = [
 ("MD", r"""# 03 Benders —— 子问题对偶与最优性割的完整推导

## 分解结构

主问题选列 $y_p\in\{0,1\}$（是否允许使用列 $p$），子问题在 $y$ 允许的列上解连续覆盖 LP。
列成本含在子问题中，故主问题目标只取 $\min\theta$（割把成本信息带回主问题）。

## 子问题 SP(y) 与它的对偶

$$\min_{x,s}\ \sum_p c_px_p+M\sum_i s_i \quad\text{s.t.}\quad \sum_p a_{ip}x_p+s_i\ge 1\ (\pi_i\ge 0),\ \sum_p x_p\le K\ (\mu\le 0),\ 0\le x_p\le y_p\ (\sigma_p\le 0),\ s_i\ge 0$$

按 02 的对偶规则逐条推导：

- 覆盖约束（≥）→ $\pi_i\ge 0$；车辆数（≤）→ $\mu\le 0$；上界 $x_p\le y_p$（≤）→ $\sigma_p\le 0$（目标系数 $y_p$）
- 变量 $x_p\ge 0$ → 对偶约束 $\sum_i a_{ip}\pi_i+\mu+\sigma_p\le c_p$
- 变量 $s_i\ge 0$（目标系数 $M$）→ 对偶约束 $\pi_i\le M$

$$\max_{\pi,\mu,\sigma}\ \sum_i\pi_i+K\mu+\sum_p\sigma_p y_p \quad\text{s.t.}\quad \sum_i a_{ip}\pi_i+\mu+\sigma_p\le c_p,\ \pi_i\le M,\ \pi_i\ge 0,\ \mu\le 0,\ \sigma_p\le 0$$

## Benders 最优性割的推导与有效性

设 $V(y)=SP(y)$ 的最优值，对偶可行点 $(\pi^k,\mu^k,\sigma^k)$ 给出对**任意** $y$ 都成立的下界：

$$V(y)\ \ge\ \sum_i\pi_i^k+K\mu^k+\sum_p\sigma_p^k y_p \qquad(\text{弱对偶：对偶可行值 ≤ 原最优值})$$

取 $\lambda_p=-\sigma_p\ge 0$，主问题引入 $\theta$（$V(y)$ 的逐次下界近似），割：

$$\theta+\sum_p\lambda_p^k y_p\ \ge\ \sum_i\pi_i^k+K\mu^k$$

- **有效性**：对偶可行点对任意 $y$ 都给出 $V(y)$ 的下界，割不过切。
- **收敛**：主问题 $\min\theta$ 被割逐次抬高，收敛到 $\min_y V(y)$ = 候选池的 LP 松弛值（=191.813620）。
- 实现细节：MathOpt 对 $x\le y$ 返回非正对偶，$\lambda=-\sigma$；车辆数对偶 $\mu\le 0$ 直接进割常数项 $K\mu$；
  人工变量使 SP 对任意 $y$ 可行（无需可行性割）。
"""),
 ("CODE", r'''# SP(y=全1) = 演示池上的完整 LP；取对偶生成第一条割
M = 10**6
yvec = [1]*P
sp = mathopt.Model()
xv = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
sv = [sp.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"s{i}") for i in range(1, n+1)]
covers = []
for i in range(1, n+1):
    covers.append(sp.add_linear_constraint(
        mathopt.fast_sum([xv[p] for p in range(P) if (masks[p] >> i) & 1]) + sv[i-1] >= 1.0, name=f"c{i}"))
veh = sp.add_linear_constraint(mathopt.fast_sum(xv) <= K, name="veh")
xub = [sp.add_linear_constraint(xv[p] <= float(yvec[p]), name=f"xub{p}") for p in range(P)]
sp.minimize(mathopt.fast_sum([costs[p]*xv[p] for p in range(P)]) + M*mathopt.fast_sum(sv))
res = mathopt.solve(sp, mathopt.SolverType.GLOP)
dv = res.dual_values()
pi = [max(0.0, dv[covers[i-1]]) for i in range(1, n+1)]
mu = dv[veh]
sigma = [dv[xub[p]] for p in range(P)]
alpha = sum(pi) + K*mu
lam = [-s for s in sigma]
V = res.objective_value()
print("SP(全1) 值 V =", round(V, 6), "（= 演示池 LP 松弛值）")
print("对偶: Σπ =", round(sum(pi), 4), "| Kμ =", round(K*mu, 4), "| 非零 λ 列数:", sum(1 for v in lam if v > 1e-9))
print("割:  theta + Σ λ_p y_p >= Σπ + Kμ =", round(alpha, 6))
# 强对偶：割在 y=全1 处取等
cut_at_ones = alpha - sum(lam)
print("割在 y=全1 处: theta >= alpha - Σλ*1 =", round(cut_at_ones, 6), "= V ->", abs(cut_at_ones - V) < 1e-6)
# 有效性验证：对另一个 y'（仅单客户列）解 SP，割仍应 ≤ V(y')
y2 = [1 if len(pool[p]) == 1 else 0 for p in range(P)]
sp2 = mathopt.Model()
xv2 = [sp2.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
sv2 = [sp2.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"s{i}") for i in range(1, n+1)]
for i in range(1, n+1):
    sp2.add_linear_constraint(
        mathopt.fast_sum([xv2[p] for p in range(P) if (masks[p] >> i) & 1]) + sv2[i-1] >= 1.0, name=f"c{i}")
sp2.add_linear_constraint(mathopt.fast_sum(xv2) <= K, name="veh")
for p in range(P):
    sp2.add_linear_constraint(xv2[p] <= float(y2[p]), name=f"xub{p}")
sp2.minimize(mathopt.fast_sum([costs[p]*xv2[p] for p in range(P)]) + M*mathopt.fast_sum(sv2))
res2 = mathopt.solve(sp2, mathopt.SolverType.GLOP)
cut_at_y2 = alpha - sum(lam[p]*y2[p] for p in range(P))
print(f"y'=仅单客户列: V(y') = {round(res2.objective_value(),4)} | 割在该 y' 处 = {round(cut_at_y2,4)} | 割有效(≤V): {cut_at_y2 <= res2.objective_value() + 1e-6}")
'''),
]

# ================= 04 =================
NB04 = [
 ("MD", r"""# 04 拉格朗日松弛 —— 对偶函数与次梯度的完整推导

## 松弛覆盖约束

原问题 $\min\{\sum_p c_px_p:\ \sum_p a_{ip}x_p\ge 1,\ \sum_p x_p\le K,\ x_p\in\{0,1\}\}$。
给覆盖约束配乘子 $\lambda_i\ge 0$（惩罚"未覆盖"），得对偶函数：

$$L(\lambda)=\min_{\Sigma x\le K,\ x\in\{0,1\}}\Big[\sum_p c_px_p+\sum_i\lambda_i\Big(1-\sum_p a_{ip}x_p\Big)\Big]
=\sum_i\lambda_i+\min_{\Sigma x\le K}\sum_p\Big(\underbrace{c_p-\sum_i a_{ip}\lambda_i}_{rc_p}\Big)x_p$$

**子问题（平凡，无需求解器）**：按 $rc_p$ 排序取最负的至多 $K$ 个列（$rc\ge 0$ 不取）。

**对偶界**：任意 $\lambda\ge 0$ 都有 $L(\lambda)\le$ 原最优值（松弛约束只会降低 min 值）；
$\max_{\lambda\ge 0}L(\lambda)$ 是对偶问题（次梯度法求解）。

## 次梯度是怎么来的

**命题**：$g_i=1-\sum_p a_{ip}x_p(\lambda)$（$x(\lambda)$ 为 $L(\lambda)$ 的子问题解）是 $L$ 在 $\lambda$ 处的**次梯度**。

**证明**（对任意 $\lambda'\ge 0$，由子问题最优性 $L(\lambda)$ 取到 $x(\lambda)$，而 $L(\lambda')$ 是其 min）：

$$L(\lambda')=\min_x\Big[\sum_p c_px_p+\sum_i\lambda'_i(1-\sum_p a_{ip}x_p)\Big]\le\sum_p c_px_p(\lambda)+\sum_i\lambda'_i(1-\sum_p a_{ip}x_p(\lambda))$$

$$=\Big[\sum_p c_px_p(\lambda)+\sum_i\lambda_i(1-\sum_p a_{ip}x_p(\lambda))\Big]+\sum_i(\lambda'_i-\lambda_i)(1-\sum_p a_{ip}x_p(\lambda))
=L(\lambda)+\sum_i g_i(\lambda'_i-\lambda_i)$$

即 $L(\lambda')\ge L(\lambda)+\langle g,\lambda'-\lambda\rangle$，$g$ 是凹函数 $L$ 的次梯度（超平面支撑）。
**注意 $g_i=1-$覆盖次数**（布尔标志不是合法次梯度——02/04 家族实测的坑）。

## 为什么对偶 = LP 松弛值（integrality property）

子问题 $\min\{\sum rc_px_p:\Sigma x\le K,\ x_p\in\{0,1\}\}$ 的 LP 松弛（$x_p\in[0,1]$）最优值相同
（都等于取最负的至多 $K$ 个 rc 之和）。子问题具有"整点性质"⇒ 拉格朗日对偶
$\max L(\lambda)$ = 原问题的 LP 松弛值 = **191.813620**（02 已证），对偶间隙为零。
"""),
 ("MD", r"""## 次梯度算法与实现要点

1. 迭代：$\lambda_i\leftarrow\max(0,\ \lambda_i+\alpha g_i)$，$\alpha=\rho\,(UB-L(\lambda))/\|g\|^2$；
   $\rho=2.0$ 起、30 轮无下界改进减半；**步长截断 100、λ 截断 1000**（$\|g\|\to 0$ 时步长发散）。
2. 修复：贪心给未覆盖客户补列（≤K 列）得上界；最后候选池 CP-SAT MIP 修复。
3. 停机：600 轮 / ρ<1e-5 / 墙钟 110s。确定性。
4. 与 02 的联系：取 $\lambda=\pi^*$（LP 最优对偶）时 $rc_p\ge 0$ 对全部列成立 ⇒ 子问题不选任何列 ⇒
   $L(\pi^*)=\sum_i\pi_i^*$ = LP 对偶目标 = LP 值 —— 强对偶点。
"""),
 ("CODE", r'''# 先解演示池 RMP 取对偶 π（作为 λ）
m = mathopt.Model()
xv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
covers = []
for i in range(1, n+1):
    covers.append(m.add_linear_constraint(
        mathopt.fast_sum([xv[p] for p in range(P) if (masks[p] >> i) & 1]) >= 1.0, name=f"c{i}"))
m.add_linear_constraint(mathopt.fast_sum(xv) <= K, name="veh")
m.minimize(mathopt.fast_sum([costs[p]*xv[p] for p in range(P)]))
res = mathopt.solve(m, mathopt.SolverType.GLOP)
dv = res.dual_values()
lam = [0.0]*(n+1)
for i in range(1, n+1):
    lam[i] = max(0.0, dv[covers[i-1]])
lp_val = res.objective_value()
print("LP 值 =", round(lp_val, 6), "| λ = LP 对偶 π")

# L(λ)：排序子问题
rcs = []
for p in range(P):
    s = 0.0
    mm = masks[p]
    while mm:
        lb = mm & -mm
        s += lam[lb.bit_length()-1]
        mm -= lb
    rcs.append((costs[p]-s, p))
rcs.sort()
S = [p for rc, p in rcs[:K] if rc < -1e-7]
L = sum(lam[1:]) + sum(costs[p]-sum(lam[i] for i in pool[p]) for p in S)
print(f"子问题选中 {len(S)} 列 | L(λ) = {round(L, 6)} | 与 LP 值相等: {abs(L-lp_val) < 1e-6}")
print("（λ=最优对偶时 rc_p≥0 对全部列成立，子问题不选列，L=Σλ=LP 对偶目标——强对偶点）")

# 次梯度（覆盖次数）与支撑性质验证（λ'=0）
cnt = [0]*(n+1)
for p in S:
    mm = masks[p]
    while mm:
        lb = mm & -mm
        cnt[lb.bit_length()-1] += 1
        mm -= lb
g = [1.0 - cnt[i] for i in range(1, n+1)]
print("次梯度 g = 1 - 覆盖次数 =", g[:6], "...")
L0 = 0.0   # λ'=0：rc=c>0，子问题不选列 -> L(0)=0
rhs = L + sum(g[i-1]*(0.0 - lam[i]) for i in range(1, n+1))
print(f"支撑性质: L(0)={L0} >= L(λ)+Σg·(0-λ) = {round(rhs,6)} -> {L0 >= rhs - 1e-9}")
'''),
]

# ================= 05 =================
NB05 = [
 ("MD", r"""# 05 LBBD —— 逻辑割（no-good / Hooker）的推导与有效性

## 主问题（CP-SAT 分配）与子问题（单车辆 TSP-TW）

$$\min\theta\ \text{s.t.}\ \sum_k a_{ik}=1,\ a_{ik}\le v_k,\ \sum_i\delta_i a_{ik}\le Q,\ \sum_k v_k=K^*,\ \theta\le\overline{UB}-1$$

子问题对车辆 $k$ 的客户集 $S_k$ 解单车辆 TSP-TW（min 距离），返回可行性与成本 $c_k$。
逻辑 Benders 不需要对偶——**割由子问题的可行性/最优性直接导出**。

## 割一：no-good（不可行集）

**事实**：若 $S_k$ 不可行（无满足时间窗+容量+返回仓库的排序），则任何把 $S_k$ 整组分配给车 $k$ 的
分配都不可行。**割**：$\sum_{i\in S_k}(1-a_{ik})\ge 1$（至少移走一个客户）。

**加强**：贪心收缩到**最小不可行核心**（逐客户尝试移除，仍不可行则移除），核心更小 ⇒ 割更强。
**有效性**：不切掉任何可行分配（只禁止了必然不可行的分配）。

## 割二：Hooker 最优性割（可行集）

设车辆 $k$ 对 $S_k$ 的最优成本为 $c_k$，定义**单客户移除节省**
$\delta_{ik}=c_k-c_k(S_k\setminus\{i\})$（子问题重解精确得到；移除后不可行时取 $\delta_{ik}=c_k$）。
对任意分配 $a'$（车 $k$ 分到 $S'_k$），在"多客户移除节省可加"的标准假设下：

$$c_k(S'_k)\ \ge\ c_k-\sum_{i\in S_k\setminus S'_k}\delta_{ik} \qquad\Longrightarrow\qquad
\theta\ \ge\ \sum_k\Big[c_k-\sum_{i\in S_k}\delta_{ik}(1-a_{ik})\Big]$$

**有效性口径**：单客户移除是精确的（重解验证）；多客户同时移除依赖标准可加性假设（TSP-TW 子问题的
经典 LBBD 口径，02/03/05 notebook 已如实声明，最终结论与 01/02/03/04/07 交叉验证一致）。

## 割三：精确成本界（严格有效）

$$\theta\ \ge\ C\cdot\Big(1-\sum_k\sum_{i\in S_k}(1-a_{ik})\Big)$$

指示器 $\mathbb 1[a=a^*]=1-\sum(1-a)\in\{0,1\}$：完全复现该分配时右端=$C$（精确），否则右端≤0≤θ（平凡成立）。
**无条件有效**，是主问题不可行性证明的"锚"。

## 最优性证明

主问题含 $\theta\le\overline{UB}-1$（必须找到比现任解更好的分配）。若主问题**不可行**，则任何可行分配
都被割覆盖且成本 ≥ UB ⇒ **UB 最优**。这就是 05 第 22 轮"主问题不可行 ⇒ 已证明最优"的数学依据。
"""),
 ("CODE", r'''from lbbd import solve_route, shrink
# no-good 演示：全部 25 客户塞一车 -> 不可行 -> 收缩核心
S_all = tuple(range(1, n+1))
st, route, cost = solve_route(S_all)
print(f"全部客户一车: {st} | 最小不可行核心 = {shrink(S_all)}")
print("no-good 割: Σ_{i∈核心}(1-a_i) >= 1（至少移走一个客户）")

# Hooker 割演示：最优路线集 S1 的成本与移除节省
S1 = ROUTES[0]
st1, r1, ck = solve_route(S1)
print()
print(f"路线 {S1}: 最优成本 c_k = {round(ck, 4)} | 状态 {st1}")
deltas = {}
for i in S1:
    S2 = tuple(j for j in S1 if j != i)
    st3, _, c2 = solve_route(S2)
    dlt = ck - c2 if st3 == "OPTIMAL" else ck
    deltas[i] = round(dlt, 4)
print("单客户移除节省 δ:", deltas)
print("Hooker 割: theta >= c_k - Σ δ_i(1-a_i)（θ 为该车成本的下界近似）")

# 精确界割（严格有效）
print("精确界割: theta >= C·(1 - Σ(1-a_i))（仅完全复现该分配时抬到 C，否则 0）")
'''),
]

# ================= 07 =================
NB07 = [
 ("MD", r"""# 07 Branch-and-Price —— 节点 RMP 的对偶、分支与剪枝的数学

## 节点 RMP 与对偶

每个 B&B 节点解（含车辆数上下界与虚拟列 $y_i$，成本 $M$）：

$$\min\ \sum_p c_px_p+M\sum_i y_i \quad\text{s.t.}\quad \sum_p a_{ip}x_p+y_i\ge 1\ (\pi_i\ge 0),\ \sum_p x_p\le K_{ub}\ (\mu_{ub}\le 0),\ \sum_p x_p\ge K_{lb}\ (\mu_{lb}\ge 0)$$

按 02 的对偶规则：**下界约束（≥）对偶 $\mu_{lb}\ge 0$**，上界约束（≤）对偶 $\mu_{ub}\le 0$，虚拟列给 $\pi_i\le M$：

$$\max\ \sum_i\pi_i+K_{ub}\mu_{ub}+K_{lb}\mu_{lb}\ \ \text{s.t.}\ \ \sum_i a_{ip}\pi_i+\mu_{ub}+\mu_{lb}\le c_p\ (\forall p),\ \pi_i\le M,\ \pi_i\ge 0,\ \mu_{ub}\le 0,\ \mu_{lb}\ge 0$$

**定价 reduced cost**（两车辆对偶合并）：
$$rc_p=c_p-\sum_{i\in p}\pi_i-(\mu_{ub}+\mu_{lb})$$

- 定价子问题：02 的 CP-SAT ESPPRC + **强制弧 $x_{ij}=1$ / 禁用弧 $x_{ij}=0$**（分支约束传播进定价）；
  池扫描同步按弧过滤。
- 虚拟列取正值 ⇒ 节点无真实覆盖 ⇒ 不可行剪枝（$\pi_i=M$ 时所有真实列 rc<0 但池中无符合弧限制的列）。

## 分支规则的数学

- **车辆数分支**：$\Sigma_p x_p=f$ 分数 → 左子 $\Sigma x\le\lfloor f\rfloor$、右子 $\Sigma x\ge\lceil f\rceil$，
  把分数车辆数区间剖开（对应 RMP 的 $K_{ub}/K_{lb}$）。
- **Ryan-Foster 弧分支**：定义弧流量 $f_{ij}=\sum_p x_p\cdot\mathbb 1[(i,j)\in p]$。
  **事实**：$x$ 分数 ⇒ 存在客户弧 $(i,j)$ 使 $f_{ij}\in(0,1)$（若全部 $f_{ij}\in\{0,1\}$，由流量守恒可证各列 $x_p\in\{0,1\}$）。
  取最接近 0.5 的弧分支：左子禁用（$x_{ij}=0$）、右子强制（$x_{ij}=1$），子代弧集不相交 ⇒ 树有限。

## 剪枝正确性（三个引理）

1. **定界剪枝**：节点 LP 值 ≥ 现任解 UB ⇒ 子树不可能出现更优整数解（节点 LP 是其下界）。
2. **不可行剪枝**：虚拟列取正 ⇒ 不存在满足分支弧限制的真实覆盖。
3. **整数节点**：$x_p\in\{0,1\}$ ⇒ LP 解即整数解，更新 UB 后剪枝（节点已解决）。
"""),
 ("CODE", r'''# 节点演示：K_lb=K_ub=3 的节点 RMP（28 列池 + 虚拟列）
M = 10**6
K_lb = K_ub = 3
m = mathopt.Model()
xv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"x{p}") for p in range(P)]
yv = [m.add_variable(lb=0.0, ub=float("inf"), is_integer=False, name=f"y{i}") for i in range(1, n+1)]
covers = []
for i in range(1, n+1):
    covers.append(m.add_linear_constraint(
        mathopt.fast_sum([xv[p] for p in range(P) if (masks[p] >> i) & 1]) + yv[i-1] >= 1.0, name=f"c{i}"))
ubc = m.add_linear_constraint(mathopt.fast_sum(xv) <= K_ub, name="vub")
lbc = m.add_linear_constraint(mathopt.fast_sum(xv) >= K_lb, name="vlb")
m.minimize(mathopt.fast_sum([costs[p]*xv[p] for p in range(P)]) + M*mathopt.fast_sum(yv))
res = mathopt.solve(m, mathopt.SolverType.GLOP)
dv = res.dual_values()
pi = [0.0]*(n+1)
for i in range(1, n+1):
    pi[i] = max(0.0, dv[covers[i-1]])
mu_ub = dv[ubc]
mu_lb = dv[lbc]
xvals = {p: res.variable_values()[xv[p]] for p in range(P)}
y_used = sum(res.variable_values()[yv[i-1]] for i in range(1, n+1))
print("节点 LP =", round(res.objective_value(), 6), "| 虚拟列使用 =", round(y_used, 6))
print("mu_ub =", round(mu_ub, 4), "| mu_lb =", round(mu_lb, 4),
      "| rc = c - Σπ - μ_ub - μ_lb（正列验证）:")
for p in range(P):
    if xvals[p] > 1e-6:
        s = 0.0
        mm = masks[p]
        while mm:
            lb = mm & -mm
            s += pi[lb.bit_length()-1]
            mm -= lb
        rc = costs[p] - s - mu_ub - mu_lb
        print(f"  列 {pool[p]}  x={round(xvals[p],3)}  rc={rc:.2e}")
# 弧流量与整数性
flows = {}
for p, val in xvals.items():
    if val <= 1e-6:
        continue
    seq = [0] + list(pool[p]) + [0]
    for a in range(len(seq)-1):
        flows[(seq[a], seq[a+1])] = flows.get((seq[a], seq[a+1]), 0.0) + val
frac_arcs = [(a, f) for a, f in flows.items() if a[0] >= 1 and a[1] >= 1 and 1e-6 < f < 1-1e-6]
integral = all(abs(v-round(v)) < 1e-6 for v in xvals.values()) and y_used <= 1e-6
print("分数客户弧数:", len(frac_arcs), "| LP 整数解:", integral,
      "-> 根节点即证明最优（无需分支）" if integral else "-> 需分支（车辆数/弧）")
'''),
]

NBS = {
  "01_direct_数学与最优性证书.ipynb": NB01,
  "02_列生成_LP对偶与reduced_cost.ipynb": NB02,
  "03_Benders_子问题对偶与割推导.ipynb": NB03,
  "04_拉格朗日_对偶函数与次梯度.ipynb": NB04,
  "05_LBBD_逻辑割推导.ipynb": NB05,
  "07_Branch_and_Price_节点对偶与分支.ipynb": NB07,
}

for fname, cells in NBS.items():
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    seen_code = False
    for kind, content in cells:
        if kind == "MD":
            nb.cells.append(nbf.v4.new_markdown_cell(content))
        else:
            code = (PREAMBLE + "\n\n" + content) if not seen_code else content
            seen_code = True
            nb.cells.append(nbf.v4.new_code_cell(code))
    nbf.write(nb, os.path.join(BASE, fname))
    print("written", fname)

idx = r"""# theory/ —— 数学建模与对偶推导详解系列（全部已执行）

按方法顺序，每个 notebook 给出：完整数学公式（LaTeX）+ 对偶/割**逐步推导** + 数值验证（小池演示，秒级运行）。

| notebook | 内容 | 对偶/割来源 |
|---|---|---|
| 01_direct_数学与最优性证书.ipynb | 3-index CP-SAT 完整模型 | 无对偶；最优性证书 = 目标=best bound |
| 02_列生成_LP对偶与reduced_cost.ipynb | RMP 模型、LP 对偶规则表、乘子法推导、rc 三命题 | 对偶约束的松弛量 rc；定价=找最违反的对偶约束 |
| 03_Benders_子问题对偶与割推导.ipynb | SP(y) 模型与对偶、割推导与有效性、数值验证 | 对偶可行点 → θ+Σλy ≥ Σπ+Kμ（弱对偶） |
| 04_拉格朗日_对偶函数与次梯度.ipynb | L(λ) 推导、平凡子问题、次梯度不等式证明、integrality property | 乘子法 + 子问题最优性不等式；g=1−覆盖次数 |
| 05_LBBD_逻辑割推导.ipynb | no-good/Hooker/精确界三种割推导与有效性 | 子问题可行性/最优性直接导出（无需对偶） |
| 07_Branch_and_Price_节点对偶与分支.ipynb | 节点 RMP 对偶（μ_ub/μ_lb）、弧流量分支、剪枝引理 | 车辆数上下界对偶进定价 rc |

配套：family 内 01–07 全部 notebook（完整求解）；02_column_generation_建模详解.md（主问题/对偶/子问题总述）；
08_主问题与子问题建模总览.md（七项一览表）。
"""
with open(os.path.join(BASE, "README.md"), "w", encoding="utf-8") as f:
    f.write(idx)
print("written theory/README.md")
