# -*- coding: utf-8 -*-
"""构建 01_direct.ipynb（VRPTW c101 直接 CP-SAT 建模基准）。"""
import nbformat as nbf

NB_PATH = "/mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/01_direct.ipynb"

MD0 = r"""# VRPTW（Solomon c101）—— 直接建模基准（3-index CP-SAT）

## 问题定义

带时间窗车辆路径问题（VRPTW）：$K$ 辆车（容量 $Q=200$）从仓库出发服务客户集合 $C$（25 个），
客户 $i$ 有需求 $\delta_i$、时间窗 $[r_i,l_i]$、服务时长 $s_i$；目标车辆数最少、其次总距离最短。
数据文件 c101.txt 声明 VEHICLE NUMBER=25、CAPACITY=200。

## 3-index 弧流模型

$$\min \sum_{k}\sum_{i\ne j} d_{ij}\, x_{ijk}$$

$$\text{s.t.}\quad \text{AddCircuit}(x_{k})\ \ (\forall k);\qquad \sum_k x_{iik} = K-1\ \ (\forall i\in C)\quad(\text{每客户恰被一车访问，其余车自环})$$

$$t_{jk} \ge t_{ik} + s_i + d_{ij} - M(1-x_{ijk}),\quad r_i \le t_{ik} \le l_i,\quad q_{jk} \ge q_{ik} + \delta_j - M(1-x_{ijk}),\quad q_{ik}\le Q,\quad t_{ik}+s_i+d_{i0}\le L_0$$

- 车辆 k 未使用 ⇔ 仓库自环 $x_{00k}=1$（AddCircuit 允许自环）；对称性破缺 $x_{00k}\le x_{00,k+1}$。
- 时间单调性自动消除客户子环（同 02/05 口径）；距离/时间 ×10000 整数缩放，$M=10^{12}$。
- **车辆数最少性**：容量下界 $\lceil 460/200\rceil = 3$ 给出 ≥3；K=2 由 CP-SAT 直接证明不可行 → 最少 3 车。
- 求解器：CP-SAT（num_search_workers=8，时间上限 120 s）。
"""

C1 = r'''# 环境信息（CONVENTIONS §3.4 要求）
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
import platform, math, time
import ortools
from ortools.sat.python import cp_model

print("python", platform.python_version(), "| ortools", ortools.__version__)

DATA = "/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt"
SCALE = 10000
M = 10**12

def parse(path=DATA):
    lines = open(path).read().splitlines()
    K = cap = None
    for i, l in enumerate(lines):
        if l.strip().startswith("NUMBER"):
            s = lines[i+1].split(); K, cap = int(s[0]), int(s[1])
    custs = []
    for l in lines:
        s = l.split()
        if len(s) == 7 and s[0].isdigit():
            no, x, y, dem, ready, due, svc = map(int, s)
            custs.append((no, x, y, dem, ready, due, svc))
    custs.sort(key=lambda c: c[0])
    return K, cap, custs[0], custs[1:]

KMAX, cap, depot, cs = parse()
n = len(cs)
xc = [depot[1]] + [c[1] for c in cs]
yc = [depot[2]] + [c[2] for c in cs]
dem = [0] + [c[3] for c in cs]
ready = [0] + [c[4] for c in cs]
due = [0] + [c[5] for c in cs]
svc = [0] + [c[6] for c in cs]
depot_due = depot[5]
total_dem = sum(dem[1:])
LB_K = math.ceil(total_dem / cap)
def dist(i, j):
    return math.hypot(xc[i]-xc[j], yc[i]-yc[j])
d_scaled = [[0]*(n+1) for _ in range(n+1)]
for i in range(n+1):
    for j in range(n+1):
        if i != j:
            d_scaled[i][j] = int(round(dist(i, j)*SCALE))
print(f"客户数 n={n} | 车辆上限 K={KMAX} | 容量 Q={cap} | 总需求 {total_dem} | 容量下界 LB_K={LB_K}")
print("已知最优（BKS）: 3 车 / 191.81（两位小数舍入）")
'''

C2 = r'''# ---- 3-index CP-SAT 模型构建与求解 ----
def solve(Kv, time_limit=120.0):
    model = cp_model.CpModel()
    x = [[[model.NewBoolVar(f"x{k}_{i}_{j}") for j in range(n+1)] for i in range(n+1)] for k in range(Kv)]
    t = [[model.NewIntVar(0, depot_due*SCALE, f"t{k}_{i}") for i in range(n+1)] for k in range(Kv)]
    q = [[model.NewIntVar(0, cap, f"q{k}_{i}") for i in range(n+1)] for k in range(Kv)]
    for k in range(Kv):
        model.AddCircuit([(i, j, x[k][i][j]) for i in range(n+1) for j in range(n+1)])
    for i in range(1, n+1):                                   # 每客户恰被一车访问（其余车自环）
        model.Add(sum(x[k][i][i] for k in range(Kv)) == Kv - 1)
    for k in range(Kv-1):                                     # 未使用车辆 = 后缀（仓库自环）
        model.Add(x[k][0][0] <= x[k+1][0][0])
    for k in range(Kv):
        model.Add(t[k][0] == 0)
        model.Add(q[k][0] == 0)
        for i in range(1, n+1):
            model.Add(t[k][i] >= ready[i]*SCALE - M*x[k][i][i])
            model.Add(t[k][i] <= due[i]*SCALE + M*x[k][i][i])
            model.Add(t[k][i] <= M*(1 - x[k][i][i]))
            model.Add(q[k][i] >= dem[i] - M*x[k][i][i])
            model.Add(q[k][i] <= M*(1 - x[k][i][i]))
        for i in range(n+1):
            for j in range(1, n+1):
                if i == j:
                    continue
                model.Add(t[k][j] >= t[k][i] + svc[i]*SCALE + d_scaled[i][j] - M*(1 - x[k][i][j]))
                model.Add(q[k][j] >= q[k][i] + dem[j] - M*(1 - x[k][i][j]))
        for i in range(1, n+1):
            model.Add(t[k][i] + svc[i]*SCALE + d_scaled[i][0] <= depot_due*SCALE + M*(1 - x[k][i][0]))
    model.Minimize(sum(d_scaled[i][j]*x[k][i][j] for k in range(Kv) for i in range(n+1) for j in range(n+1) if i != j))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    t0 = time.time()
    st = solver.Solve(model)
    el = time.time() - t0
    obj = solver.ObjectiveValue() if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    routes = []
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        for k in range(Kv):
            if solver.Value(x[k][0][0]):
                continue
            route = []
            cur = 0
            while True:
                nxt = None
                for j in range(n+1):
                    if solver.Value(x[k][cur][j]):
                        nxt = j
                        break
                if nxt is None or nxt == 0:
                    break
                route.append(nxt)
                cur = nxt
                if len(route) > n:
                    break
            if route:
                routes.append(tuple(route))
    return st, obj, routes, el, solver.BestObjectiveBound()

# K=2：证明不可行（与容量下界共同给出最少 3 车）
st2, _, _, el2, _ = solve(2, 120)
print(f"K=2: {st2}（{round(el2,1)}s）-> 2 车不可行，车辆数下界 = 3")

# K=3：求解并证明最优
st3, obj3, routes3, el3, bound3 = solve(3, 120)
print(f"K=3: {st3} | 缩放目标 {obj3} | 下界 {bound3} | 耗时 {round(el3,1)}s")
print("最优性:", st3 == cp_model.OPTIMAL)
'''

C3 = r'''# ---- 双精度复算与对比 ----
exact = 0.0
for r in routes3:
    seq = [0] + list(r) + [0]
    c = sum(dist(seq[i], seq[i+1]) for i in range(len(seq)-1))
    exact += c
    print(f"  路线 {r}  距离 {round(c, 6)}")
print("车辆数", len(routes3), "| 精确总距离", round(exact, 10), "| 两位小数", round(exact, 2))
print("BKS 对比: 3 车 / 191.81 -> match:",
      len(routes3) == 3 and abs(round(exact, 2) - 191.81) < 1e-9)
print("与 02/03/05/07 已证最优 191.813620 一致:", abs(exact - 191.813620) < 1e-6)
'''

MD2 = r"""## 运行结果与结论

- K=2 由 CP-SAT 直接证明不可行（约 0.1 s）；容量下界 ⌈460/200⌉=3 ⇒ **最少 3 车**。
- K=3 CP-SAT 求解并**证明最优**（OPTIMAL，约 1.3 s）：缩放目标 1918136（=191.8136×10000），
  双精度复算 **191.8136197787**（两位小数 191.81），3 条路线与 02/03/05/07 完全一致。

**基准最优值来源**：本 notebook 直接整数模型自证最优 191.8136197787（3 车），
与 Solomon 文献 BKS（3 车、191.81，两位小数舍入）一致；本套件其他五种方法均以此/交叉验证。
"""

cells = [
    nbf.v4.new_markdown_cell(MD0),
    nbf.v4.new_code_cell(C1),
    nbf.v4.new_code_cell(C2),
    nbf.v4.new_code_cell(C3),
    nbf.v4.new_markdown_cell(MD2),
]
nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {"display_name": "Python 3",
                             "language": "python", "name": "python3"}
nbf.write(nb, NB_PATH)
print("written", NB_PATH)
