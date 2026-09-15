# 列生成与分解算法求解套件 —— 总览

> 本套件对下载的经典用例集合中的 6 个问题家族，各取 1 个代表性实例，逐一用
> **直接建模 → 列生成 → Benders 分解 → 拉格朗日松弛 → LBBD** 五种方法求解，
> 并给出完整建模与求解文档。

## 运行环境

- WSL：Ubuntu-20.04
- Conda 环境：jupyter-env（Python 3.10.20）
- 求解器栈：OR-Tools 9.15.6755
  - CP-SAT（约束规划，直接建模 / LBBD 主问题 / 定价子问题）
  - MathOpt + HiGHS（MIP/LP）、GLOP（LP 主问题）、GSCIP（备选 MIP）

## 家族与实例

| # | 家族 | 实例 | 问题 | 目标 | 基准最优值 |
|---|---|---|---|---|---|
| 1 | csp_1d_waescher | Waescher_TEST0005（m=57, L=10000） | 1D 切割下料 | min 辊数 | 待填入 |
| 2 | bpp_falkenauer | u120_00（n=120, C=150） | 1D 装箱 | min 箱数 | 48 |
| 3 | scp_beasley | scp41（200×1000） | 集合覆盖 | min 成本 | 429 |
| 4 | crew_scheduling | csp50（50 任务） | 机组排班 | min 成本 | 27 crew / 3139（已证明） |
| 5 | vrptw_solomon25 | c101（25 客户） | VRPTW | min 车数→距离 | 3 车 / 191.813620（=191.81，已证明） |
| 6 | cutting_2d_cgcut | cgcut1（7 件，15×10） | 2D guillotine 切割 | max 价值 | 244 |

## 结果汇总（全部完成后自动填入）

| 家族 | 直接建模 | 列生成 | Benders | 拉格朗日 | LBBD | 基准 | 全部最优? |
|---|---|---|---|---|---|---|---|
| csp_1d_waescher | | | | | | | |
| bpp_falkenauer | | | | | | | |
| scp_beasley | | | | | | | |
| crew_scheduling | 3139 / 0.05s / 0.0% | 3139 / 0.5s / 0.0% | 3139 / 0.5s / 0.0% | 3139 / 15s / 0.0% | 3139 / 0.2s / 0.0% | 27 crew / 3139 | 是（另含 B&P 0.05s 证明） |
| vrptw_solomon25 | 191.81 / 0.9s / 0.0% | 191.81 / 10s / 0.0% | 191.81 / 11s / 0.0% | 191.81 / 20s / 0.0% | 191.81 / 3.6s / 0.0% | 3车 / 191.81 | 是（另含 B&P 4.5s 证明） |
| cutting_2d_cgcut | | | | | | | |

（格式：目标值 / 时间s / gap%，如 "48 / 3.2s / 0.0%"；未证明最优标记 †）

## 目录结构

```
column-generation-solvers/
├── CONVENTIONS.md          # 共享规范（环境/目录/notebook 要求/方法范式）
├── 00_overview.md          # 本文件
├── results_summary.csv     # 汇总表（由 scripts/aggregate_results.py 生成）
├── scripts/
│   └── aggregate_results.py
├── csp_1d_waescher/        # 每族：5 个已执行 notebook + README.md + results.json
├── bpp_falkenauer/
├── scp_beasley/
├── crew_scheduling/
├── vrptw_solomon25/
└── cutting_2d_cgcut/
```

## 方法说明

- **直接建模**：完整整数模型（HIGHS MIP 或 CP-SAT），给出证明最优基准。
- **列生成**：受限主问题（GLOP LP）+ 定价子问题（背包 DP / 池扫描 / 路径生成），再恢复整数解。
- **Benders 分解**：模式/路径选择（整数主问题）+ 连续 LP 子问题，对偶生成最优性割。
- **拉格朗日松弛**：松弛耦合约束（需求/覆盖），次梯度法求对偶界，贪心修复可行解。
- **LBBD**：CP-SAT 主问题做分配/组合决策 + CP-SAT 可行性子问题，不可行回传逻辑割。

## 每个家族文件夹内容

- `01_direct.ipynb` ... `05_lbbd.ipynb`：问题定义、LaTeX 模型、原理、实现、真实运行结果、结论（已执行、含输出）。
- `README.md`：家族说明 + 五方法结果对比表。
- `results.json`：机器可读结果。
