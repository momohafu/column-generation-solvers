# vrptw_solomon25：VRPTW c101（25 客户）—— 统一报告

实例：Solomon **c101**（25 客户，K=25，Q=200，仓库时间窗 [0,1236]，服务时长 90）。目标：车辆数最少，其次总距离最短。

- 数据文件：/mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt
- 环境：WSL Ubuntu-20.04，Python 3.10.20，OR-Tools 9.15.6755（CP-SAT + MathOpt 的 GLOP/HIGHS）
- 基准最优值：**3 车 / 191.813620（两位小数 191.81）**，由本家族 02/03/05 三种分解方法**独立证明**，
  与 Solomon 文献 BKS（3 车、191.81，两位小数舍入）一致。

## 统一结果对比

| 方法 | 目标值 | 时间(s) | 迭代/割数 | gap% | 是否证明最优 | 备注 |
|---|---|---|---|---|---|---|
| 01 直接建模（3-index CP-SAT） | 191.813620 | ≈0.9（K=2 不可行 0.1 + K=3 最优 0.8） | — | 0.0 | 是 | K=2 不可行 + 容量下界 ⇒ 最少 3 车；K=3 证明最优，本套件基准 |
| 02 列生成（CP-SAT 定价） | 191.813620 | ≈10（枚举2.8+CG6.2+修复1.1） | 8 轮 | 0.0 | 是 | GLOP 覆盖 LP 主问题 + CP-SAT 定价 ESPPRC；LP 下界=整数目标；CP-SAT 列与池最负列每轮一致 |
| 03 Benders（对偶割） | 191.813620 | ≈11（枚举2.0+列池6.6+迭代1.2+修复1.1） | 3 轮 / 3 割 | 0.0 | 是 | 二元选列主问题 min θ + 连续覆盖 LP 子问题；下界=整数修复 |
| 04 拉格朗日松弛 | 191.813620 | ≈20（池10+次梯度9+修复1） | 600 轮次梯度 | 0.0 | 是 | 松弛覆盖约束；对偶下界=191.813620；贪心修复与 MIP 修复均达最优 |
| 05 LBBD（逻辑割） | 191.813620 | ≈3（扫描0.3+迭代1.3~3.0） | 10~22 轮 / no-good 15~39 + Hooker 2 | 0.0 | 是 | 分配主问题 + 单车辆 TSP-TW 子问题；θ≤UB−1 主问题不可行证明最优 |
| 07 Branch-and-Price | 191.813620 | ≈4.5（根节点 CG 5 轮） | 1 节点 / 0 分支 | 0.0 | 是 | 根节点 LP 即整数（3 车）；车辆数/弧分支规则经受控实验验证 |
说明：

- 时间(s) 为 notebook 实际墙钟（主循环 + 枚举/列池构建 + 整数修复）；不同运行间略有波动，以各 notebook 输出为准。
- gap% 相对基准最优 191.813620 计算，三个已完成方法均为 0.0%。
- LBBD 迭代/割数随 CP-SAT 多线程调度有波动（no-good 数 15~39），最终结果与证明不变。
- results.json 与 results_summary.csv（顶层）已生成。

## 最优解（三方法输出完全一致）

| 车 | 路线 | 距离 |
|---|---|---|
| 1 | 0→13→17→18→19→15→16→14→12→0 | 95.884709 |
| 2 | 0→20→24→25→23→22→21→0 | 36.440680 |
| 3 | 0→5→3→7→8→10→11→9→6→4→2→1→0 | 59.488231 |
| 合计 | 3 车 | **191.813620**（=191.81 两位小数） |

**最优性交叉验证**：三种证明机制彼此独立——

1. 列生成：收敛时对偶可行（全池无负 reduced cost 列）⇒ LP 下界 191.813620 = 整数目标；
2. Benders：对偶最优性割收敛，下界 = 整数修复 191.813620；
3. LBBD：no-good/Hooker 逻辑割 + θ≤UB−1 使主问题不可行 ⇒ 不存在更优分配。

## 交付文件

- 01_direct.ipynb（直接建模基准，3-index CP-SAT 证明最优）
- 02_column_generation.ipynb（列生成，CP-SAT 定价子问题）
- 03_benders.ipynb（Benders 分解，选列主问题 + 覆盖 LP 子问题）
- 04_lagrangian.ipynb（拉格朗日松弛，次梯度 + 修复）
- 05_lbbd.ipynb（LBBD，分配主问题 + TSP-TW 子问题 + 逻辑割）
- 07_branch_and_price.ipynb（Branch-and-Price，节点 CG + 车辆数/弧分支，已执行）
- 06_unified_report.ipynb（三方法统一报告，已执行）
- 02_column_generation_建模详解.md（主问题/对偶/定价子问题完整推导）
- 00_进度说明.md（家族进度与状态）
- scripts/：cg_cpsat.py、benders.py、lbbd.py、branch_and_price.py（独立可运行，含 __main__ 保护）、
  verify_consistency.py（三方法一致性校验）、probe_dual.py、probe_arc_branch.py、
  build_nb_02/03/05/06/07.py（notebook 构建器）及早期草稿

## 数学建模与对偶推导详解（theory/ 系列，全部已执行）

| notebook | 内容 |
|---|---|
| 01_direct_数学与最优性证书.ipynb | 3-index 完整模型；无对偶的最优性证书（目标=best bound） |
| 02_列生成_LP对偶与reduced_cost.ipynb | LP 对偶规则表 + 乘子法推导对偶；rc=对偶松弛量三命题；数值验证（强对偶/互补松弛） |
| 03_Benders_子问题对偶与割推导.ipynb | SP(y) 对偶逐步推导；θ+Σλy≥Σπ+Kμ 割的弱对偶有效性；数值验证 |
| 04_拉格朗日_对偶函数与次梯度.ipynb | L(λ) 推导；次梯度支撑不等式证明；integrality property；数值验证 |
| 05_LBBD_逻辑割推导.ipynb | no-good/Hooker/精确界三种逻辑割推导与有效性；数值验证 |
| 07_Branch_and_Price_节点对偶与分支.ipynb | 节点 RMP 对偶（μ_ub/μ_lb）；弧流量分支；三类剪枝引理；数值验证 |

（配套总览：08_主问题与子问题建模总览.md）

## 执行与验证方式

每个 notebook 均用 nbconvert 实际执行并保留输出：

    MSYS_NO_PATHCONV=1 wsl.exe -e bash -lc '/home/zkjqw/miniconda3/envs/jupyter-env/bin/jupyter nbconvert --to notebook --execute --inplace /mnt/d/exactTest/column-generation-solvers/vrptw_solomon25/<nb>.ipynb'

## 结论

- **六种方法全部独立证明最优**（gap 均为 0.0%），目标 191.813620（=191.81）
- **Branch-and-Price** 直接求解该问题：根节点 LP 即整数（3 车 191.813620），0 分支证明最优；
  车辆数分支与 Ryan-Foster 弧分支两种机制经受控实验验证。B&P 是列生成的直接扩展。
- **列生成**是 VRPTW 最自然的分解：定价子问题即 ESPPRC，可直接扩展为 100 节点 branch-and-price。
- **LBBD** 的「分配 + 排序可行性」分层与 VRPTW 结构天然契合，逻辑割无需对偶信息，实现直观。
- **Benders** 因覆盖模型缺少天然连续 recourse，属可运行变体（第 1 轮 y=全1 时子问题即整个 LP），
  定位为展示对偶割机制。
- 六种方法（01/02/03/04/05/07）均在约 20 秒内独立证明最优 **191.813620**（=191.81），最优性被六重确认。
