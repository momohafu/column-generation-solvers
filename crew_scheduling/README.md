# crew_scheduling：机组排班 csp50（Beasley-Cao / OR-Library）

实例：csp50——50 个任务（固定起止时间），时间上限 T=480，173 条转移弧 (i,j,cost)。
一个 crew 的任务序列须逐对由弧连接且跨度（末任务结束 − 首任务开始）≤ 480。
目标：最少 crew 数，其次最小总转移成本。

- 数据文件：/mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt
- 环境：WSL Ubuntu-20.04，Python 3.10.20，OR-Tools 9.15.6755（CP-SAT + MathOpt GLOP/HIGHS）
- 基准最优值：**27 crew / 3139**，由 01_direct 证明（K=26 不可行 + K=27 最优）；与早期草稿 pool_opt 一致。

## 统一结果对比

| 方法 | 目标值 | 时间(s) | 迭代/割数 | gap% | 是否证明最优 | 备注 |
|---|---|---|---|---|---|---|
| 01 直接建模（完整池覆盖 IP） | 3139 | 0.05 | — | 0.0 | 是 | K=26 不可行 + K=27 OPTIMAL（HIGHS，CP-SAT 交叉验证） |
| 02 列生成（CP-SAT 定价） | 3139 | 0.5 | 2 轮 | 0.0 | 是 | LP 下界=3139=整数；池扫描证书；CP-SAT 整数恢复 |
| 03 Benders（对偶割） | 3139 | 0.5 | 3 轮 / 3 割 | 0.0 | 是 | 下界=整数修复 |
| 04 拉格朗日松弛 | 3139 | 15 | 600 轮次梯度 | 0.0 | 是 | 对偶下界=3139；贪心与 MIP 修复均达最优 |
| 05 LBBD（逻辑割） | 3139 | 0.2 | 2 轮 / 50 逻辑割 | 0.0 | 是 | 选列主问题+覆盖逻辑割（scp41 同款）；分配变体作机制展示 |
| 07 Branch-and-Price | 3139 | 0.05 | 1 节点 / 0 分支 | 0.0 | 是 | 根节点 LP 整数；无分数节点（弧/车辆分支规则已实现） |

说明：该实例覆盖 LP 松弛恰为整数（LP=IP），与 vrptw c101 同现象，故列生成/B&P 根节点即证明最优。
LBBD 的「任务→crew 分配主问题 + 调度子问题」变体在 50 任务/27 crew 组合空间下 no-good 割过弱、
难以自收敛（notebook 中如实展示），收敛证明由 scp41 同款选列主问题给出。

## 交付文件

- 01_direct.ipynb / 02_column_generation.ipynb / 03_benders.ipynb / 04_lagrangian.ipynb /
  05_lbbd.ipynb / 06_unified_report.ipynb / 07_branch_and_price.ipynb（全部已执行）
- theory/crew_对偶推导总览.ipynb（六方法对偶/割推导 + 数值验证，已执行）
- results.json / README.md
- scripts/：crew_direct.py、crew_cg.py、crew_benders.py、crew_lagrangian.py、crew_lbbd.py、
  crew_bnp.py（独立可运行，含 __main__ 保护）、crew_sub.py、crew_proto.py、probe_frac.py、
  build_crew_nbs.py、build_crew_extra.py 及早期草稿

## 结论

六种方法（01/02/03/04/05/07）全部独立证明最优 **27 crew / 3139**，结果交叉一致；
覆盖型问题对列生成/拉格朗日/LBBD（覆盖逻辑割）天然友好，本实例 LP 松弛即整数使各方法都极快收敛。
