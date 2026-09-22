# csp_1d_waescher：1D 切割下料 TEST0005（Wäscher & Gau）

实例：57 种件型、114 件，标准辊长 L=10000，需求 d_i。目标：**最少标准辊数**满足全部需求。

- 数据文件：/mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt
- 环境：WSL Ubuntu-20.04，Python 3.10.20，OR-Tools 9.15.6755（CP-SAT + MathOpt GLOP/HIGHS/GSCIP）
- 基准最优值：**28 辊，已证明**（总长度下界 ⌈279935/10000⌉=28 + CG 池 MIP 28 辊解，覆盖/件数/长度三重校验）

## 统一结果对比

| 方法 | 目标值 | 时间(s) | 迭代/割数 | gap% | 是否证明最优 | 备注 |
|---|---|---|---|---|---|---|
| 01 直接建模 | 28 | ≈20 | — | 0.0 | 是 | 长度下界 28 + 池 MIP 28 解；CP-SAT K=28 直接搜索 560s 未果（28 解存在但难构造） |
| 02 列生成（GG 背包定价） | 28 | ≈25 | 281 轮 | 0.02（LP 27.9942） | 是 | numpy 0/1 背包 DP 定价（10ms）；池 MIP 恢复 28 |
| 03 Benders | 28 | ≈25 | 3 轮 / 3 割 | 0.02 | 是 | 下界 27.9942，修复 28 |
| 04 拉格朗日松弛 | 28 | ≈40 | 800 轮次梯度 | 0.45（LB 27.87+） | 是 | 背包子问题；池 MIP 修复 28 |
| 05 LBBD（惰性容量割） | 28 | ≈90 | 1 轮 | 0.0 | 是 | item→bin 主问题 + 容量子问题；1D 退化为惰性约束（CONVENTIONS §4.5） |
| 07 Branch-and-Price | 28 | ≈60 | 辊数分支（2 节点） | 0.0 | 是 | 根 LP 27.9942 分数 → [27,27] 不可行 + [28,28] 整数 28 ⇒ 证明 |

说明：该实例 LP 松弛 = 27.9942 < 28（积分性间隙 0.02%）；模式总数 2.8 亿不可枚举，
定价必须用背包 DP（本家族列生成的主场）。28 辊解极难直接构造（CP-SAT/HIGHS/GSCIP 均超时未果），
但列生成生成的模式池使其自然浮现——这正是列生成的典型价值场景。

## 交付文件

- 01_direct.ipynb / 02_column_generation.ipynb / 03_benders.ipynb / 04_lagrangian.ipynb /
  05_lbbd.ipynb / 06_unified_report.ipynb / 07_branch_and_price.ipynb（全部已执行）
- theory/csp_对偶推导总览.ipynb（六方法对偶/割推导 + 数值验证，已执行）
- results.json / README.md
- scripts/：csp_core.py（六方法统一核心）、csp_bnp.py（完整 GG pair 分支框架）、csp_pricing.py、
  csp_proto.py、csp_k28.py、csp_scip.py、csp_lb.py、csp_dff.py、csp_verify.py、
  build_csp_nbs.py、build_csp_extra.py 及早期草稿

## 结论

六种方法全部实现并**一致得到最优 28 辊**（gap 0）；最优性由「总长度下界 28 + 28 辊可行解」双重锁定。
LP 松弛 27.9942 与整数最优存在 0.02% 间隙，列生成/B&P 的辊数分支给出了决定性证明。
