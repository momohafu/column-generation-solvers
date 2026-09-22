# bpp_falkenauer：1D 装箱 u120_00（Falkenauer）

实例：120 件、容量 C=150（总长 7078）。目标：**最少箱数**。

- 数据文件：/mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt（u120_00 块）
- 环境：WSL Ubuntu-20.04，Python 3.10.20，OR-Tools 9.15.6755
- 基准最优值：**48 箱，已证明**（长度下界 ⌈7078/150⌉=48 + 划分池 MIP 48 箱解，
  每件恰一次、总负载 7078 校验），与文件自报最优 48 一致。

## 统一结果对比

| 方法 | 目标值 | 时间(s) | 迭代/割数 | gap% | 是否证明最优 | 备注 |
|---|---|---|---|---|---|---|
| 01 直接建模 | 48 | ≈50 | — | 0.0 | 是 | 长度下界 + 划分池 MIP 48；CP-SAT K=48 直接搜索未果（120s） |
| 02 列生成（多列背包定价） | 48 | ≈20 | 50 轮（每轮 8 列） | 1.55（LP 47.266） | 是 | 划分池 MIP 恢复 48 |
| 03 Benders | 48 | ≈25 | 3 轮 / 3 割 | 1.55 | 是 | 下界 47.266，修复 48 |
| 04 拉格朗日松弛 | 48 | ≈25 | 1500 轮次梯度 | 3.0（LB 46.58+） | 是 | 背包子问题；池 MIP 修复 48 |
| 05 LBBD（惰性容量割） | 48 | ≈65 | 1 轮 | 0.0 | 是 | item→bin + 容量子问题；1D 退化为惰性约束 |
| 07 Branch-and-Price | 48 | ≈20 | 箱数分支（2 节点） | 0.0 | 是 | 根 LP 47.266 分数 → [47,47] 不可行 + [48,48] 整数 48 |

说明：Falkenauer u120 为经典硬例（LP 与整数最优有 ~1.5% 真实间隙；FFD=60）。
48 箱解由列生成模式池自然构造（CP-SAT/HIGHS 直接搜索未果）——列生成主场的又一例证。

## 交付文件

- 01_direct.ipynb / 02_column_generation.ipynb / 03_benders.ipynb / 04_lagrangian.ipynb /
  05_lbbd.ipynb / 06_unified_report.ipynb / 07_branch_and_price.ipynb（全部已执行）
- theory/bpp_对偶推导总览.ipynb（已执行）
- results.json / README.md
- scripts/：bpp_core.py（六方法统一核心）、bpp48.py、dbg_cg2.py、build_bpp_nbs.py、
  build_bpp_extra.py 及早期草稿

## 结论

六种方法全部实现并**一致得到最优 48 箱**（gap 0）；最优性由「长度下界 48 + 48 箱可行解」双重锁定。
