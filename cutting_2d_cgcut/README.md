# cutting_2d_cgcut：2D guillotine 切割 cgcut1（Christofides-Whitlock）

实例：单张板材 15×10，7 种件型（8×4×2=66、3×7×1=35、8×2×3=24、3×4×5=17、3×3×2=11、3×2×2=8、2×1×1=2）。
切割必须为 guillotine（贯通横/竖切）；**固定朝向**（不允许旋转，文献约定）。目标：最大化切出总价值。

- 数据文件：/mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt
- 环境：WSL Ubuntu-20.04，Python 3.10.20，OR-Tools 9.15.6755
- 基准最优值：**244，已证明**（候选全枚举 2119 → 精确 guillotine 检验 → 完整池 1748 模式 → 池 IP 244），
  与文献 cgcut1 最优 244 一致。注意：允许旋转时另有 260（本家族固定朝向与文献一致，已文档化）。

## 统一结果对比

| 方法 | 目标值 | 时间(s) | 迭代/割数 | gap% | 是否证明最优 | 备注 |
|---|---|---|---|---|---|---|
| 01 直接建模（完整池 IP） | 244 | ≈55（首次建池） | — | 0.0 | 是 | 池完整 ⇒ 精确；模式 (2,1,1,2,1,1,0) |
| 02 列生成（池定价） | 244 | ≈1（池缓存后） | 1 轮 | 0.0 | 是 | LP(conv)=244=IP |
| 03 Benders（max 版本） | 244 | ≈1 | 1 轮（σ 全零收敛） | 0.0 | 是 | 上界=244 |
| 04 拉格朗日松弛 | 244 | ≈5 | 800 轮 | 0.0 | 是 | 对偶上界 244（首轮即达） |
| 05 LBBD（逻辑 Benders） | 244 | ≈30 | 91 轮 / 90 no-good | 0.0 | 是 | 主问题候选最大化 + guillotine 检验子问题 |
| 07 Branch-and-Price | 244 | ≈1 | 1 节点 / 0 分支 | 0.0 | 是 | 根 LP 整数 |

说明：单张板材结构使 LP(conv) = 最大模式价值 = 244 = IP（与 vrptw/crew 同现象），
多数方法根节点即证明；LBBD 的主问题无几何约束、需 90 个 no-good 排除不可切候选后收敛——展示逻辑割机制。

## 交付文件

- 01_direct.ipynb / 02_column_generation.ipynb / 03_benders.ipynb / 04_lagrangian.ipynb /
  05_lbbd.ipynb / 06_unified_report.ipynb / 07_branch_and_price.ipynb（全部已执行）
- theory/cg2_对偶推导总览.ipynb（已执行）
- results.json / README.md
- scripts/：cg2_core.py（六方法核心，含精确 guillotine 检验递归）、pool_norot.pkl（缓存完整池）、
  test_norot.py、insp_pool.py、build_cg2_nbs.py、build_cg2_extra.py 及早期草稿

## 结论

六种方法全部实现并**一致得到最优 244**；最优性由「完整池枚举 + 精确 guillotine 检验 + 池 IP」证明，
与文献一致。
