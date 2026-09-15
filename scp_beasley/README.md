# scp_beasley：集合覆盖 scp41（OR-Library / Beasley）

实例：`scp41`，m=200 行、n=1000 列，列成本 1..100。目标最小化选中列总成本，使每行至少被一列覆盖。

- 数据文件：`/mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt`
- 环境：WSL Ubuntu-20.04，Python 3.10.20，OR-Tools 9.15.6755
- 求解器：MathOpt 的 HIGHS（MIP/LP）与 GLOP（LP）
- API 注意：`math_opt` 无 `mathopt.inf`，无穷大用 `float("inf")`；本环境 `SolveParameters.time_limit` 必须传 `datetime.timedelta(seconds=...)`（传 `120.0` 会在 solve 阶段报 AttributeError），GLOP 只接受连续变量。
- 基准最优值：**429**，由 `01_direct.ipynb` 直接 MIP 自证最优，与 OR-Library 文献值一致。

## 统一结果对比

| 方法 | 目标值 | 时间(s) | 迭代/割数 | gap% | 是否证明最优 | 备注 |
|---|---|---|---|---|---|---|
| 01 直接 MIP (HIGHS) | 429 | 0.25 | — | 0.0 | 是 | 66 列；solve_time≈0.070s |
| 02 列生成 (GLOP+HIGHS) | 429 | 3.02 | 166 次定价迭代 | 0.0 | 是 | LP 下界=429；池扫描定价；整数修复 0.112s |
| 03 Benders (GLOP+HIGHS) | 429 | 0.98 | 3 条最优性割 | 0.0 | 是 | LP 松弛下界=429；整数修复 0.153s |
| 04 拉格朗日松弛 | 429 | 1.54 | 600 次次梯度 | 0.0 | 是 | LB=428.9885；贪心 UB=434；MIP 修复 429 |
| 05 LBBD (逻辑 Benders) | 429 | 0.20 | 2 轮 / 200 条逻辑割 | 0.0 | 是 | 覆盖逻辑割=原覆盖约束；退化为惰性约束 MIP |

说明：

- 时间(s) 为 notebook 实际墙钟时间（方法主循环 + 最终整数/修复步骤）；各 notebook 单独运行时略有波动，以 notebook 输出为准。
- gap% 相对基准最优 429 计算，均为 0.0%。
- 列生成的定价为“完整池扫描”：池=全部 1000 列，扫描全部候选列取最小 reduced cost，与动态定价在池完整时等价（见 notebook 02）。
- 拉格朗日方法报告的是最终可行解目标 429；其下界为 428.98847765，对偶间隙相对 429 为 2.7e-05（约 0.0027%）。

## 交付文件

- `01_direct.ipynb`
- `02_column_generation.ipynb`
- `03_benders.ipynb`
- `04_lagrangian.ipynb`
- `05_lbbd.ipynb`
- `results.json`
- `README.md`
- `scripts/`（调试与 notebook 构建脚本）

## 执行与验证方式

每个 notebook 均用 nbconvert 实际执行并保留输出：

```bash
wsl.exe -d Ubuntu-20.04 bash -lc '/home/zkjqw/miniconda3/envs/jupyter-env/bin/jupyter nbconvert --to notebook --execute --inplace /mnt/d/exactTest/column-generation-solvers/scp_beasley/<nb>.ipynb'
```

执行后所有 notebook 无 kernel 错误；关键输出（目标值、时间、迭代数、最优性）见各 notebook 与上方结果表。

## 结论

- SCP 是**列生成与拉格朗日松弛的天然主场**：行覆盖约束的对偶直接给出列定价 reduced cost；松弛覆盖约束后子问题按列完全可分离。
- 本实例 LP 松弛最优值恰为整数最优 429，因此列生成 LP 下界与整数修复上界相等，可证明最优。
- **Benders** 对纯覆盖问题没有天然连续 recourse，本套件实现为“二元选列主问题 + 连续覆盖 LP 子问题”的对偶割变体，3 轮即得到 LP 下界 429，属可运行但非最自然的分解。
- **LBBD** 在覆盖问题中的逻辑割 `Σ_{j 覆盖 r} y_j ≥ 1` 就是原覆盖约束本身，因此 LBBD 退化为带惰性约束的 MIP，2 轮收敛、实现简单，是覆盖类问题很自然的分解形式。
