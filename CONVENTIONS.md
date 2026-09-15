# 列生成 / 分解算法求解套件 —— 共享规范 (CONVENTIONS)

> 本文件是所有家族求解子任务必须遵守的统一规范。每个家族子任务开工前必须先完整阅读本文件。

## 1. 运行环境（所有 Python 代码必须在 WSL 内运行）

- WSL 发行版：`Ubuntu-20.04`（默认发行版）
- Python 解释器：`/home/zkjqw/miniconda3/envs/jupyter-env/bin/python`（Python 3.10.20）
- OR-Tools 9.15.6755（pip 安装于 jupyter-env）
  - CP-SAT：`from ortools.sat.python import cp_model`
  - MathOpt：`from ortools.math_opt.python import mathopt`
    - 可用求解器（已实测）：`HIGHS`(MIP/LP)、`GLOP`(纯LP，变量必须连续)、`GSCIP`(MIP)、`CP_SAT`(MathOpt 接口)、`PDLP`
    - 勘误（9.15 实测）：math_opt 没有 `mathopt.inf` 常量，无穷大用 `float("inf")`；结果枚举用 `res.termination.reason.name` 取字符串
    - 勘误 2（scp_beasley 家族实测）：`SolveParameters(time_limit=120.0)` 会报 `AttributeError: 'float' object has no attribute 'seconds'`，必须用 `time_limit=datetime.timedelta(seconds=120)`
- 从 Windows 侧调用 WSL 的模板：
  ```
  wsl.exe -d Ubuntu-20.04 bash -lc '/home/zkjqw/miniconda3/envs/jupyter-env/bin/python /mnt/d/exactTest/.../script.py'
  ```
- 执行 notebook 并保留输出：
  ```
  wsl.exe -d Ubuntu-20.04 bash -lc '/home/zkjqw/miniconda3/envs/jupyter-env/bin/jupyter nbconvert --to notebook --execute --inplace /mnt/d/exactTest/column-generation-solvers/<family>/<nb>.ipynb'
  ```
- Windows 路径 `D:\exactTest` 在 WSL 内为 `/mnt/d/exactTest`。
- 数据目录：`D:\exactTest\column-generation-testcases\...`（下载好的用例集合，只读使用，不要修改）。

## 2. 交付目录结构（在 D:\exactTest\column-generation-solvers 下）

每族一个文件夹，命名如下：

```
column-generation-solvers/
├── CONVENTIONS.md            # 本文件（主控编写）
├── 00_overview.md            # 总览（主控在全部完成后编写）
├── csp_1d_waescher/
│   ├── README.md
│   ├── results.json
│   ├── 01_direct.ipynb
│   ├── 02_column_generation.ipynb
│   ├── 03_benders.ipynb
│   ├── 04_lagrangian.ipynb
│   └── 05_lbbd.ipynb
├── bpp_falkenauer/           # 同上结构
├── scp_beasley/              # 同上结构
├── crew_scheduling/          # 同上结构
├── vrptw_solomon25/          # 同上结构
└── cutting_2d_cgcut/         # 同上结构
```

## 3. 每个 notebook 的统一要求

1. 文件名与顺序固定：`01_direct.ipynb`（直接建模）→ `02_column_generation.ipynb`（列生成）→ `03_benders.ipynb`（Benders 分解）→ `04_lagrangian.ipynb`（拉格朗日松弛）→ `05_lbbd.ipynb`（LBBD）。
2. 所有 notebook 必须用 nbconvert 实际执行过并保留 cell 输出（含打印的结果、耗时、迭代信息）。
3. Markdown cell（中文）必须包含：
   - 问题定义（一句话 + 集合/参数）
   - 完整数学模型（LaTeX：目标、约束、变量域；分解方法要写出主问题/子问题模型）
   - 方法原理要点（3-5 条）
   - 实现要点（用到的 OR-Tools API、数据结构、定价/割生成细节）
   - 运行结果（与基准最优值对比、gap、时间、迭代次数）
   - 结论与适用性讨论（该方法对该问题是否自然、优缺点）
4. 代码 cell 顶部统一打印环境信息：
   ```python
   import platform, time, ortools
   print("python", platform.python_version(), "| ortools", ortools.__version__)
   ```
5. 每个求解方法必须设置合理的停机条件（时间上限 ≤ 120 s、迭代上限、gap 阈值），不允许无限循环；未证最优时必须如实说明。
6. 所有随机行为必须可复现（固定种子或确定性算法）。
7. 代码要解析原始数据文件（不要硬编码实例数值），解析代码写在 notebook 中。
8. 每个 notebook 的最后一段 Markdown 给出该方法的结论（1-2 句），并注明"基准最优值来源"（直接 MIP 求得的证明最优值 / 文献值 / 仅上界等）。

## 4. 五类方法的统一实现范式（按问题家族适用性落地）

### 4.1 直接建模（01_direct.ipynb）
- 首选 OR-Tools 完整模型：
  - 线性整数模型用 MathOpt + `HIGHS`（若 HIGHS 在个别实例上过慢，可换 `GSCIP`，并在文档中注明）。
  - 含序列/调度/几何约束的模型用 `cp_model`（CP-SAT）。
- 给出该实例的证明最优解（作为其他四个方法的基准）；若 120s 内无法证明最优，报告 best bound 与 gap，并标注。

### 4.2 列生成（02_column_generation.ipynb）
- 主问题（受限主问题 RMP）用 MathOpt + `GLOP` 求解 LP 松弛。
- 定价子问题：
  - 1D 切割/装箱：0/1 背包 → 动态规划定价（负 reduced cost 模式）。
  - 集合覆盖/机组排班（列池已知或可枚举）：池扫描定价（对全部候选列计算 reduced cost，取最负者；注明与逐列动态定价的等价性——池完整时二者等价）。
  - VRPTW：先枚举可行路径池（CP-SAT 或 DFS + 可行性检查），再池定价。
  - 2D 切割：枚举 guillotine 模式池，池定价。
- 迭代到无负 reduced cost 列（容差 1e-7），记录：迭代次数、LP 最优值、定价时间。
- 整数解恢复：用生成的列池 + MathOpt MIP（HIGHS）求解整数版本，或说明 rounding 方案；给出整数目标与 gap。

### 4.3 Benders 分解（03_benders.ipynb）
- 结构（覆盖/划分/切割型问题通用）：把"模式（列）是否启用"作为整数主问题变量 y∈{0,1}，把"连续用量 x"作为 LP 子问题：
  - 主问题：min Σc_p·y_p + θ，y 二元，含已生成的最优性割 θ ≥ …。
  - 子问题：固定 y 后 LP（需求覆盖/件数约束 + x ≤ M·y），由对偶解生成 Benders 最优性割。
- 迭代至主问题下界与子问题上界收敛（gap 阈值 1e-4 或迭代上限），记录割数量与迭代轨迹。
- 最终在最优割对应的模式集上解整数 MIP 得整数解。
- 若问题天然不适合 Benders，仍需给出一个能运行的 Benders 变体并在文档中说明"此处 Benders 与该问题的契合度低、性能一般"的原因（不得跳过）。

### 4.4 拉格朗日松弛（04_lagrangian.ipynb）
- 松弛耦合约束（需求/覆盖约束），得到可分解子问题；用次梯度法（或替代梯度/体积法）求解对偶，迭代 200~1000 次，记录下界（min 问题）/上界（max 问题）轨迹。
- 用每轮 Lagrange 解恢复原始可行解（贪心修复 + 局部改进），报告可行解目标与对偶 gap。
- 说明松弛的是哪组约束、子问题结构、步长公式。

### 4.5 LBBD（05_lbbd.ipynb）
- 主问题：CP-SAT（或 MIP）决策"组合/分配"变量（物品→箱、任务→crew、客户→车、件类型→数量等）。
- 子问题：CP-SAT 可行性检查（容量/时间窗/2D 装得下/排序可行）。
- 不可行时回传逻辑割（no-good / 组合割），加回主问题，迭代至主问题最优且子问题可行，或达到迭代上限。
- 说明割的具体形式（如"分配到箱 b 的物品集合 S 容量超限 ⇒ 至少移走一个"）。
- 若子问题退化为线性容量检查，也要实现并说明（LBBD 退化为带惰性约束的主问题，这是 1D 情形下的合理形式）。

## 5. 家族任务清单与数据

| 家族文件夹 | 实例文件（WSL 路径） | 问题 | 目标 |
|---|---|---|---|
| csp_1d_waescher | /mnt/d/exactTest/column-generation-testcases/bpplib_1d/3_Wäscher_CSP/Waescher_CSP/Waescher_TEST0005.txt | 1D 切割下料（CSP），格式：第一行 m=57，第二行 L=10000，随后 m 行 "长度 需求" | 最小化标准辊数 |
| bpp_falkenauer | /mnt/d/exactTest/column-generation-testcases/bin_packing_falkenauer/binpack1.txt 中的第 1 个实例 u120_00（解析：首行 P=20；实例块：名字行 / "150 120 48" / 120 行物品尺寸） | 1D 装箱（BPP） | 最小箱数（文件自报最优 48） |
| scp_beasley | /mnt/d/exactTest/column-generation-testcases/set_covering/scp41.txt（格式：首行 m n=200 1000；第二行 1000 个列成本；随后 200 行 "覆盖列数 列号..."） | 集合覆盖（SCP） | 最小成本（文献最优 429） |
| crew_scheduling | /mnt/d/exactTest/column-generation-testcases/crew_scheduling/csp50.txt（格式见 OR-Library cspinfo：首行 N 时间上限；N 行任务起止时间；随后转移弧 i,j,cost；需自行检查是否含 depot 0/N+1） | 机组排班（Beasley-Cao） | 用最少可行 crew 数 K（自行论证下界并说明）使总成本最小 |
| vrptw_solomon25 | /mnt/d/exactTest/column-generation-testcases/vrptw_solomon_25/c101.txt（Solomon 25 客户，容量 200，时间窗 0..1236，服务时间 90） | VRPTW | 最小车辆数、其次最小距离（BKS：3 车，距离 191.81，需自行复算验证） |
| cutting_2d_cgcut | /mnt/d/exactTest/column-generation-testcases/cutting_stock_2d_guillotine/cgcut1.txt（格式：m=7；板材 15×10；每件 "长 宽 最大数量 价值"） | 2D guillotine 切割（单张板材，件数上限） | 最大化价值（文献最优 244） |

- 数据解析代码必须写进 notebook；若解析发现格式与上述说明不符，以实际文件为准并注明。
- 实例运行结果记录到 `<family>/results.json`，字段：
  ```json
  {"instance": "...", "best_known": 123, "proved_optimal": true,
   "methods": {"direct": {"objective": 123, "time_s": 1.2, "gap": 0.0, "iterations": 0, "note": ""}, ...}}
  ```

## 6. 质量红线

1. 不得伪造输出：所有结果必须来自 notebook 实际运行；超时/未收敛要如实记录。
2. 每个 notebook 独立可运行（数据路径存在、环境一致）。
3. 公式必须与代码一致。
4. 完成前必须用 nbconvert 重新执行全部 5 个 notebook 并确认无异常（kernel 报错必须修复）。
5. README.md 必须含统一结果对比表：
   ```
   | 方法 | 目标值 | 时间(s) | 迭代/割数 | gap% | 是否证明最优 | 备注 |
   ```
6. 全部完成后向主控回报一段总结（含 results.json 内容与遇到的坑）。
