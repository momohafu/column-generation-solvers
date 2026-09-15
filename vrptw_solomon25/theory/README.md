# theory/ —— 数学建模与对偶推导详解系列（全部已执行）

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
